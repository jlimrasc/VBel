#' Compute the Full-Covariance Gaussian VB Empirical Likelihood Posterior
#' 
#' Requires a given data set, moment conditions and parameter values and returns
#' a list of the final mean and variance-covariance along with an array of the 
#' in-between calculations at each iteration for analysis of convergence
#' 
#'
#' @param mu0           p x 1 initial vector of Gaussian VB mean
#' @param C0            p x p initial lower triangular matrix of Gaussian VB Cholesky
#' @param h             User-defined moment-condition function. Note that output should be an (n-1) x K matrix where K is necessarily \eqn{\geq}{<=} p
#' @param delthh        User-defined first-order derivative of moment-condition function. Note that output should be a K x p  matrix of h(zi,th) with respect to theta
#' @param delth_logpi   User-defined first-order derivative of log-prior function. Note that output should be a p x 1 vector
#' @param z             Data matrix, n-1 x d matrix
#' @param lam0          Initial vector for Lagrange multiplier lambda
#' @param rho           Scalar numeric beteen 0 to 1. ADADELTA accumulation constant
#' @param epsil         Positive numeric scalar stability constant. Should be specified with a small value
#' @param a             Positive scalar adjustment constant. For more accurate calculations, small values are recommended
#' @param SDG_iters     Number of Stochastic Gradient-Descent iterations for optimising mu and C. Default: 10,000
#' @param AEL_iters     Number of iterations using Newton-Raphson for optimising AEL lambda. Default: 500
#' @param fullCpp       Bool whether to calculate the main section in cpp (TRUE) or only partially (FALSE, doing all the AEL calculations in R before handing values to cpp). Default: TRUE
#' @param verbosity     Integer for how often to print updates on current iteration number. Default:500
#'
#' @returns A list containing:  \enumerate{
#'              \item mu_FC: VB Posterior Mean at final iteration. A vector of 
#'              size p x 1
#'              \item C_FC: VB Posterior Variance-Covariance (Cholesky) at 
#'              final iteration. A lower-triangular matrix of size p x p
#'              \item mu_FC_arr: VB Posterior Mean for each iteration. A matrix 
#'              of size p x (SDG_iters + 1)
#'              \item C_FC_arr: VB Posterior Variance-Covariance (Cholesky) for 
#'              each iteration. An array of size p x p x (SDG_iters + 1)
#'              }
#' 
#' @export
#' 
#' @author Weichang Yu, Jeremy Lim
#' @references Yu, W., & Bondell, H. D. (2023). Variational Bayes for Fast and 
#' Accurate Empirical Likelihood Inference. Journal of the American Statistical 
#' Association, 1–13. \doi{doi:10.1080/01621459.2023.2169701}
#' 
#' @examples
#' # -----------------------------
#' # Initialise Inputs
#' # -----------------------------
#' # Generating 30 data points from a simple linear-regression model
#' set.seed(1)
#' x    <- runif(30, min = -5, max = 5)
#' vari <- rnorm(30, mean = 0, sd = 1)
#' y    <- 0.75 - x + vari
#' lam0 <- matrix(c(0,0), nrow = 2)
#' z    <- cbind(x, y)
#' 
#' # Specify moment condition functions for linear regression and its corresponding derivative
#' h    <- function(z, th) {
#'     xi     <- z[1]
#'     yi     <- z[2]
#'     h_zith <- c(yi - th[1] - th[2] * xi, xi*(yi - th[1] - th[2] * xi))
#'     matrix(h_zith, nrow = 2)
#' }
#' 
#' delthh <- function(z, th) {
#'     xi <- z[1]
#'     matrix(c(-1, -xi, -xi, -xi^2), 2, 2)
#' }
#' 
#' # Specify derivative of log prior
#' delth_logpi <- function(theta) { -0.0001 * mu0 }
#' 
#' # Specify AEL constant and Newton-Rhapson iteration
#' a         <- 0.00001
#' AEL_iters <- 10
#' 
#' # Specify initial values for GVA mean vector and Cholesky
#' reslm <- lm(y ~ x)
#' mu0   <- matrix(unname(reslm$coefficients),2,1)
#' C0    <- unname(t(chol(vcov(reslm))))
#' 
#' # Specify details for ADADELTA (Stochastic Gradient-Descent)
#' SDG_iters <- 50
#' epsil     <- 10^-5
#' rho       <- 0.9
#' 
#' # -----------------------------
#' # Main
#' # -----------------------------
#' resultHalfR <-compute_GVA(mu0, C0, h, delthh, delth_logpi, z, lam0, 
#' rho, epsil, a, SDG_iters, AEL_iters)
#' resultPureC <-compute_GVA(mu0, C0, h, delthh, delth_logpi, z, lam0, 
#' rho, epsil, a, SDG_iters, AEL_iters, fullCpp = TRUE)
#' 
compute_GVA <- function(mu0, C0, h, delthh, delth_logpi, z, lam0, rho, epsil, a, 
                        SDG_iters = 10000, AEL_iters = 500, fullCpp = TRUE, 
                        verbosity = 500) {
    # Initialise values
    returnAll   <- FALSE

    p           <- nrow(C0)
    
    if (fullCpp) {
        res <- compute_GVA_Rcpp_inner_full(mu0, C0, h, delthh, delth_logpi, z, lam0, 
                                           rho, epsil, a, SDG_iters, AEL_iters, p, verbosity)
        res$mu_FC   <- matrix(res$mu_FC, nrow = p, ncol = 1)
        res$Egmu    <- matrix(res$Egmu, nrow = p, ncol = 1)
        res$delmu   <- matrix(res$delmu, nrow = p, ncol = 1)
        res$Edelmu  <- matrix(res$Edelmu, nrow = p, ncol = 1)
        res$gmu     <- matrix(res$gmu, nrow = p, ncol = 1)

        
        # Store
        mu_arr  <- res$mu_arr
        C_arr   <- array(unlist(res$C_arr), dim = c(dim(C0), SDG_iters+1))
        res$C_arr <- C_arr

    } else {
        Egmu        <- numeric(p)
        Edelmu      <- numeric(p)
        EgC         <- matrix(0, nrow = p, ncol = p)
        EdelC       <- matrix(0, nrow = p, ncol = p)
        mu_t        <- mu0
        mu_arr      <- matrix(0, nrow = p, ncol = SDG_iters + 1)#array(dim = c(dim(mu_t), SDG_iters+1))
        mu_arr[,1]  <- mu_t
        C_t         <- C0        # Covariance Cholesky
        C_arr       <- array(dim = c(dim(C_t), SDG_iters + 1))
        C_arr[,,1]  <- C_t
        M           <- matrix(1, p, p)
        n           <- nrow(z) + 1
        xi          <- matrix(stats::rnorm(SDG_iters * p), SDG_iters, p)        # I     - Draw xi
        
        for (i in 1:(SDG_iters)) {
            th      <- mu_t + C_t %*% xi[i,]                                    # II    - Set theta
            gmu     <- compute_nabmu_ELBO_RcppfromR(delth_logpi, delthh, 
                                                    th, h, lam0, z, 
                                                    n, a, AEL_iters)            # III   - Compute g_{mu}^{t+1}
            res <- compute_GVA_Rcpp_inner_IVtoXII(rho, epsil, Egmu, Edelmu, EgC, 
                                                  EdelC, gmu, mu_t, C_t, xi, M, 
                                                  p, i - 1)                     # IV-XII
            
            # Overwrite with new values
            mu_t    <- res[[1]]
            C_t     <- res[[2]]
            Egmu    <- res[[3]]
            delmu   <- res[[4]]
            Edelmu  <- res[[5]]
            gC_t    <- res[[6]]
            EgC     <- res[[7]]
            delC    <- res[[8]]
            gmu     <- res[[9]]
            # Store
            mu_arr[,i + 1]   <- mu_t
            C_arr[,,i + 1]   <- C_t
            if (verbosity && i %% verbosity == 0) { cat("Iteration:", i, "\n") }
        }
        res <- list(
            "mu_FC"  = mu_t,
            "C_FC"   = C_t,
            "mu_arr" = mu_arr,
            "C_arr"  = C_arr,
            "gmu"    = gmu,
            "Egmu"   = Egmu,
            "delmu"  = delmu,
            "Edelmu" = Edelmu,
            "gC_t"   = gC_t, 
            "EgC"    = EgC, 
            "delC"   = delC
        )
    }
    
    if (!returnAll) {
        res2 <- list(
            "mu_FC"  = res$mu_FC,
            "C_FC"   = res$C_FC,
            "mu_arr" = res$mu_arr,
            "C_arr"  = res$C_arr
        )
        return(res2)
    } else {
        return(res)
    }
}

compute_nabmu_ELBO_RcppfromR <- function(delth_logpi, delthh, theta, h, lam0, z, n, a, AEL_iters) {
    res <- compute_AEL(theta, h, lam0, a, z, AEL_iters, returnH = TRUE) 
    # Returns a list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth)
    lambda <- res$"lambda"
    h_arr <- res$"h_arr"
    hznth <- h_arr[, , n]
    
    # Calculate gradient LogAEL with respect to theta
    nabth_logAEL <- 0 # Matrix
    for (i in 1:(n-1)) {
        nabth_logAEL <- nabth_logAEL - (1/(1 + t(lambda) %*% h_arr[,,i]) - (a/(n-1)) / (1 + t(lambda) %*% hznth))[1] * (t(delthh(t(z[i,]), theta)) %*% lambda)
    }
    nabmu_ELBO <- nabth_logAEL + delth_logpi(theta)
}