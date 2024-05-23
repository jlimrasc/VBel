#' @title Full-Covariance Gaussian VB Empirical Likelihood Posterior
#' 
#' @description
#' Function for computing the Full-Covariance Gaussian VB Empirical Likelihood Posterior
#' 
#'
#' @param mu            Column vector, initial value of Gaussian VB mean
#' @param C             Lower Triangular Matrix, initial value of Gaussian VB Cholesky
#' @param h             User-defined moment-condition function, outputs a k x 1 matrix containing the kth row of h. Function must take two arguments: zi and theta for h(zi,th)
#' @param delthh        User defined function, outputs k x p Jacobian matrix of h(zi,th) with respect to theta
#' @param delth_logpi   User-defined function, outputs p x 1 matrix, derivative of log prior function
#' @param z             Data matrix, n-1 x d matrix
#' @param rho           Scalar (0 <~ 1) ADADELTA accumulation constant
#' @param elip          Scalar numeric stability constant. Should be specified with a small value
#' @param lam0          Initial guess for lambda, k x 1 matrix
#' @param a             Scalar AEL constant, must be >0, small values recommended
#' @param T             Number of iterations for GVA (default:10,000)
#' @param T2            Number of iterations for Log AEL (default:500)
#' @param fullCpp       Bool whether to calculate the main section in cpp (TRUE) or only partially (FALSE, doing all the AEL calculations in R before handing values to cpp) (default: TRUE)
#' @param verbosity     Integer for how often to print updates on current iteration number (default:500)
#'
#' @returns A list containing:  \enumerate{
#'              \item mu_FC: VB Posterior Mean at final iteration. A vector of 
#'              size p x 1
#'              \item C_FC: VB Posterior Variance-Covariance (Cholesky) at 
#'              final iteration. A lower-triangular matrix of size p x p
#'              \item mu_FC_arr: VB Posterior Mean for each iteration. A matrix 
#'              of size p x T+1
#'              \item C_FC_arr: VB Posterior Variance-Covariance (Cholesky) for 
#'              each iteration. An array of size p x p x T
#'              }
#' 
#' @export
#' 
#' @author Wei Chang Yu, Jeremy Lim
#' @references Yu, W., & Bondell, H. D. (2023). Variational Bayes for Fast and 
#' Accurate Empirical Likelihood Inference. Journal of the American Statistical 
#' Association, 1–13. \url{https://doi.org/10.1080/01621459.2023.2169701}
#' 
#' @examples
#' set.seed(1)
#' x    <- runif(30, min = -5, max = 5)
#' elip <- rnorm(30, mean = 0, sd = 1)
#' y    <- 0.75 - x + elip
#' lam0 <- matrix(c(0,0), nrow = 2)
#' th   <- matrix(c(0.8277, -1.0050), nrow = 2)
#' a <- 0.00001
#' z    <- cbind(x, y)
#' h    <- function(z, th) {
#'     xi <- z[1]
#'     yi <- z[2]
#'     h_zith <- c(yi - th[1] - th[2] * xi, xi*(yi - th[1] - th[2] * xi))
#'     matrix(h_zith, nrow = 2)
#' }
#' 
#' delthh    <- function(z, th) {
#'     xi <- z[1]
#'     matrix(c(-1, -xi, -xi, -xi^2), 2, 2)
#' }
#' 
#' n <- 31
#' reslm <- lm(y ~ x)
#' mu <- matrix(unname(reslm$coefficients),2,1)
#' C_0 <- unname(t(chol(vcov(reslm))))
#' 
#' delth_logpi <- function(theta) {-0.0001 * mu}
#' elip <- 10^-5
#' T <- 10
#' T2 <- 50
#' rho <- 0.9
#' 
#' # -----------------------------
#' # Main
#' # -----------------------------
#' options(digits = 20)
#' set.seed(1)
#' ansGVARcppHalf <-compute_GVA(mu, C_0, h, delthh, delth_logpi, z, lam0, 
#' rho, elip, a, T, T2, fullCpp = FALSE)
#' set.seed(1)
#' ansGVARcppPure <-compute_GVA(mu, C_0, h, delthh, delth_logpi, z, lam0, 
#' rho, elip, a, T, T2, fullCpp = TRUE)
#' 
compute_GVA <- function(mu, C, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2, fullCpp, verbosity) {
    # Initialise values
    if (missing(T)) { T <- 10000 }
    if (missing(T2)) { T2 <- 500 }
    if (missing(fullCpp)) { fullCpp <- TRUE }
    if (missing(verbosity)) { verbosity <- 500 }
    returnAll   <- FALSE

    p           <- nrow(C)
    xi          <- matrix(stats::rnorm(T*p),T,p)                # I     - Draw xi
    
    if (fullCpp) {
        res <- compute_GVA_Rcpp_inner_full(mu, C, h, delthh, delth_logpi, z, lam0, 
                                           rho, elip, a, T, T2, p, verbosity)
        res$mu_FC   <- matrix(res$mu_FC, nrow = p, ncol = 1)
        res$Egmu    <- matrix(res$Egmu, nrow = p, ncol = 1)
        res$delmu   <- matrix(res$delmu, nrow = p, ncol = 1)
        res$Edelmu  <- matrix(res$Edelmu, nrow = p, ncol = 1)
        res$gmu     <- matrix(res$gmu, nrow = p, ncol = 1)

        
        # Store
        mu_arr  <- res$mu_arr
        C_arr   <- array(unlist(res$C_arr), dim = c(dim(C),T+1))
        res$C_arr <- C_arr
        
    } else {
        Egmu        <- numeric(p)
        Edelmu      <- numeric(p)
        EgC         <- matrix(0, nrow = p, ncol = p)
        EdelC       <- matrix(0, nrow = p, ncol = p)
        mu_t        <- mu
        mu_arr      <- matrix(0,nrow = p, ncol = T+1)#array(dim = c(dim(mu_t), T+1))
        mu_arr[,1]  <- mu_t
        C_t         <- C        # Covariance Cholesky
        C_arr       <- array(dim = c(dim(C_t), T+1))
        C_arr[,,1]  <- C_t
        M           <- matrix(1,p,p)
        n           <- nrow(z) + 1
        
        for (i in 1:(T)) {
            th      <- mu_t + C_t %*% xi[i,]                    # II    - Set theta
            gmu     <- compute_nabmu_ELBO_RcppfromR(delth_logpi, delthh, 
                                                    th, h, lam0, z, 
                                                    n, a, T2)             # III   - Compute g_{mu}^{t+1}
            res <- compute_GVA_Rcpp_inner_IVtoXII(rho, elip, Egmu, Edelmu, EgC, EdelC, gmu, mu_t, C_t, xi, M, p, i-1)
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
            mu_arr[,i+1]   <- mu_t
            C_arr[,,i+1]    <- C_t
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

compute_nabmu_ELBO_RcppfromR <- function(delth_logpi, delthh, theta, h, lam0, z, n, a, T2) {
    res <- compute_AEL(theta, h, lam0, a, z, T2, returnH = TRUE) # list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth)
    lambda <- res$"lambda"
    h_arr <- res$"h_arr"
    hznth <- h_arr[,,n]
    
    # Calculate gradient LogAEL with respect to theta
    nabth_logAEL <- 0 # Matrix
    for (i in 1:(n-1)) {
        nabth_logAEL <- nabth_logAEL - (1/(1 + t(lambda) %*% h_arr[,,i]) - (a/(n-1)) / (1 + t(lambda) %*% hznth))[1] * (t(delthh(t(z[i,]), theta)) %*% lambda)
    }
    nabmu_ELBO <- nabth_logAEL + delth_logpi(theta)
}