#' Computes Full-Covariance Gaussian VB posterior in R
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
#' @param verbosity     Integer for how often to print updates on current iteration number (default:500)
#' @param returnAll     Bool whether to return result for every line of the last iteration (default:FALSE)
#'
#' @returns A list containing:  \enumerate{
#'              \item A vector mu_FC
#'              \item A matrix C_FC 
#'              \item An array mu_FC_arr 
#'              \item An array C_FC_arr. 
#'              } Access using those names. If returnAll is TRUE, also inludes gmu, Egmu, delmu, Edelmu, gC_t, EgC, delC
#'
#' @export
#'
#' @seealso [compute_GVA_Rcpp()] for mix of R and C++ computation
#' 
#' @examples
#' set.seed(1)
#' x    <- runif(30, min = -5, max = 5)
#' elip <- rnorm(30, mean = 0, sd = 1)
#' y    <- 0.75 - x + elip
#' lam0 <- matrix(c(0,0), nrow = 2)
#' th   <- matrix(c(0.8277, -1.0050), nrow = 2)
#' a    <- 0.00001
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
#' T2 <- 500
#' rho <- 0.9
#' 
#' # -----------------------------
#' # Main
#' # -----------------------------
#' options(digits = 20)
#' set.seed(1)
#' ansGVA <-compute_GVA_R(mu, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2)
#' 
compute_GVA_R <- function(mu, C, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2, verbosity, returnAll) {
    # Initialise values
    if (missing(T)) { T <- 10000 }
    if (missing(T2)) { T2 <- 500 }
    if (missing(verbosity)) { verbosity <- 500 }
    if (missing(returnAll)) { returnAll <- FALSE }
    
    p           <- nrow(C)
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

    xi          <- matrix(stats::rnorm(T*p),T,p)                   # I     - Draw xi
    
    for (i in 1:(T)) {
        th      <- mu_t + C_t %*% xi[i,]                    # II    - Set theta
        gmu     <- compute_nabmu_ELBO(delth_logpi, delthh, 
                                      th, h, lam0, z, 
                                      n, a, T2)             # III   - Compute g_{mu}^{t+1}
        Egmu    <- rho * Egmu + (1 - rho) * gmu^2           # IV    - Accumulate gradients
        delmu   <- sqrt(Edelmu + elip * rep(1,p)) / 
            sqrt(Egmu + elip * rep(1,p)) * gmu              # V     - Compute update
        mu_t    <- mu_t + delmu                             # VI    - Update mean
        Edelmu  <- rho * Edelmu + (1 - rho) * delmu^2       # VII   - Accumulate updates
        gC_t    <- compute_nabC_ELBO(gmu, xi[i,], C_t)      # VIII  - Compute g_C^{t+1}
        gC_t[upper.tri(gC_t)] <- 0                          #       - Set gC_t to lower triag matx
        EgC     <- rho * EgC + (1 - rho) * gC_t^2           # IX    - Accumulate gradients
        delC    <- sqrt(EdelC + elip * M) / 
            sqrt(EgC + elip * M) * gC_t                     # X     - Compute update
        C_t     <- C_t + delC                               # XI    - Update covariance Cholesky
        EgC     <- rho * EgC + (1 - rho) * delC^2           # XII   - Accumulate updates
        # Store
        mu_arr[,i+1]   <- mu_t
        C_arr[,,i+1]    <- C_t
        
        if (verbosity && i %% verbosity == 0) { cat("Iteration:", i, "\n") }
    }

    if (!returnAll) {
        return(list(
            "mu_FC"  = mu_t,
            "C_FC"   = C_t,
            "mu_arr" = mu_arr,
            "C_arr"  = C_arr
        ))
    } else {
        return(list(
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
        ))
    }
}

compute_nabmu_ELBO <- function(delth_logpi, delthh, theta, h, lam0, z, n, a, T2) { 
    res <- compute_AEL_R(theta, h, lam0, a, z, T2, returnH = TRUE) # list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth)
    lambda <- res$"lambda"
    h_arr <- res$"h_arr"
    hznth <- h_arr[,,n]

    # Calculate gradient LogAEL with respect to theta
    nabth_logAEL <- 0 # Vector
    for (i in 1:(n-1)) {
        nabth_logAEL <- nabth_logAEL - (1/(1 + t(lambda) %*% h_arr[,,i]) - (a/(n-1)) / (1 + t(lambda) %*% hznth))[1] * (t(delthh(t(z[i,]), theta)) %*% lambda)
    }
    nabmu_ELBO <- nabth_logAEL + delth_logpi(theta)
}


compute_nabC_ELBO <- function(gmu, xi, C_t) {
    nabC_ELBO <- gmu %*% t(xi) + diag(1/diag(C_t))
}