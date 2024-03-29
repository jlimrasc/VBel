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

#'
#' @returns             List of mu_FC and C_FC. Access using those names.
#' @export
#'
#' @examples
#' # example code
#' 
compute_GVA_Rcpp <- function(mu, C, h, delthh, delth_logpi, z, rho, elip, lam0, a, T, T2) {
    # Initialise values
    if (missing(T)) { T <- 10000 }
    if (missing(T2)) { T2 <- 500 }
    
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
    
    xi          <- matrix(rnorm(T*p),T,p)                   # I     - Draw xi
    
    for (i in 1:(T)) {
        th      <- mu_t + C_t %*% xi[i,]                    # II    - Set theta
        gmu     <- compute_nabmu_ELBO2(delth_logpi, delthh, 
                                      th, n, h, lam0, 
                                      a, z, T2)             # III   - Compute g_{mu}^{t+1}
        output <- compute_GVA_Rcpp_inner(rho, elip, Egmu, Edelmu, EgC, EdelC, gmu, mu_t, C_t, xi, M, p, i-1)
        # Egmu    <- rho * Egmu + (1 - rho) * gmu^2           # IV    - Accumulate gradients
        # delmu   <- sqrt(Edelmu + elip * rep(1,p)) /
        #     sqrt(Egmu + elip * rep(1,p)) * gmu              # V     - Compute update
        # mu_t    <- mu_t + delmu                             # VI    - Update mean
        # Edelmu  <- rho * Edelmu + (1 - rho) * delmu^2       # VII   - Accumulate updates
        # gC_t    <- compute_nabC_ELBO2(gmu, xi[i,], C_t)      # VIII  - Compute g_C^{t+1}
        # gC_t[upper.tri(gC_t)] <- 0                          #       - Set gC_t to lower triag matx
        # EgC     <- rho * EgC + (1 - rho) * gC_t^2           # IX    - Accumulate gradients
        # delC    <- sqrt(EdelC + elip * M) /
        #     sqrt(EgC + elip * M) * gC_t                     # X     - Compute update
        # C_t     <- C_t + delC                               # XI    - Update covariance Cholesky
        # EgC     <- rho * EgC + (1 - rho) * delC^2           # XII   - Accumulate updates
        # if (any(output[[1]] != mu_t)) {
        #     warning(sprintf("mu_t deviated at line %i",i))
        #     browser()
        # }
        # if (any(output[[2]] != C_t)) {
        #     warning(sprintf("C_t deviated at line %i",i))
        #     browser()
        # }
        # if (any(output[[3]] != Egmu)) {
        #     warning(sprintf("Egmu deviated at line %i",i))
        #     browser()
        # }
        # if (any(output[[4]] != delmu)) {
        #     warning(sprintf("delmu deviated at line %i",i))
        #     browser()
        # }
        # if (any(output[[5]] != Edelmu)) {
        #     warning(sprintf("Edelmu deviated at line %i",i))
        #     browser()
        # }
        # if (any(output[[6]] != gC_t)) {
        #     warning(sprintf("gC_t deviated at line %i",i))
        #     browser()
        # }
        # if (any(output[[7]] != EgC)) {
        #     warning(sprintf("EgC deviated at line %i",i))
        #     browser()
        # }
        # if (any(output[[8]] != delC)) {
        #     warning(sprintf("delC deviated at line %i",i))
        #     browser()
        # }
        mu_t <- output[[1]]
        C_t <- output[[2]]
        Egmu <- output[[3]]
        delmu <- output[[3]]
        Edelmu <- output[[3]]
        gC_t <- output[[3]]
        EgC <- output[[3]]
        delC <- output[[3]]
        # print(mu_t)
        # print(C_t)
        # Store
        mu_arr[,i+1]   <- mu_t
        C_arr[,,i+1]    <- C_t
    }
    
    return(list("mu_FC" = mu_t, "C_FC" = C_t, "mu_FC_arr" = mu_arr, "C_FC_arr" = C_arr))
}

compute_nabmu_ELBO2 <- function(delth_logpi, delthh, theta, n, h, lam0, a, z, T2) { 
    res <- compute_AEL_R(theta, h, lam0, a, z, T2, 1) # list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth)
    lambda <- res$"lambda"
    h_arr <- res$"h_arr"
    hznth <- h_arr[,,n]
    
    # Calculate gradient LogAEL with respect to theta
    nabth_logAEL <- 0 # Matrix
    for (i in 1:(n-1)) {
        nabth_logAEL <- nabth_logAEL - (1/(1 + t(lambda) %*% h_arr[,,i]) - (a/(n-1)) / (1 + t(lambda) %*% hznth))[1] * (t(delthh(z[i], theta)) %*% lambda)
    }
    nabmu_ELBO <- nabth_logAEL + delth_logpi(theta)
}


compute_nabC_ELBO2 <- function(gmu, xi, C_t) {
    nabC_ELBO <- gmu %*% t(xi) + diag(1/diag(C_t))
}