#' Compute the Full-Covariance Gaussian VB Empirical Likelihood Posterior
#' 
#' Requires a given data set, moment conditions and parameter values and returns
#' a list of the final mean and variance-covariance along with an array of the 
#' in-between calculations at each iteration for analysis of convergence
#' 
#'
#' @param mu0           p x 1 initial vector of Gaussian VB mean
#' @param C0            p x p initial lower triangular matrix of Gaussian VB Cholesky
#' @param h             User-defined moment-condition function. 
                     #' Note that output should be an n x K matrix where K is necessarily \eqn{\geq}{>=} p. 
                     #' Input format for h should be (zi, th) were zi corresponds to the ith observation's data and th is the parameter vector
#' @param delthh        User-defined first-order derivative of moment-condition function. 
                     #' Note that output should be a K x p  matrix of h(zi,th) with respect to theta. 
                     #' Input format for h should be (zi, th) were zi corresponds to the ith observation's data and th is the parameter vector
#' @param delth_logpi   User-defined first-order derivative of log-prior function. 
                     #' Note that output should be a p x 1 vector. Input format for h should be (th) the parameter vector
#' @param z             Data matrix, n x d matrix
#' @param lam0          Initial vector for Lagrange multiplier lambda
#' @param rho           Scalar numeric beteen 0 to 1. ADADELTA accumulation constant
#' @param epsil         Positive numeric scalar stability constant. Should be specified with a small value
#' @param a             Positive scalar adjustment constant. For more accurate calculations, small values are recommended
#' @param SGD_iters     Number of Stochastic Gradient-Descent iterations for optimising mu and C. Default: 10,000
#' @param AEL_iters     Number of iterations using Newton-Raphson for optimising AEL lambda. Default: 500
#' @param verbosity     Integer for how often to print updates on current iteration number. To turn off, set to 0. Default:500
#'
#' @returns A list containing:  \enumerate{
#'              \item mu_FC: VB Posterior Mean at final iteration. A vector of 
#'              size p x 1
#'              \item C_FC: VB Posterior Variance-Covariance (Cholesky) at 
#'              final iteration. A lower-triangular matrix of size p x p
#'              \item mu_FC_arr: VB Posterior Mean for each iteration. A matrix 
#'              of size p x (SGD_iters + 1)
#'              \item C_FC_arr: VB Posterior Variance-Covariance (Cholesky) for 
#'              each iteration. An array of size p x p x (SGD_iters + 1)
#'              }
#' 
#' @export
#' 
#' @author Weichang Yu, Jeremy Lim
#' @references Yu, W., & Bondell, H. D. (2023). "Variational Bayes for Fast and 
#' Accurate Empirical Likelihood Inference", Journal of the American Statistical 
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
#' h    <- function(zi, th) {
#'     xi     <- zi[1]
#'     yi     <- zi[2]
#'     h_zith <- c(yi - th[1] - th[2] * xi, xi*(yi - th[1] - th[2] * xi))
#'     matrix(h_zith, nrow = 2)
#' }
#' 
#' delthh <- function(zi, th) {
#'     xi <- zi[1]
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
#' SGD_iters <- 50
#' epsil     <- 10^-5
#' rho       <- 0.9
#' 
#' # -----------------------------
#' # Main
#' # -----------------------------
#' result <-compute_GVA(mu0, C0, h, delthh, delth_logpi, z, lam0, 
#' rho, epsil, a, SGD_iters, AEL_iters)
#' 
compute_GVA <- function(mu0, C0, h, delthh, delth_logpi, z, lam0, rho, epsil, a, 
                        SGD_iters = 10000, AEL_iters = 500, verbosity = 500) {
    # Initialise values

    p           <- nrow(C0)
    
    res <- compute_GVA_Rcpp_inner_full(mu0, C0, h, delthh, delth_logpi, z, lam0, 
                                       rho, epsil, a, SGD_iters, AEL_iters, p, verbosity)
    
    # Return necessary values
    res2 <- list(
        "mu_FC"  = matrix(res$mu_FC, nrow = p, ncol = 1),
        "C_FC"   = res$C_FC,
        "mu_arr" = res$mu_arr,
        "C_arr"  = array(unlist(res$C_arr), dim = c(dim(C0), SGD_iters+1))
    )
    return(res2)
}