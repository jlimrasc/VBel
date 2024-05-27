#' @title Adjusted Empirical Likelihood
#' 
#' @description 
#' Function for evaluating the Adjusted Empirical Likelihood for a given 
#' data set, moment conditions and parameter values
#' 
#' @param th        Vector or scalar theta
#' @param h         User-defined function, outputs array
#' @param lam0      Initial vector for lambda
#' @param a         Scalar constant
#' @param z         n-1 by d matrix
#' @param T         Number of iterations using Newton-Raphson for estimation of lambda (default: 500)
#' @param useR_forz Bool whether to calculate the function first in R (True) or call the function in C (False) (default: True)
#' @param returnH   Whether to return calculated values of h, H matrix and lambda
#'
#' @return A numeric value for the Adjusted Empirical Likelihood function 
#' computed evaluated at a given theta value
#' 
#' @useDynLib VBel, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#' 
#' @author Wei Chang Yu, Jeremy Lim
#' @references Yu, W., & Bondell, H. D. (2023). Variational Bayes for Fast and 
#' Accurate Empirical Likelihood Inference. Journal of the American Statistical 
#' Association, 1–13. \doi{doi:10.1080/01621459.2023.2169701}
#' 
#' @examples
#' # Generate toy variables
#' set.seed(1)
#' x    <- runif(30, min = -5, max = 5)
#' elip <- rnorm(30, mean = 0, sd = 1)
#' y    <- 0.75 - x + elip
#' 
#' # Set initial values for AEL computation
#' lam0 <- matrix(c(0,0), nrow = 2)
#' th   <- matrix(c(0.8277, -1.0050), nrow = 2)
#' a    <- 0.00001
#' T    <- 10
#' 
#' # Define Dataset and h-function
#' z    <- cbind(x, y)
#' h    <- function(z, th) {
#'     xi <- z[1]
#'     yi <- z[2]
#'     h_zith <- c(yi - th[1] - th[2] * xi, xi*(yi - th[1] - th[2] * xi))
#'     matrix(h_zith, nrow = 2)
#' }
#' ansAELRcpp <- compute_AEL(th, h, lam0, a, z, T, useR_forz = TRUE)
compute_AEL <- function(th, h, lam0, a, z, T, useR_forz, returnH) {
    
    # -----------------------------
    # Default values
    # -----------------------------
    if (missing(T)) { T <- 500 }
    if (missing(useR_forz)){ useR_forz <- TRUE }
    if (missing(returnH)){ returnH <- FALSE }
    
    if (!useR_forz) {
        res <- compute_AEL_Rcpp_inner_wrap(th, h, lam0, a, z, T)
        
    } else if (useR_forz) {
        p <- ncol(z)
        n <- nrow(z) + 1
        h_sum <- 0
        H_Zth <- c()

        for (i in 1:(n - 1)) {
            zi <- matrix(z[i, ], nrow = p) # Row of z as vertical vector
            h_zith <- h(zi, th)
            
            h_sum <- h_sum + h_zith # For h(zn,th)
            H_Zth <- rbind(H_Zth, t(h_zith)) # Build up H(Z,th)
        }
        
        h_znth <- -a / (n - 1) * h_sum
        H_Zth <- rbind(H_Zth, t(h_znth)) # Last row of H is h(zn,th)
        res <- compute_AEL_Rcpp_inner_prez(th, H_Zth, lam0, a, z, T)
        
    } else {
        warning("Error: Incorrect input for useR_forz")
        return()
    }
    
    if (!returnH) {
        res$log_AEL
    } else if (!useR_forz) {
        return(list(
            "log_AEL" = res[[1]],
            "lambda" = res[[2]],
            "h_arr" = array(unlist(res[[3]]), dim = c(1, ncol(z), nrow(z) + 1)),
            "H" = res[[4]]
        ))
    } else {
        return(list(
            "log_AEL" = res[[1]],
            "lambda" = res[[2]],
            "h_arr" = array(t(H_Zth), dim = c(1, ncol(H_Zth), n)),
            "H" = H_Zth
        ))
    }
}