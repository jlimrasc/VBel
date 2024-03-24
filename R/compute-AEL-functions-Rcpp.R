#' Computes Accurate Empirical Likelihood Interference for a data set using cpp
#' with RcppEigen
#' 
#' @description 
#' `compute_AEL_Rcpp` returns the result of running AEL with the values present in 
#' the arguments.
#' 
#' @param th        Vector or scalar theta
#' @param h         User-defined function, outputs array
#' @param lam0      Initial vector for lambda
#' @param a         Scalar constant
#' @param z         n-1 by d matrix
#' @param T         Number of iterations using Newton-Raphson for estimation of lambda (default: 500)
#' @param useR_forz Bool whether to calculate the function first in R (True) or call the function in C (False) (default: True)
#'
#' @return The AEL of the data set
#' @export
#'
#' @examples compute_AEL_Rcpp(matrix(c(0.8277, -1.0050), nrow = 2), function(z, th) {matrix(c(z[2] - th[1] - th[2] * z[1], z[1]*(z[2] - th[1] - th[2] * z[1])), nrow = 2)}, matrix(c(0,0), nrow = 2), 0.001, cbind(runif(30, min = -5, max = 5), 0.75 - runif(30, min = -5, max = 5) + rnorm(30, mean = 0, sd = 1)))

compute_AEL_Rcpp <- function(th, h, lam0, a, z, T, useR_forz) {
    
    # -----------------------------
    # Default values
    # -----------------------------
    if (missing(T)) {
        T <- 500
    }
    if (missing(useR_forz)){
        useR_forz <- TRUE
    }
    
    if (!useR_forz) {
        compute_AEL_Rcpp_inner(th, h, lam0, a, z, T)
    } else if (useR_forz) {
        n <- nrow(z) + 1
        
        for (i in 1:(n - 1)) {
            zi <- t(z[i,]) # Row of z as vertical vector
            h_zith <- h(zi,th)
            
            h_sum <- h_sum + h_zith # For h(zn,th)
            H_Zth <- rbind(H_Zth, t(h_zith)) # Build up H(Z,th)
        }
        
        h_znth <- -a / (n - 1) * h_sum
        H_Zth <- rbind(H_Zth, t(h_znth)) # Last row of H is h(zn,th)
        compute_AEL_Rcpp_inner_prez(th, H_Zth, lam0, a, z, T)
    } else {
        warning("Error: Incorrect input for useR_forz")
    }
}