#' Computes Adjusted Empirical Likelihood Interference for a data set using R and C++
#' with RcppEigen
#' 
#' @description 
#' `compute_AEL_Rcpp` returns the result of running AEL with the values present in the arguments.
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
#' @return The AEL of the data set
#' @useDynLib VBel, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#'
#' @seealso [compute_AEL_R()] for purely R computation
#' 
#' @examples compute_AEL_Rcpp(matrix(c(0.8277, -1.0050), nrow = 2), function(z, th) {matrix(c(z[2] - th[1] - th[2] * z[1], z[1]*(z[2] - th[1] - th[2] * z[1])), nrow = 2)}, matrix(c(0,0), nrow = 2), 0.001, cbind(runif(30, min = -5, max = 5), 0.75 - runif(30, min = -5, max = 5) + rnorm(30, mean = 0, sd = 1)))

compute_AEL_Rcpp <- function(th, h, lam0, a, z, T, useR_forz, returnH) {
    
    # -----------------------------
    # Default values
    # -----------------------------
    if (missing(T)) { T <- 500 }
    if (missing(useR_forz)){ useR_forz <- TRUE }
    if (missing(returnH)){ returnH <- FALSE }
    
    if (!useR_forz) {
        res <- compute_AEL_Rcpp_inner(th, h, lam0, a, z, T)
    } else if (useR_forz) {
        n <- nrow(z) + 1
        h_sum <- 0
        H_Zth <- c()
        
        for (i in 1:(n - 1)) {
            zi <- t(z[i,]) # Row of z as vertical vector
            h_zith <- h(zi,th)
            
            h_sum <- h_sum + h_zith # For h(zn,th)
            H_Zth <- rbind(H_Zth, t(h_zith)) # Build up H(Z,th)
        }
        
        h_znth <- -a / (n - 1) * h_sum
        H_Zth <- rbind(H_Zth, t(h_znth)) # Last row of H is h(zn,th)
        res <- compute_AEL_Rcpp_inner_prez(th, H_Zth, lam0, a, z, T)
    } else {
        warning("Error: Incorrect input for useR_forz")
        return(NULL)
    }
    
    if (!returnH) {
        res$log_AEL
    } else if(!useR_forz) {
        return(list("log_AEL" = res[[1]], "lambda" = res[[2]], "h_arr" = array(unlist(res[[3]]),dim = c(1,ncol(z),nrow(z)+1)), "H" = res[[4]]))
    } else {
        return(list("log_AEL" = res[[1]], "lambda" = res[[2]], "h_arr" = array(t(H_Zth), dim = c(1, ncol(H_Zth), n)), "H" = H_Zth))
    }
}