#' Computes Accurate Empirical Likelihood Interference for a data set using R
#' 
#' @description 
#' `compute_AEL` returns the result of running AEL with the values present in 
#' the arguments.
#' 
#' @param th    Vector or scalar
#' @param h     User-defined function, outputs array
#' @param lam0  Initial vector for lambda
#' @param a     Scalar constant
#' @param z     n-1 by d matrix
#' @param T     Number of iterations using Newton-Raphson for estimation of lambda (default: 100)
#'
#' @return The AEL of the data set
#' @export
#'
#' @examples compute_AEL_R(matrix(c(0.8277, -1.0050), nrow = 2), function(z, th) {matrix(c(z[2] - th[1] - th[2] * z[1], z[1]*(z[2] - th[1] - th[2] * z[1])), nrow = 2)}, matrix(c(0,0), nrow = 2), 0.001, cbind(runif(30, min = -5, max = 5), 0.75 - runif(30, min = -5, max = 5) + rnorm(30, mean = 0, sd = 1)))

compute_AEL_R <- function(th, h, lam0, a, z, T) {

  # -----------------------------
  # Default values
  # -----------------------------
  if (missing(T)) {
    T <- 100
  }
  
  # -----------------------------
  # Starting variables (h(zi,th), h(zn,th), H)
  # -----------------------------
  h_list <- list() # Initialise
  h_sum <- 0
  H_Zth <- c()
  
  n <- nrow(z) + 1
  
  for (i in 1:(n - 1)) {
    zi <- t(z[i,]) # Row of z as vertical vector
    h_zith <- h(zi,th)
    
    h_list[[i]] <- h_zith # Collect all h(zi,th) so no need to calculate again
    
    h_sum <- h_sum + h_zith # For h(zn,th)
    H_Zth <- rbind(H_Zth, t(h_zith)) # Build up H(Z,th)
  }
  
  h_znth <- -a / (n - 1) * h_sum
  H_Zth <- rbind(H_Zth, t(h_znth)) # Last row of H is h(zn,th)
  
  h_list[[n]] <- h_znth
  
  # -----------------------------
  # Lambda Calculation
  # -----------------------------
  lambda <- compute_lambda(h_list, H_Zth, lam0, a, T, n)

  # -----------------------------
  # AEL Calculation
  # -----------------------------
  temp_accu <- 0
  for (i in 1:n) {
    temp_accu <- temp_accu + log(1 + t(lambda) %*% h_list[[i]])
  }
  log_AEL <- -(temp_accu + n * log(n))
  
  log_AEL[1,1]
}

compute_lambda <- function(h_list, H_Zth, lam0, a, T, n) {
    #' Computes lambda<sup>T</sup> for AEL for a specific set of user-inputted variables
    #' 
    #' @description 
    #' `compute_lambda` returns the result of running a modified Newton-Raphson to 
    #' find the lambda variable needed in calculating AEL
    #' 
    #' @details 
    #' juicy deets
    #' 
    #' @param h_list  List of h(z1,th) to h(zn,th)
    #' @param H_Zth   Hessian matrix H
    #' @param lam0    Initial vector for lambda
    #' @param a       Scalar constant
    #' @param T       Number of iterations using Newton-Raphson for estimation of lambda
    #' @param n       Height/Width of z + 1
    #' 
    #' @keywords internal
    #' @noRd
    
    ## Functions  
    
    get_wi_arr <- function(tild_lam) {
        #' @param tild_lam  Vector, current iteration of lambda
        
        wi_arr <- c() # Store for use in D in P
        for (i in 1:n) {
            # store as vector so can use in multiplication later, if not it fails
            wi_arr[i] <- as.vector((1 + t(tild_lam) %*% h_list[[i]]) ^-1) 
        }
        
        wi_arr
    }
    # -----------------------------
    # Compute dF
    # -----------------------------
    get_dF <- function(wi_arr) {
        #' @param wi_arr Array of vectors of wi(th, lam~)
        dF <- 0
        for (i in 1:n) {
            wi <- wi_arr[i]# Store in vector
            
            # Evaluate vi
            if (wi^-1 >= 1 / n) {
                vi <- wi
            } else {
                vi <- 2 * n - n^2 / wi
            }
            
            # Calculate sum for this iteration of dF
            dF <- dF + vi * h_list[[i]]
        }
        
        dF
    }
    
    # -----------------------------
    # Compute d2F
    # -----------------------------
    get_d2F <- function(wi_arr) {
        #' @param wi_arr    Array of vectors of wi(th, lam~)
        
        # Build diagonal matrix
        # Diagonal matrix changes with respect to lambda's current guess so need 
        # to recompute each iteration
        D_arr <- c()
        for (i in 1:n) {
            if (wi_arr[i]^-1 >= 1/n) {
                vi2 <- wi_arr[i]
            } else {
                vi2 <- n
            }
            D_arr[i] <- vi2^2
        }
        D <- diag(D_arr)
        
        # Find P i.e. d2F
        P <- -t(H_Zth) %*% D %*% H_Zth
        
        P
    }
    
    # -----------------------------
    # Compute lambda using modified Newton-Raphson
    # -----------------------------
    lam_prev <- lam0 # Initial guess
    for (i in 1:T) {
        # wi
        wi_arr <- get_wi_arr(lam_prev)
        
        # dF
        dF <- get_dF(wi_arr)
        
        # P
        P <- get_d2F(wi_arr)
        lam_prev <- lam_prev - solve(P) %*% dF
    }
    lambdaT <- lam_prev # Final lambda
    
    lambdaT
}  