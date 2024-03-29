#' Computes Accurate Empirical Likelihood Interference for a data set using R
#'
#' @description
#' `compute_AEL` returns the result of running AEL with the values present in
#' the arguments.
#'
#' @param th        Vector or scalar theta, k x 1 matrix
#' @param h         User-defined moment-condition function, outputs a k x 1 matrix containing the kth row of h. Function must take two arguments: zi and theta for h(zi,th)
#' @param lam0      Initial guess for lambda, k x 1 matrix
#' @param a         Scalar AEL constant, must be >0, small values recommended
#' @param z         Data matrix, n-1 x d matrix
#' @param T         Number of iterations using Newton-Raphson for estimation of lambda (default: 500)
#' @param returnH   Whether to return calculated values of h, H matrix and lambda
#'
#' @return The AEL of the data set
#' @export
#'
#' @examples
#' set.seed(1)
#' x    <- runif(30, min = -5, max = 5)
#' elip <- rnorm(30, mean = 0, sd = 1)
#' y    <- 0.75 - x + elip
#' lam0 <- matrix(c(0,0), nrow = 2)
#' th   <- matrix(c(0.8277, -1.0050), nrow = 2)
#' a <- 0.001
#' z    <- cbind(x, y)
#' h    <- function(z, th) {
#'     xi <- z[1]
#'     yi <- z[2]
#'     h_zith <- c(yi - th[1] - th[2] * xi)
#'     h_zith[2] <- xi * h_zith[1]
#'     matrix(h_zith, nrow = 2)
#' }
#' ans <- compute_AEL_R(th, h, lam0, a, z)

compute_AEL_R <- function(th, h, lam0, a, z, T, returnH) {
    # -----------------------------
    # Default values
    # -----------------------------
    if (missing(T)) { T <- 500 }
    if (missing(returnH)) { returnH <- 0 }
    # -----------------------------
    # Starting variables (h(zi,th), h(zn,th), H)
    # -----------------------------
    n <- nrow(z) + 1
    h_arr <- array(dim = c(dim(h(t(z[1, ]),th)),n)) # Initialise
    h_sum <- 0
    H_Zth <- c()
    
    for (i in 1:(n - 1)) {
        zi <- t(z[i, ]) # Row of z as vertical vector
        h_zith <- h(zi, th)
        
        h_arr[,,i] <-
            h_zith # Collect all h(zi,th) so no need to calculate again
        
        h_sum <- h_sum + h_zith # For h(zn,th)
        H_Zth <- rbind(H_Zth, t(h_zith)) # Build up H(Z,th)
    }
    
    h_znth <- -a / (n - 1) * h_sum
    H_Zth <- rbind(H_Zth, t(h_znth)) # Last row of H is h(zn,th)
    
    h_arr[,,n] <- h_znth
    
    # -----------------------------
    # Lambda Calculation
    # -----------------------------
    lambda <- compute_lambda(h_arr, H_Zth, lam0, a, T, n)

    # -----------------------------
    # AEL Calculation
    # -----------------------------
    temp_accu <- 0
    for (i in 1:n) {
        temp_accu <- temp_accu + log(1 + t(lambda) %*% h_arr[,,i])
    }
    log_AEL <- -(temp_accu + n * log(n))

    if (!returnH) {
        log_AEL
    } else {
        return(list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth))
    }
}

compute_lambda <- function(h_arr, H_Zth, lam0, a, T, n) {
    #' Computes lambda<sup>T</sup> for AEL for a specific set of user-inputted variables
    #'
    #' @description
    #' `compute_lambda` returns the result of running a modified Newton-Raphson to
    #' find the lambda variable needed in calculating AEL
    #'
    #' @details
    #' juicy deets
    #'
    #' @param h_arr  List of h(z1,th) to h(zn,th)
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
            wi_arr[i] <- as.vector((1 + t(tild_lam) %*% h_arr[,,i]) ^ -1)
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
                vi <- 2 * n - n ^ 2 / wi
            }
            
            # Calculate sum for this iteration of dF
            dF <- dF + vi * h_arr[,,i]
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
            if (wi_arr[i]^-1 >= 1 / n) {
                vi2 <- wi_arr[i]
            } else {
                vi2 <- n
            }
            D_arr[i] <- vi2 ^ 2
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
}  