compute_lambda <- function(th, h, lam0, a, z, T, n) {
  #' Computes lambda T for AEL for a specific set of user-inputted variables
  #' 
  #' @description 
  #' `compute_lambda` returns the result of running a modified Newton-Raphson to 
  #' find the lambda variable needed in calculating AEL
  #' 
  #' @details 
  #' juicy deets
  #' 
  #' @param th    Vector or scalar
  #' @param h     User-defined function, outputs array
  #' @param lam0  Initial vector for lambda
  #' @param a     Scalar constant
  #' @param z     n-1 by d matrix
  #' @param T     Number of iterations using Newton-Raphson for estimation of lambda
  #' @param n     Height/Width of z + 1

  ## Functions  
  # -----------------------------
  # Compute dF
  # -----------------------------
  get_dF <- function(tild_lam) {
    #' @param tild_lam  Vector, current iteration of lambda
    dF <- 0
    wi_arr <- c() # Store for use in D in P
    vi_arr <- c() # Store for use in D in P
    for (i in 1:n) {
      wi <- as.vector((1 + t(tild_lam) %*% h_list[[i]]) ^-1)
      wi_arr[i] <- wi # Store in vector
      
      # Evaluate vi
      if (wi^-1 >= 1 / n) {
        vi <- wi
      } else {
        vi <- 2 * n - n^2 / wi
      }
      
      vi_arr[i] <- vi # Store in vector
      dF <- dF + vi * h_list[[i]]
    }
    #browser()
    return(list(dF = dF, wi_arr = wi_arr, vi_arr = vi_arr))
  }
  
  # -----------------------------
  # Compute d2F
  # -----------------------------
  get_d2F <- function(wi_arr) {
    #' @param wi_arr    Array of wi(th, lam~)
    
    # Build diagonal matrix
    # Diagonal matrix changes with respect to lambda's current guess so need 
    # to recompute each iteration
    D_arr <- c()
    for (i in 1:n) {
      if (wi_arr[i]^-1 >= 1/n) {
        vi2 <- wi_arr[i]
      } else {
        vi2 <- n#^2
      }
      D_arr[i] <- vi2^2
    }
    D <- diag(D_arr)
    
    # Find P i.e. d2F
    P <- -t(H_Zth) %*% D %*% H_Zth
    P
  }
  
  # -----------------------------
  # Starting variables (h(zi,th), h(zn,th), H)
  # -----------------------------
  h_list <- list() # Initialise
  h_sum <- 0
  H_Zth <- c()
  for (i in 1:(n - 1)) {
    zi <- t(z[i,]) # Row of z as vertical vector
    #browser()
    h_zith <- h(zi,th)
    #browser()
    h_list[[i]] <- h_zith # Collect all h(zi,th) so no need to calculate again
    #browser()
    h_sum <- h_sum + h_zith # For h(zn,th)
    #browser()
    H_Zth <- rbind(H_Zth, t(h_zith)) # Build up H(Z,th)
  }
  h_znth <- -a / (n - 1) * h_sum
  H_Zth <- rbind(H_Zth, t(h_znth)) # Last row of H is h(zn,th)
  h_list[[n]] <- h_znth
  #browser()
  # -----------------------------
  # Compute lambda using modified Newton-Raphson
  # -----------------------------
  lam_prev <- lam0
  for (i in 1:T) {
    # dF
    dF_res <- get_dF(lam_prev)
    dF <- dF_res$dF
    print(dF_res$wi_arr)
    print(dF_res$vi_arr)
    print("")
    
    # P
    P <- get_d2F(dF_res$wi_arr)
    lam_prev <- lam_prev - P^-1 %*% dF
    
  }
  lambdaT <- lam_prev
  
  return(list(lambda = lambdaT, h_store = h_list))
}  
  
compute_AEL <- function(th, h, lam0, a, z, T) {
  #' Computes __??__
  #' 
  #' @description 
  #' `compute_AEL` returns the result of running AEL with the values present in 
  #' the arguments.
  #' 
  #' @details 
  #' juicy deets
  #' 
  #' @param th    Vector or scalar
  #' @param h     User-defined function, outputs array
  #' @param lam0  Initial vector for lambda
  #' @param a     Scalar constant
  #' @param z     n-1 by d matrix
  #' @param T     Number of iterations using Newton-Raphson for estimation of lambda

  # -----------------------------
  # Default values
  # -----------------------------
  if (missing(T)) {
    T <- 100
  }
  
  # -----------------------------
  # Lambda Calculation
  # -----------------------------
  n <- nrow(z) + 1
  lam_res <- compute_lambda(th, h, lam0, a, z, T, n)
  lambda <- lam_res$lambda
  h_list <- lam_res$h_store
  
  # -----------------------------
  # AEL Calculation
  # -----------------------------
  temp_accu <- 0
  print(lambda)
  for (i in 1:n) {
    # print(c(i,t(lambda) %*% h_list[[i]]))
    temp <- log(1 + t(lambda) %*% h_list[[i]]) # FIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    temp_accu <- temp_accu + temp
  }
  log_AEL <- -(temp_accu + n * log(n))
}
