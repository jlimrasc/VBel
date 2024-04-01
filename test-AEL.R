# source("compute-AEL-functions.R")

# -----------------------------
# Main
# -----------------------------
set.seed(1)
x    <- runif(30, min = -5, max = 5)
elip <- rnorm(30, mean = 0, sd = 1)
y    <- 0.75 - x + elip
lam0 <- matrix(c(0,0), nrow = 2)
th   <- matrix(c(0.8277, -1.0050), nrow = 2)
a <- 0.001
z    <- cbind(x, y)
h    <- function(z, th) {
    xi <- z[1]
    yi <- z[2]
    h_zith <- c(yi - th[1] - th[2] * xi)
    h_zith[2] <- xi * h_zith[1]
    matrix(h_zith, nrow = 2)
}

ansAEL          <- compute_AEL_Rcpp(th, h, lam0, a, z, returnH = TRUE)
ansAELRcppprez  <- compute_AEL_Rcpp(th, h, lam0, a, z, returnH = TRUE)
ansAELRcpp      <- compute_AEL_Rcpp(th, h, lam0, a, z, useR_forz = TRUE, returnH = TRUE)

any(ansAEL$log_AEL==ansAELRcppprez$log_AEL)
any(ansAEL$lambda==ansAELRcppprez$lambda)
any(ansAEL$h_arr==ansAELRcppprez$h_arr)
any(ansAEL$H==ansAELRcppprez$H)
any(ansAELRcpp$log_AEL==ansAELRcppprez$log_AEL)
any(ansAELRcpp$lambda==ansAELRcppprez$lambda)
any(ansAELRcpp$h_arr==ansAELRcppprez$h_arr)
any(ansAELRcpp$H==ansAELRcppprez$H)


# tic("R")
# compute_AEL_R(th, h, lam0, a, z)
# toc()
# tic("Rcpp")
# compute_AEL_Rcpp(th, h, lam0, a, z)
# toc()
# tic("Rcpp_prez")
# compute_AEL_Rcpp(th, h, lam0, a, z, TRUE)
# toc()
