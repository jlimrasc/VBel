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

ans <-compute_AEL_R(th, h, lam0, a, z)
ans

# tic("R")
# compute_AEL_R(th, h, lam0, a, z)
# toc()
# tic("Rcpp")
# compute_AEL_Rcpp(th, h, lam0, a, z)
# toc()
# tic("Rcpp_prez")
# compute_AEL_Rcpp(th, h, lam0, a, z, TRUE)
# toc()
