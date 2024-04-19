# source("R/compute-GVA-functions.R")
# source("R/compute-AEL-functions.R")
# REQUIRES TIC-TOC PACKAGE TO MEASURE TIMES, COMMENT OUT NECESSARY LINES IF NOT USING THEM

# -----------------------------
# Initialise Variables
# -----------------------------
set.seed(1)
x    <- runif(30, min = -5, max = 5)
elip <- stats::rnorm(30, mean = 0, sd = 1)
y    <- 0.75 - x + elip
lam0 <- matrix(c(0,0), nrow = 2)
th   <- matrix(c(0.8277, -1.0050), nrow = 2)
a    <- 0.00001
z    <- cbind(x, y)
h    <- function(z, th) {
    xi <- z[1]
    yi <- z[2]
    h_zith <- c(yi - th[1] - th[2] * xi, xi*(yi - th[1] - th[2] * xi))
    matrix(h_zith, nrow = 2)
}

delthh    <- function(z, th) {
    xi <- z[1]
    matrix(c(-1, -xi, -xi, -xi^2), 2, 2)
}

n <- 31
reslm <- lm(y ~ x)
mu <- matrix(unname(reslm$coefficients),2,1)
C_0 <- unname(t(chol(vcov(reslm))))

delth_logpi <- function(theta) {-0.0001 * mu}
elip <- 10^-5
T <- 10
T2 <- 500
rho <- 0.9

# -----------------------------
# Main
# -----------------------------
options(digits = 20)
set.seed(1)
# tic("R")
ansGVA <-compute_GVA_R(mu, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2)
# toc()
set.seed(1)
# tic("RcppHalf")
ansGVARcppHalf <-compute_GVA_Rcpp(mu, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2, fullCpp = FALSE)
# toc()
set.seed(1)
# tic("RcppPure")
ansGVARcppPure <-compute_GVA_Rcpp(mu, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2, fullCpp = TRUE)
# toc()

# Testing for discrepencies, FALSE if the same (might be floating point errors even with rounding)
any(round(ansGVA$mu_arr,5) != round(ansGVARcppPure$mu_arr,5))
any(round(ansGVA$C_arr,5) != round(ansGVARcppPure$C_arr,5))
any(round(ansGVA$mu_arr,5) != round(ansGVARcppHalf$mu_arr,5))
any(round(ansGVA$C_arr,5) != round(ansGVARcppHalf$C_arr,5))

# Results
ansGVA$mu_FC
ansGVARcppHalf$mu_FC
ansGVARcppPure$mu_FC
ansGVA$C_FC
ansGVARcppHalf$C_FC
ansGVARcppPure$C_FC