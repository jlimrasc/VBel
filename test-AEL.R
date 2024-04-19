# source("compute-AEL-functions.R")
# REQUIRES TIC-TOC PACKAGE TO MEASURE TIMES, COMMENT OUT NECESSARY LINES IF NOT USING THEM

# -----------------------------
# Initialise Variables
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

# -----------------------------
# Main
# -----------------------------
# tic("R")
ansAEL          <- compute_AEL_Rcpp(th, h, lam0, a, z)
# toc()
# tic("Rcpp")
ansAELRcppprez  <- compute_AEL_Rcpp(th, h, lam0, a, z, useR_forz = FALSE)
# toc()
# tic("Rcpp_prez")
ansAELRcpp      <- compute_AEL_Rcpp(th, h, lam0, a, z, useR_forz = TRUE)
# toc()

# Testing for discrepencies, FALSE if the same (might be floating point errors even with rounding)
round(ansAEL,5) == round(ansAELRcpp,5)
round(ansAEL,5) == round(ansAELRcppprez,5)

