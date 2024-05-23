test_that("GVA outputs right length, 2x2", {
    # -----------------------------
    # Initialise Variables
    # -----------------------------
    # Generate toy variables
    set.seed(1)
    x    <- runif(30, min = -5, max = 5)
    elip <- rnorm(30, mean = 0, sd = 1)
    y    <- 0.75 - x + elip
    
    # Set initial values for AEL computation
    lam0 <- matrix(c(0,0), nrow = 2)
    a    <- 0.00001
    
    # Define Dataset and h-function
    z    <- cbind(x, y)
    h    <- function(z, th) {
        xi <- z[1]
        yi <- z[2]
        h_zith <- c(yi - th[1] - th[2] * xi, xi*(yi - th[1] - th[2] * xi))
        matrix(h_zith, nrow = 2)
    }
    
    # Define h-gradient function
    delthh    <- function(z, th) {
        xi <- z[1]
        matrix(c(-1, -xi, -xi, -xi^2), 2, 2)
    }
    
    # Set initial values for GVA computation
    n       <- 31 # Number of rows in z
    reslm   <- lm(y ~ x)
    mu      <- matrix(unname(reslm$coefficients),2,1)
    C_0     <- unname(t(chol(vcov(reslm))))
    rho     <- 0.9
    
    # Set other variables for GVA
    delth_logpi <- function(theta) {-0.0001 * mu}
    elip    <- 10^-5
    T       <- 2 # Number of iterations for GVA
    T2      <- 5 # Number of iterations for AEL
    p       <- 2
    
    # -----------------------------
    # Main
    # -----------------------------
    set.seed(1)
    ansGVARcppHalf <-compute_GVA(mu, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2, fullCpp = FALSE)
    set.seed(1)
    ansGVARcppPure <-compute_GVA(mu, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2, fullCpp = TRUE)

    # Testing for length
    # (floating point errors and different random number generation between 
    # R & C++ make it difficult to test for number similarities)
    expect_length(ansGVARcppPure$mu_FC, p)
    expect_length(ansGVARcppHalf$mu_FC, p)
    expect_length(ansGVARcppPure$mu_arr, (T+1)*p)
    expect_length(ansGVARcppHalf$mu_arr, (T+1)*p)
    expect_length(ansGVARcppPure$C_FC, p*p)
    expect_length(ansGVARcppHalf$C_FC, p*p)
    expect_length(ansGVARcppPure$C_arr, (T+1)*p*p)
    expect_length(ansGVARcppHalf$C_arr, (T+1)*p*p)

    set.seed(NULL) # Reset seed
})