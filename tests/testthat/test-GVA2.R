test_that("GVA outputs right length, 10x10", {
    # -----------------------------
    # Initialise Variables
    # -----------------------------
    # Generate toy variables
    seedNum <- 100
    set.seed(seedNum)
    n       <- 100
    p       <- 10
    lam0    <- matrix(0, nrow = p)
    
    # Calculate z
    mean    <- rep(1, p)
    sigStar <- matrix(0.4, p, p) + diag(0.6, p)
    z       <- mvtnorm::rmvnorm(n = n-1, mean = mean, sigma = sigStar)
    
    # Calculate intermediate variables
    zbar    <- 1/(n-1) * matrix(colSums(z), nrow = p)
    sumVal  <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
        zi      <- matrix(z[i,], nrow = p)
        sumVal  <- sumVal + (zi - zbar) %*% matrix(zi - zbar, ncol = p)
    }
    sigHat  <- 1/(n-2) * sumVal
    
    # Initial values for GVA
    mu_0    <- matrix(zbar, p, 1)
    C_0     <- 1/sqrt(n) * t(chol(sigHat))
    
    # Define h-function
    h       <- function(zi, th) { matrix(zi - th, nrow = 10) }
    
    # Define h-gradient function
    delthh  <- function(z, th) { -diag(p) }
    
    # Set other initial values
    delth_logpi <- function(theta) {-0.0001 * theta}
    elip    <- 10^-5
    T       <- 5 # Number of iterations for GVA
    T2      <- 5 # Number of iterations for AEL
    rho     <- 0.9
    a       <- 0.00001
    
    # -----------------------------
    # Main
    # -----------------------------
    options(digits = 20)
    set.seed(seedNum)
    ansGVARcppHalf <-compute_GVA(mu_0, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2, fullCpp = FALSE)
    ansGVARcppPure <-compute_GVA(mu_0, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, a, T, T2, fullCpp = TRUE)
    
    # Testing for length
    # (floating point errors and different random number generation between 
    # R & C++ make it difficult to test for number similarities)
    expect_length(ansGVARcppPure$mu_FC, p)
    expect_length(ansGVARcppHalf$mu_FC, p)
    
    expect_length(ansGVARcppPure$mu_arr, (T+1)*p)
    expect_length(ansGVARcppHalf$mu_arr, (T+1)*p)
    
    expect_length(ansGVARcppPure$C_FC, p*p)
    expect_length(ansGVARcppHalf$C_FC, p*p)

    # C_FC is upper triangular
    expect_equal(ansGVARcppHalf$C_FC[upper.tri(ansGVARcppHalf$C_FC)],rep(0,p*(p-1)/2))
    expect_equal(ansGVARcppPure$C_FC[upper.tri(ansGVARcppPure$C_FC)],rep(0,p*(p-1)/2))

    
    expect_length(ansGVARcppPure$C_arr, (T+1)*p*p)
    expect_length(ansGVARcppHalf$C_arr, (T+1)*p*p)
    
    expect_length(ansGVARcppHalf, 4)
    expect_length(ansGVARcppPure, 4)
    
    set.seed(NULL)
})