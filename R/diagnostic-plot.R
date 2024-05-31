#' Check the convergence of a data set computed by `compute_GVA`
#' 
#' Plots mu and variance in a time series plot to check for convergence of the 
#' computed data (i.e. Full-Covariance Gaussian VB Empirical Likelihood 
#' Posterior)
#' 
#'
#' @param dataList Named list of data generated from \link{compute_GVA}
#' @param muList Array of indices of mu_arr to plot. (default:all)
#' @param cList Matrix of indices of variance to plot, 2xn matrix, each row is 
#' (col,row) of variance matrix
#'
#' @return Matrix of variance of C_FC
#' @export
#'
#' @examples
#' # Generate toy variables
#' seedNum <- 100
#' set.seed(seedNum)
#' n       <- 100
#' p       <- 10
#' lam0    <- matrix(0, nrow = p)
#' 
#' # Calculate z
#' mean    <- rep(1, p)
#' sigStar <- matrix(0.4, p, p) + diag(0.6, p)
#' z       <- mvtnorm::rmvnorm(n = n-1, mean = mean, sigma = sigStar)
#' 
#' # Calculate intermediate variables
#' zbar    <- 1/(n-1) * matrix(colSums(z), nrow = p)
#' sumVal  <- matrix(0, nrow = p, ncol = p)
#' for (i in 1:p) {
#' zi      <- matrix(z[i,], nrow = p)
#' sumVal  <- sumVal + (zi - zbar) %*% matrix(zi - zbar, ncol = p)
#' }
#' sigHat  <- 1/(n-2) * sumVal
#' 
#' # Initial values for GVA
#' mu_0    <- matrix(zbar, p, 1)
#' C_0     <- 1/sqrt(n) * t(chol(sigHat))
#' 
#' # Define h-function
#' h       <- function(zi, th) { matrix(zi - th, nrow = 10) }
#' 
#' # Define h-gradient function
#' delthh  <- function(z, th) { -diag(p) }
#' 
#' # Set other initial values
#' delth_logpi <- function(theta) {-0.0001 * theta}
#' elip    <- 10^-5
#' T       <- 5 # Number of iterations for GVA
#' T2      <- 5 # Number of iterations for AEL
#' rho     <- 0.9
#' a       <- 0.00001
#' 
#' ansGVA <-compute_GVA(mu_0, C_0, h, delthh, delth_logpi, z, lam0, rho, elip, 
#' a, T, T2, fullCpp = TRUE)
#' 
#' diagnostic_plot(ansGVA)
#' diagnostic_plot(ansGVA, muList = c(1,4))
#' diagnostic_plot(ansGVA, cList = matrix(c(1,1, 5,6, 3,3), ncol = 2))
diagnostic_plot <- function(dataList, muList, cList) {
    # Function to generate cList
    gen_cList <- function() {
        if (p <= 3) {
            cList <- matrix(c(1, 1, 1, p), ncol = 2)
        } else {
            cList <- matrix(c(1, p, floor(p / 2), 
                              1, floor(p / 2), floor(p / 2) + 1), 
                            ncol = 2)
        }
    }
    
    # Function to generate muList
    gen_muList <- function() {
        if (p <= 3) {
            muList <- 1:p
        } else {
            muList <- c(1, floor(p / 2), p)
        }
    }
    
    # Function to check dataList structure
    check_strutcture <- function(dataList) {
        return(c("mu_FC", "mu_arr", "C_arr") %in% names(dataList))
    }
    
    # Validate inputs
    if (!is.list(dataList)) {
        try(stop("List not provided", call. = FALSE))
        return()
    } else if (!all(check_strutcture(dataList))) {
        try(stop(
            "List structure incorrect. Requires \"mu_FC\", \"mu_arr\", \"C_FC\" and \"C_arr\"",
            call. = FALSE
        ))
        return()
    }
    
    # Initialise variables and presets
    p <- length(dataList$mu_FC)         # Resolution of data
    T <- length(dataList$mu_arr) / p    # Number of iterations
    
    # Generate default lists
    if (missing(muList)) {
        muList <- gen_muList()
    } else if (is.null(muList)) {
        warning("muList is empty. Using default list")
        muList <- gen_muList()
    }
    if (missing(cList)) {
        cList <- gen_cList()
    } else if (is.null(cList)) {
        warning("cList empty. Using default list")
        cList <- gen_cList()
    } else if (ncol(cList) != 2) {
        warning("cList wrong dimensions - requires 2 columns. Using default cList.")
        cList <- gen_cList()
    }
    
    # Calculate variances
    variance_arr <- array(0, dim = c(p, p, T))
    for (t in 1:(T)) {
        variance_arr[, , t] = dataList$C_arr[, , t] %*% t(dataList$C_arr[, , t])
    }
    
    # Plot mu
    for (i in muList) {
        stats::plot.ts(
            dataList$mu_arr[i, ],
            main = sprintf("muFC %d", i),
            xlab = "Iterations",
            ylab = sprintf("mu_arr[%d,]", i)
        )
    }
    
    # Plot C
    for (i in 1:nrow(cList)) {
        stats::plot.ts(
            variance_arr[cList[i, 1], cList[i, 2], ],
            main = sprintf("Variance %d,%d", cList[i, 1], cList[i, 2]),
            xlab = "Iterations",
            ylab = sprintf("Variance[%d,%d,]", cList[i, 1], cList[i, 2])
        )
    }
    
    return(variance_arr)
}