# # # # Diagnostic check - R
# p <- 10
# plot.ts(ansGVA3$mu_arr[1,], main = "muFC1 R")
# plot.ts(ansGVA3$mu_arr[2,], main = "muFC2 R")
# variance_arr <- array(0,dim=c(p,p,T+1))
# for (t in 1:(T+1)) {
#     variance_arr[,,t] = ansGVA$C_arr[,,t] %*% t(ansGVA$C_arr[,,t])
# }
# plot.ts(variance_arr[1,1,], main = "variance 1,1 R")
# plot.ts(variance_arr[1,2,], main = "variance 1,2 R")
# plot.ts(variance_arr[2,2,], main = "variance 2,2 R")
# 
# # Diagnostic check - RcppHalf
# plot.ts(ansGVARcppHalf3$mu_arr[1,], main = "muFC1 RcppHalf")
# plot.ts(ansGVARcppHalf3$mu_arr[2,], main = "muFC2 RcppHalf")
# variance_arr <- array(0,dim=c(p,p,T+1))
# for (t in 1:(T+1)) {
#     variance_arr[,,t] = ansGVARcppHalf2$C_arr[,,t] %*% t(ansGVARcppHalf$C_arr[,,t])
# }
# plot.ts(variance_arr[1,1,], main = "variance 1,1 RcppHalf")
# plot.ts(variance_arr[1,2,], main = "variance 1,2 RcppHalf")
# plot.ts(variance_arr[2,2,], main = "variance 2,2 RcppHalf")

# Diagnostic check - RcppPure
# for (i in 1:(p)) {
#     plot.ts(ansGVARcppPure3$mu_arr[i,], main = sprintf("muFC%d RcppPure", i))
# }
# variance_arr <- array(0,dim=c(p,p,T+1))
# for (t in 1:(T+1)) {
#     variance_arr[,,t] = ansGVARcppPure3$C_arr[,,t] %*% t(ansGVARcppPure3$C_arr[,,t])
# }
# plot.ts(variance_arr[1,1,], main = "variance 1,1 RcppPure")
# plot.ts(variance_arr[1,10,], main = "variance 1,10 RcppPure")
# plot.ts(variance_arr[5,5,], main = "variance 5,5 RcppPure")


#' @title Plotting mu and variance to check for convergence
#'
#' @param dataList Named list of data generated from \link{compute_GVA}
#' @param muList Array of indices of mu_arr to plot. (default:all)
#' @param cList Matrix of indices of variance to plot, 2xn matrix, each row is (col,row) of variance matrix
#'
#' @return Matrix of variance of C_FC
#' @export
#'
#' @examples
#' \dontrun{
#' diagnostic_plot(ansGVA)
#' diagnostic_plot(ansGVA, muList = c(1,4))
#' diagnostic_plot(ansGVA, cList = matrix(c(1,1, 5,6, 3,3), ncol = 2))}
diagnostic_plot <- function(dataList, muList, cList) {
    # Function to generate cList
    gen_cList <- function() {
        if (p <= 3) { cList <- matrix(c(1,1,1,p), ncol = 2) }
        else { cList <- matrix(c(1,p,floor(p/2), 1,floor(p/2),floor(p/2)+1), ncol = 2) }
    }
    
    # Function to generate muList
    gen_muList <- function() {
        if (p <= 3) { muList <- 1:p }
        else { muList <- c(1,floor(p/2),p) }
    }
    
    # Function to check dataList structure
    check_strutcture <- function(dataList) {
        return(c("mu_FC","mu_arr","C_arr") %in% names(dataList))
    }

    # Validate inputs
    if (!is.list(dataList)) {
        try(stop("List not provided", call. = FALSE))
        return()
    } else if (!all(check_strutcture(dataList))) {
        try(stop("List structure incorrect. Requires \"mu_FC\", \"mu_arr\", \"C_FC\" and \"C_arr\"", call. = FALSE))
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
    variance_arr <- array(0,dim=c(p,p,T))
    for (t in 1:(T)) {
        variance_arr[,,t] = dataList$C_arr[,,t] %*% t(dataList$C_arr[,,t])
    }
    
    # Plot mu
    for (i in muList) {
        stats::plot.ts(dataList$mu_arr[i,],
                main = sprintf("muFC %d", i),
                xlab = "Iterations", ylab = sprintf("mu_arr[%d,]",i))
    }
    
    # Plot C
    for (i in 1:nrow(cList)) {
        stats::plot.ts(variance_arr[cList[i,1],cList[i,2],], 
                main = sprintf("Variance %d,%d",cList[i,1],cList[i,2]),
                xlab = "Iterations", ylab = sprintf("Variance[%d,%d,]",cList[i,1],cList[i,2]))
    }
    
    return(variance_arr)
}