# Diagnostic check - R
plot.ts(ansGVA$mu_arr[1,], main = "muFC1 R")
plot.ts(ansGVA$mu_arr[2,], main = "muFC2 R")
variance_arr <- array(0,dim=c(2,2,T+1))
for (t in 1:(T+1)) {
    variance_arr[,,t] = ansGVA$C_arr[,,t] %*% t(ansGVA$C_arr[,,t])
}
plot.ts(variance_arr[1,1,], main = "variance 1,1 R")
plot.ts(variance_arr[1,2,], main = "variance 1,2 R")
plot.ts(variance_arr[2,2,], main = "variance 2,2 R")

# Diagnostic check - RcppHalf
plot.ts(ansGVARcppHalf$mu_arr[1,], main = "muFC1 RcppHalf")
plot.ts(ansGVARcppHalf$mu_arr[2,], main = "muFC2 RcppHalf")
variance_arr <- array(0,dim=c(2,2,T+1))
for (t in 1:(T+1)) {
    variance_arr[,,t] = ansGVARcppHalf$C_arr[,,t] %*% t(ansGVARcppHalf$C_arr[,,t])
}
plot.ts(variance_arr[1,1,], main = "variance 1,1 RcppHalf")
plot.ts(variance_arr[1,2,], main = "variance 1,2 RcppHalf")
plot.ts(variance_arr[2,2,], main = "variance 2,2 RcppHalf")

# Diagnostic check - RcppPure
plot.ts(ansGVARcppPure$mu_arr[1,], main = "muFC1 RcppPure")
plot.ts(ansGVARcppPure$mu_arr[2,], main = "muFC2 RcppPure")
variance_arr <- array(0,dim=c(2,2,T+1))
for (t in 1:(T+1)) {
    variance_arr[,,t] = ansGVARcppPure$C_arr[,,t] %*% t(ansGVARcppPure$C_arr[,,t])
}
plot.ts(variance_arr[1,1,], main = "variance 1,1 RcppPure")
plot.ts(variance_arr[1,2,], main = "variance 1,2 RcppPure")
plot.ts(variance_arr[2,2,], main = "variance 2,2 RcppPure")