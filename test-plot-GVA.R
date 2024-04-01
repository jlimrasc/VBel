# Diagnostic check
plot.ts(ansGVA$mu_FC_arr[1,])
plot.ts(ansGVA$mu_FC_arr[2,])
variance_arr <- array(0,dim=c(2,2,T+1))
for (t in 1:(T+1)) {
    variance_arr[,,t] = ansGVA$C_FC_arr[,,t] %*% t(ansGVA$C_FC_arr[,,t])
}
plot.ts(variance_arr[1,1,])
plot.ts(variance_arr[1,2,])
plot.ts(variance_arr[2,2,])
