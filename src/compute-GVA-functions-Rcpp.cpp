// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <stdio.h>
#include <iostream>
#include <cmath>
#include "compute-AEL-functions-Rcpp.h"
// using Eigen::MatrixXd;      // variable size martix, double precision
// using Eigen::VectorXd;      // variable size martix, double precision
// using Eigen::VectorXi;      // variable size martix, integer precision
// using Eigen::DiagonalMatrix;// diagonal matrix
// using namespace std;
using namespace Rcpp;
using namespace RcppEigen;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// Functions
Eigen::MatrixXd compute_nabC_ELBO_Rcpp(Eigen::VectorXd gmu, Eigen::VectorXd xi, Eigen::MatrixXd C_t) {
    // Store wi for use in D in P
    Eigen::MatrixXd nabC_ELBO;
    nabC_ELBO = gmu * xi.transpose() + C_t.diagonal().cwiseInverse().asDiagonal().toDenseMatrix(); // .row changes it to a row vector rather than column but VectorXd makes it a column again so transpose
    return nabC_ELBO;
}
// compute_nabC_ELBO <- function(gmu, xi, C_t) {
//     nabC_ELBO <- gmu %*% t(xi) + diag(1/diag(C_t))
// }

Eigen::VectorXd compute_nabmu_ELBO_Rcpp(Rcpp::Function delth_logpi, Rcpp::Function delthh, 
                                        Eigen::VectorXd theta, Rcpp::Function h, 
                                        Eigen::VectorXd lam0, Eigen::MatrixXd z, 
                                        int n, double a, int p, int T2) {
    // Rcpp::Rcout << "before";
    Rcpp::List res = compute_AEL_Rcpp_inner(theta, h, lam0, a, z, T2); // list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth)
    // Rcpp::Rcout << " after";
    Eigen::MatrixXd lambda              = res["lambda"];
    std::vector<Eigen::VectorXd> h_arr  = res["h_arr"];
    Eigen::VectorXd hznth               = h_arr[n-1];

    // Calculate gradient LogAEL with respect to theta
    // Rcpp::Rcout << " calc";
    Eigen::VectorXd nabth_logAEL = Eigen::VectorXd::Zero(p); // Vector
    Eigen::MatrixXd delthh_zith;
    for (int i = 0; i < n-1; i++) {
        nabth_logAEL = nabth_logAEL - ((1/(1 + (lambda.transpose() * h_arr[i]).array()) - (a/(n-1)) / (1 + (lambda.transpose() * hznth).array()))(0,0) * (Rcpp::as<Eigen::MatrixXd>(delthh(z.row(i).transpose(),theta)).transpose() * lambda).array()).matrix();
        //Rcpp::Rcout << "\nlamT" << lambda.rows() << " " << lambda.cols() << " " << lambda << "\nharr" << h_arr[i].rows() << " " << h_arr[i].cols() << " " << h_arr[i] << "\nhznth" << hznth.rows() << " " << hznth.cols() << " " << hznth 
        // Rcpp::Rcout << "\nTerm1\n" << ((1/(1 + (lambda.transpose() * h_arr[i]).array()) - (a/(n-1)) / (1 + (lambda.transpose() * hznth).array()))(0,0) * (Rcpp::as<Eigen::MatrixXd>(delthh(z.row(i).transpose(),theta)).transpose() * lambda).array()).matrix() << "\nTerm2\n" << nabth_logAEL << "\n";
    }

    // Rcpp::Rcout << "\n" << nabth_logAEL << std::endl;
    // Rcpp::Rcout << " ELBO";
    Eigen::VectorXd nabmu_ELBO = (nabth_logAEL.array() + Rcpp::as<Eigen::VectorXd>(delth_logpi(theta)).array()).matrix();
    // Rcpp::Rcout.precision(20);
    // Rcpp::Rcout << nabmu_ELBO << std::fixed << " Done" << std::endl;
    return nabmu_ELBO;
}
// compute_nabmu_ELBO <- function(delth_logpi, delthh, theta, h, lam0, z, n, a, T2) { 
//     res <- compute_AEL_R(theta, h, lam0, a, z, T2, 1) # list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth)
//     lambda <- res$"lambda"
//     h_arr <- res$"h_arr"
//     hznth <- h_arr[,,n]
//     
// # Calculate gradient LogAEL with respect to theta
//     nabth_logAEL <- 0 # Matrix
//     for (i in 1:(n-1)) {
//         nabth_logAEL <- nabth_logAEL - (1/(1 + t(lambda) %*% h_arr[,,i]) - (a/(n-1)) / (1 + t(lambda) %*% hznth))[1] * (t(delthh(z[i], theta)) %*% lambda)
//     }
//     nabmu_ELBO <- nabth_logAEL + delth_logpi(theta)
// }


// [[Rcpp::export]]
std::vector<Eigen::MatrixXd> compute_GVA_Rcpp_inner_IVtoXII(const double rho, const double elip, Eigen::VectorXd Egmu, Eigen::VectorXd Edelmu, 
                                       Eigen::MatrixXd EgC, Eigen::MatrixXd EdelC, Eigen::VectorXd gmu, 
                                       Eigen::VectorXd mu_t, Eigen::MatrixXd C_t, 
                                       const Eigen::MatrixXd xi, Eigen::MatrixXd M, const int p, const int i) {
    
    Egmu = rho * Egmu + (1 - rho) * gmu.cwiseProduct(gmu);      // IV   - Accumulate gradients
    Eigen::VectorXd delmu = (((Edelmu + elip * Eigen::VectorXd::Ones(p)).cwiseSqrt().array() / (Egmu + elip * Eigen::VectorXd::Ones(p)).cwiseSqrt().array()) * gmu.array()).matrix(); // V     - Compute update
    mu_t = mu_t + delmu;                                        // VI    - Update mean
    Edelmu = rho * Edelmu + (1 - rho) * delmu.cwiseProduct(delmu);      // VII   - Accumulate updates
    Eigen::MatrixXd gC_t = (compute_nabC_ELBO_Rcpp(gmu, xi.row(i), C_t)).triangularView<Eigen::Lower>();          // VIII  - Compute g_C^{t+1}
    EgC = rho * EgC + (1 - rho) * gC_t.cwiseProduct(gC_t);           // IX    - Accumulate gradients
    Eigen::MatrixXd delC    = ((EdelC + elip * M).cwiseSqrt().array() / (EgC + elip * M).cwiseSqrt().array() * gC_t.array()).matrix();     // X     - Compute update
    C_t     = C_t + delC;                            // XI    - Update covariance Cholesky
    EgC     = rho * EgC + (1 - rho) * delC.cwiseProduct(delC);        // XII   - Accumulate updates
        
    // return Rcpp::List::create(
    //         Rcpp::Named("mu_t") = mu_t,
    //         Rcpp::Named("C_t") = C_t);
    // cout << gC_t << endl;
    std::vector<Eigen::MatrixXd> res = {mu_t, C_t, Egmu, delmu, Edelmu, gC_t, EgC, delC, gmu};
    // std::vector<int> v = {(gmu * xi.row(i)).rows(), (gmu * xi.row(i)).cols(), gC_t.rows(), gC_t.cols(), gmu.rows(), gmu.cols(), xi.row(i).rows(), xi.row(i).cols()};
    return res;
}
    
// [[Rcpp::export]]
Rcpp::List compute_GVA_Rcpp_inner_full(
        Eigen::VectorXd mu, Eigen::MatrixXd C, Rcpp::Function h, Rcpp::Function delthh,
        Rcpp::Function delth_logpi, Eigen::MatrixXd z, Eigen::VectorXd lam0, Eigen::MatrixXd xi, 
        double rho, double elip, double a, int T, int T2, int p, bool rFuncs){

    Eigen::VectorXd Egmu    = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd gmu     = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd Edelmu  = Eigen::VectorXd::Zero(p);
    Eigen::MatrixXd EgC     = Eigen::MatrixXd::Zero(p,p);
    Eigen::MatrixXd EdelC   = Eigen::MatrixXd::Zero(p,p);
    Eigen::VectorXd mu_t    = mu;
    Eigen::MatrixXd mu_arr  = Eigen::MatrixXd::Zero(p,T+1);
    mu_arr.col(0)           = mu_t; // Can save row vector to column?
    Eigen::MatrixXd C_t     = C;        // Covariance Cholesky
    std::vector<Eigen::MatrixXd> C_arr = {C_t}; // No preallocation?
    Eigen::MatrixXd M       = Eigen::MatrixXd::Ones(p,p);
    int n                   = z.rows() + 1;
    Eigen::VectorXd delmu;
    Eigen::MatrixXd gC_t, delC;
    
    for (int i = 0; i < T; i++) {
        Eigen::VectorXd th      = mu_t + C_t * xi.row(i).transpose();                    // II    - Set theta
        gmu = compute_nabmu_ELBO_Rcpp(delth_logpi, delthh, th, h, 
                                       lam0, z, n, a, p, T2);
        // Rcpp::List res = Rcpp::List::create(_["gmu"] = gmu);
        // return(res);
        Egmu = rho * Egmu + (1 - rho) * gmu.cwiseProduct(gmu);      // IV   - Accumulate gradients
        delmu = (((Edelmu + elip * Eigen::VectorXd::Ones(p)).cwiseSqrt().array() / (Egmu + elip * Eigen::VectorXd::Ones(p)).cwiseSqrt().array()) * gmu.array()).matrix(); // V     - Compute update
        mu_t = mu_t + delmu;                                        // VI    - Update mean
        Edelmu = rho * Edelmu + (1 - rho) * delmu.cwiseProduct(delmu);      // VII   - Accumulate updates
        gC_t = (compute_nabC_ELBO_Rcpp(gmu, xi.row(i), C_t)).triangularView<Eigen::Lower>();          // VIII  - Compute g_C^{t+1}
        EgC = rho * EgC + (1 - rho) * gC_t.cwiseProduct(gC_t);           // IX    - Accumulate gradients
        delC    = ((EdelC + elip * M).cwiseSqrt().array() / (EgC + elip * M).cwiseSqrt().array() * gC_t.array()).matrix();     // X     - Compute update
        C_t     = C_t + delC;                            // XI    - Update covariance Cholesky
        EgC     = rho * EgC + (1 - rho) * delC.cwiseProduct(delC);        // XII   - Accumulate updates
        
        // Store
        mu_arr.col(i+1) = mu_t;
        C_arr.push_back(C_t);
        
        if ((i+1)%500 == 0) { Rcpp::Rcout << "Iteration: " << i + 1 << std::endl;}
    }
    
    // return Rcpp::List::create(
    //         Rcpp::Named("mu_t") = mu_t,
    //         Rcpp::Named("C_t") = C_t);
    // cout << gC_t << endl;
    Rcpp::List res = Rcpp::List::create(
        _["mu_FC"]   = mu_t,
        _["C_FC"]    = C_t,
        _["mu_arr"] = mu_arr,
        _["C_arr"]  = C_arr,
        _["gmu"]    = gmu,
        _["Egmu"]   = Egmu,
        _["delmu"]  = delmu, 
        _["Edelmu"] = Edelmu, 
        _["gC_t"]   = gC_t, 
        _["EgC"]    = EgC, 
        _["delC"]   = delC
    );
    // std::vector<Eigen::MatrixXd> res = {mu_t, C_t, mu_arr, C_arr, Egmu, delmu, Edelmu, gC_t, EgC, delC};
    // std::vector<int> v = {(gmu * xi.row(i)).rows(), (gmu * xi.row(i)).cols(), gC_t.rows(), gC_t.cols(), gmu.rows(), gmu.cols(), xi.row(i).rows(), xi.row(i).cols()};
    return res;
}
    
// compute_GVA_Rcpp <- function(mu, C, h, delthh, delth_logpi, z, rho, elip, lam0, a, T, T2) {
// # Initialise values
// if (missing(T)) { T <- 10000 }
// if (missing(T2)) { T2 <- 500 }
// 
// p           <- nrow(C)
// Egmu        <- numeric(p)
// Edelmu      <- numeric(p)
// EgC         <- matrix(0, nrow = p, ncol = p)
// EdelC       <- matrix(0, nrow = p, ncol = p)
// mu_t        <- mu
// mu_arr      <- matrix(0,nrow = p, ncol = T+1)#array(dim = c(dim(mu_t), T+1))
// mu_arr[,1]  <- mu_t
// C_t         <- C        # Covariance Cholesky
// C_arr       <- array(dim = c(dim(C_t), T+1))
// C_arr[,,1]  <- C_t
// M           <- matrix(1,p,p)
// n           <- nrow(z) + 1

// xi          <- matrix(rnorm(T*p),T,p) 
// Egmu    <- rho * Egmu + (1 - rho) * gmu^2           # IV    - Accumulate gradients
// delmu   <- sqrt(Edelmu + elip * rep(1,p)) / 
//     sqrt(Egmu + elip * rep(1,p)) * gmu              # V     - Compute update
// mu_t    <- mu_t + delmu                             # VI    - Update mean
// Edelmu  <- rho * Edelmu + (1 - rho) * delmu^2       # VII   - Accumulate updates
// gC_t    <- compute_nabC_ELBO(gmu, xi[i,], C_t)      # VIII  - Compute g_C^{t+1}
// gC_t[upper.tri(gC_t)] <- 0                          #       - Set gC_t to lower triag matx
// EgC     <- rho * EgC + (1 - rho) * gC_t^2           # IX    - Accumulate gradients
// delC    <- sqrt(EdelC + elip * M) / 
//     sqrt(EgC + elip * M) * gC_t                     # X     - Compute update
// C_t     <- C_t + delC                               # XI    - Update covariance Cholesky
// EgC     <- rho * EgC + (1 - rho) * delC^2           # XII   - Accumulate updates



/*** R
# -----------------------------
# Main
# -----------------------------
# set.seed(1)
# x    <- runif(30, min = -5, max = 5)
# elip <- rnorm(30, mean = 0, sd = 1)
# y    <- 0.75 - x + elip
# lam0 <- matrix(c(0,0), nrow = 2)
# th   <- matrix(c(0.8277, -1.0050), nrow = 2)
# a <- 0.001
# z    <- cbind(x, y)
# h    <- function(z, th) {
#     xi <- z[1]
#     yi <- z[2]
#     h_zith <- c(yi - th[1] - th[2] * xi)
#     h_zith[2] <- xi * h_zith[1]
#     matrix(h_zith, nrow = 2)
# }
# 
# delthh    <- function(z, th) {
#     xi <- z[1]
#     matrix(c(-1, -xi, -xi, -xi^2), 2, 2)
# }
# 
# n <- 31
# reslm <- lm(y ~ x)
# mu <- matrix(unname(reslm$coefficients),2,1)
# C_0 <- unname(t(chol(vcov(reslm))))
# 
# # warning("Treating delth_logpi as column vector")
# delth_logpi <- function(theta) {-0.0001 * mu}
# elip <- 10^-5
# T <- 10
# T2 <- 500
# rho <- 0.9
# 
# p           <- nrow(C)
# Egmu        <- numeric(p)
# Edelmu      <- numeric(p)
# EgC         <- matrix(0, nrow = p, ncol = p)
# EdelC       <- matrix(0, nrow = p, ncol = p)
# mu_t        <- mu
# mu_arr      <- matrix(0,nrow = p, ncol = T+1)#array(dim = c(dim(mu_t), T+1))
# mu_arr[,1]  <- mu_t
# C_t         <- C        # Covariance Cholesky
# C_arr       <- array(dim = c(dim(C_t), T+1))
# C_arr[,,1]  <- C_t
# M           <- matrix(1,p,p)
# n           <- nrow(z) + 1
# 
# xi          <- matrix(rnorm(T*p),T,p)                   # I     - Draw xi
# 
# i <- 1
# compute_nabmu_ELBO <- function(delth_logpi, delthh, theta, n, h, lam0, a, z, T2) { 
#     res <- compute_AEL_R(theta, h, lam0, a, z, T2, 1) # list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth)
#     lambda <- res$"lambda"
#     h_arr <- res$"h_arr"
#     hznth <- h_arr[,,n]
#     
#     # Calculate gradient LogAEL with respect to theta
#     nabth_logAEL <- 0 # Matrix
#     for (i in 1:(n-1)) {
#         nabth_logAEL <- nabth_logAEL - (1/(1 + t(lambda) %*% h_arr[,,i]) - (a/(n-1)) / (1 + t(lambda) %*% hznth))[1] * (t(delthh(z[i], theta)) %*% lambda)
#     }
#     nabmu_ELBO <- nabth_logAEL + delth_logpi(theta)
# }
# 
# th      <- mu_t + C_t %*% xi[i,]                    # II    - Set theta
# gmu     <- compute_nabmu_ELBO(delth_logpi, delthh, 
#                               th, n, h, lam0, 
#                               a, z, T2)             # III   - Compute g_{mu}^{t+1}
# output <- compute_GVA_Rcpp_inner(rho, elip, Egmu, Edelmu, EgC, EdelC, gmu, mu_t, C_t, xi, p, i-1)
# Egmu    <- rho * Egmu + (1 - rho) * gmu^2           # IV    - Accumulate gradients
# delmu   <- sqrt(Edelmu + elip * rep(1,p)) /
#     sqrt(Egmu + elip * rep(1,p)) * gmu              # V     - Compute update
# mu_t    <- mu_t + delmu                             # VI    - Update mean
# Edelmu  <- rho * Edelmu + (1 - rho) * delmu^2       # VII   - Accumulate updates
# browser()
# gC_t    <- compute_nabC_ELBO(gmu, xi[i,], C_t)      # VIII  - Compute g_C^{t+1}
# browser()
# gC_t[upper.tri(gC_t)] <- 0   
# 
# # ans <-compute_GVA_Rcpp_inner(th, h, lam0, a, z)
# 

ans
*/
