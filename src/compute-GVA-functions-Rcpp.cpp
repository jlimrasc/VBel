// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <random>
#include "compute-AEL-functions-Rcpp.h"

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
Eigen::MatrixXd compute_nabC_ELBO_Rcpp(const Eigen::VectorXd &gmu, const Eigen::VectorXd &xi, const Eigen::MatrixXd &C_t) {
    // Store wi for use in D in P
    Eigen::MatrixXd nabC_ELBO;
    nabC_ELBO = gmu * xi.transpose() + C_t.diagonal().cwiseInverse().asDiagonal().toDenseMatrix(); // .row changes it to a row vector rather than column but VectorXd makes it a column again so transpose
    return nabC_ELBO;
}

Eigen::VectorXd compute_nabmu_ELBO_Rcpp(Rcpp::Function delth_logpi, Rcpp::Function delthh, 
                                        const Eigen::VectorXd &theta, Rcpp::Function h, 
                                        const Eigen::VectorXd &lam0, const Eigen::MatrixXd &z, 
                                        int n, double a, int p, int T2, int i_out) {
    Rcpp::List res = compute_AEL_Rcpp_inner(theta, h, lam0, a, z, T2); // list("log_AEL" = log_AEL[1, 1], "lambda" = lambda, "h_arr" = h_arr, "H" = H_Zth)
    Eigen::MatrixXd lambda              = res["lambda"];
    std::vector<Eigen::VectorXd> h_arr  = res["h_arr"];
    Eigen::VectorXd hznth               = h_arr[n-1];
    if (lambda.hasNaN()) {
        Rcpp::Rcout << "\nLambda Iteration: " << i_out << " has NaN:\n" << lambda << std::endl;
    }
    

    // Calculate gradient LogAEL with respect to theta
    Eigen::VectorXd nabth_logAEL = Eigen::VectorXd::Zero(p); // Vector
    Eigen::MatrixXd delthh_zith;
    for (int i = 0; i < n-1; i++) {
        nabth_logAEL = nabth_logAEL - ((1/(1 + (lambda.transpose() * h_arr[i]).array()) - (a/(n-1)) / (1 + (lambda.transpose() * hznth).array()))(0,0) * (Rcpp::as<Eigen::MatrixXd>(delthh(z.row(i).transpose(),theta)).transpose() * lambda).array()).matrix();
    }
    
    Eigen::VectorXd nabmu_ELBO = (nabth_logAEL.array() + Rcpp::as<Eigen::VectorXd>(delth_logpi(theta)).array()).matrix();
    return nabmu_ELBO;
}

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
        
    std::vector<Eigen::MatrixXd> res = {mu_t, C_t, Egmu, delmu, Edelmu, gC_t, EgC, delC, gmu};
    return res;
}
    
// [[Rcpp::export]]
Rcpp::List compute_GVA_Rcpp_inner_full(
        Eigen::VectorXd mu, Eigen::MatrixXd C, Rcpp::Function h, Rcpp::Function delthh,
        Rcpp::Function delth_logpi, Eigen::MatrixXd z, Eigen::VectorXd lam0, 
        double rho, double elip, double a, int T, int T2, int p, int verbosity){

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
    
    // Set up for normal distribution
    // From: https://eigen.tuxfamily.org/dox/classEigen_1_1DenseBase.html#a3340c9b997f5b53a0131cf927f93b54c
    std::default_random_engine generator{static_cast<long unsigned int>(time(0))};
    std::normal_distribution<double> distribution(0,1);
    auto normal_dist = [&] (double) {return distribution(generator);};

    Eigen::MatrixXd xi = Eigen::MatrixXd::NullaryExpr(T, p, normal_dist );                  // I    - Draw xi
    
    for (int i = 0; i < T; i++) {
        Eigen::VectorXd th = mu_t + C_t * xi.row(i).transpose();                           // II   - Set theta
        gmu = compute_nabmu_ELBO_Rcpp(delth_logpi, delthh, th, h, 
                                       lam0, z, n, a, p, T2, i);
        // Rcpp::List res = Rcpp::List::create(_["gmu"] = gmu);
        // return(res);
        Egmu = rho * Egmu + (1 - rho) * gmu.cwiseProduct(gmu);                              // IV   - Accumulate gradients
        delmu = (((Edelmu + elip * Eigen::VectorXd::Ones(p)).cwiseSqrt().array() / 
            (Egmu + elip * Eigen::VectorXd::Ones(p)).cwiseSqrt().array()) * 
            gmu.array()).matrix();                                                          // V     - Compute update
        mu_t = mu_t + delmu;                                                                // VI   - Update mean
        Edelmu = rho * Edelmu + (1 - rho) * delmu.cwiseProduct(delmu);                      // VII  - Accumulate updates
        gC_t = (compute_nabC_ELBO_Rcpp(gmu, xi.row(i), C_t)).triangularView<Eigen::Lower>();// VIII - Compute g_C^{t+1}
        EgC = rho * EgC + (1 - rho) * gC_t.cwiseProduct(gC_t);                              // IX   - Accumulate gradients
        delC    = ((EdelC + elip * M).cwiseSqrt().array() / 
            (EgC + elip * M).cwiseSqrt().array() * gC_t.array()).matrix();                  // X     - Compute update
        C_t     = C_t + delC;                                                               // XI   - Update covariance Cholesky
        EgC     = rho * EgC + (1 - rho) * delC.cwiseProduct(delC);                          // XII  - Accumulate updates
        
        // Store
        mu_arr.col(i+1) = mu_t;
        C_arr.push_back(C_t);
        
        if (verbosity && (i+1) % verbosity == 0) { Rcpp::Rcout << "Iteration: " << i + 1 << std::endl;}
        if (mu_t.hasNaN() || C_t.hasNaN()) {
            Rcpp::Rcout << "NaN found in calculation. Terminating early. Please run function again.\n";
            break;
        }
    }
    
    Rcpp::List res = Rcpp::List::create(
        Rcpp::_["mu_FC"]   = mu_t,
        Rcpp::_["C_FC"]    = C_t,
        Rcpp::_["mu_arr"] = mu_arr,
        Rcpp::_["C_arr"]  = C_arr,
        Rcpp::_["gmu"]    = gmu,
        Rcpp::_["Egmu"]   = Egmu,
        Rcpp::_["delmu"]  = delmu, 
        Rcpp::_["Edelmu"] = Edelmu, 
        Rcpp::_["gC_t"]   = gC_t, 
        Rcpp::_["EgC"]    = EgC, 
        Rcpp::_["delC"]   = delC
    );
    return res;
}
    



/*** R
# R code not needed, too much dependencies on other functions
*/
