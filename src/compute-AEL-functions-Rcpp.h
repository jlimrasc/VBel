#ifndef COMPUTE_AEL_FUNCTIONS_RCPP
#define COMPUTE_AEL_FUNCTIONS_RCPP

#include <stdio.h>
#include <cmath>
#include <RcppEigen.h>

Eigen::VectorXd get_wi_arr_Rcpp(std::vector<Eigen::VectorXd> h_list, Eigen::VectorXd tild_lam, int n);
Eigen::VectorXd get_dF_Rcpp(std::vector<Eigen::VectorXd> h_list, Eigen::VectorXd wi_arr, int n, int d);
Eigen::MatrixXd get_d2F_Rcpp(std::vector<Eigen::VectorXd> h_list, Eigen::MatrixXd H_Zth, Eigen::VectorXd wi_arr, int n, int d);
Eigen::MatrixXd compute_lambda_Rcpp(std::vector<Eigen::VectorXd> h_list, Eigen::MatrixXd H_Zth, Eigen::VectorXd lam0, double a, int T, int n, int d);
Rcpp::List compute_AEL_Rcpp_inner_main(std::vector<Eigen::VectorXd> h_list, Eigen::MatrixXd H_Zth, Eigen::VectorXd lam0, double a, int T, int n, int d);
Rcpp::List compute_AEL_Rcpp_inner(Eigen::VectorXd th, Rcpp::Function h, Eigen::VectorXd lam0, double a, Eigen::MatrixXd z, int T = 100);

#endif