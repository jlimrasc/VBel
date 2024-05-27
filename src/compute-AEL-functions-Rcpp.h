#ifndef COMPUTE_AEL_FUNCTIONS_RCPP
#define COMPUTE_AEL_FUNCTIONS_RCPP

#include <stdio.h>
#include <cmath>
#include <RcppEigen.h>

Eigen::VectorXd get_wi_arr_Rcpp(const std::vector<Eigen::VectorXd> &h_list, const Eigen::VectorXd &tild_lam, int n);
Eigen::VectorXd get_dF_Rcpp(const std::vector<Eigen::VectorXd> &h_list, const Eigen::VectorXd &wi_arr, int n, int d);
Eigen::MatrixXd get_d2F_Rcpp(const std::vector<Eigen::VectorXd> &h_list, const Eigen::MatrixXd &H_Zth, const Eigen::VectorXd &wi_arr, int n, int d);
Eigen::MatrixXd compute_lambda_Rcpp(const std::vector<Eigen::VectorXd> &h_list, const Eigen::MatrixXd &H_Zth, const Eigen::VectorXd &lam0, double a, int T, int n, int d);
Rcpp::List compute_AEL_Rcpp_inner_main(const std::vector<Eigen::VectorXd> &h_list, const Eigen::MatrixXd &H_Zth, const Eigen::VectorXd &lam0, double a, int T, int n, int d);
Rcpp::List compute_AEL_Rcpp_inner(const Eigen::VectorXd &th, Rcpp::Function h, const Eigen::VectorXd &lam0, double a, const Eigen::MatrixXd &z, int T = 100);

#endif