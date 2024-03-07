// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <stdio.h>
#include <cmath>
using Eigen::MatrixXd;      // variable size martix, double precision
using Eigen::VectorXd;      // variable size martix, double precision
using Eigen::DiagonalMatrix;// diagonal matrix
using namespace std;
// using namespace Rcpp;

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
VectorXd get_wi_arr_Rcpp(vector<VectorXd> h_list, VectorXd tild_lam, int n) {
    // Store wi for use in D in P
    VectorXd wi_arr(n);
    for (int i = 0; i < n; i++) {
        wi_arr(i) = pow((1 + tild_lam.transpose() * h_list[i]),-1);
    }
    return wi_arr;
}

// -----------------------------
// Compute dF
// -----------------------------
VectorXd get_dF_Rcpp(vector<VectorXd> h_list, VectorXd wi_arr, int n, int d) {
    // #' @param wi_arr Array of vectors of wi(th, lam~)
    double wi, vi;
    VectorXd dF = VectorXd::Zero(d);
    
    for (int i = 0; i < n; i++) {
        
        // Evaluate vi
        wi = wi_arr[i];
        if (pow(wi,-1) >= 1.0/(double)n) {
            vi = wi;
        } else {
            vi = 2 * n - pow(n,2) / wi;
        }
        dF += vi * h_list[i]; // Calculate sum for this iteration of dF
    }
    return dF;
}

// -----------------------------
// Compute d2F
// -----------------------------
MatrixXd get_d2F_Rcpp(vector<VectorXd> h_list, MatrixXd H_Zth, VectorXd wi_arr, int n, int d) {
    // Build diagonal matrix
    // Diagonal matrix changes with respect to lambda's current guess so need 
    // to recompute each iteration
    VectorXd D_arr(n);
    MatrixXd P(d,d), D(n,n);
    double wi, vi2;
    
    // Evaluate vi2
    for (int i = 0; i < n; i++) {
        wi = wi_arr[i];
        if (pow(wi,-1) >= 1.0/(double)n) {
            vi2 = wi;
        } else {
            vi2 = n;
        }
        D_arr[i] = pow(vi2,2);
    }
    D = D_arr.asDiagonal();
    
    // Calculate P (i.e. d2F)
    return -H_Zth.transpose() * D * H_Zth;
}


// [[Rcpp::export]]
MatrixXd compute_lambda_Rcpp(vector<VectorXd> h_list, MatrixXd H_Zth, VectorXd lam0, double a, int T, int n, int d) {
    // # -----------------------------
    // # Compute lambda using modified Newton-Raphson
    // # -----------------------------
    VectorXd lam_prev = lam0;
    VectorXd wi_arr(n), dF(d);
    MatrixXd P(d,d);
    
    for (int i = 0; i < T; i++) {
        wi_arr = get_wi_arr_Rcpp(h_list, lam_prev, n); // Wi
        
        dF = get_dF_Rcpp(h_list, wi_arr, n, d); // dF
        
        P = get_d2F_Rcpp(h_list, H_Zth, wi_arr, n, d); // d2F
        
        lam_prev = lam_prev - P.inverse() * dF;
        
        // cout << i << " " << dF << " " << P << " " << lam_prev << endl;
    }
    
    return lam_prev;
}

// [[Rcpp::export]]
double compute_AEL_Rcpp(VectorXd th, Rcpp::Function h, 
                        VectorXd lam0, double a, MatrixXd z, int T = 100) {
    // # -----------------------------
    // # Default values
    // # -----------------------------
    //       if (missing(T)) {
    //         T <- 100
    //       }
    //       
    // # -----------------------------
    // # Starting variables (h(zi,th), h(zn,th), H)
    // # -----------------------------
    int n = z.rows() + 1;
    int d = z.cols();
    vector<VectorXd> h_list(n);
    MatrixXd H_Zth(n,d), h_znth(1, 2);
    VectorXd h_zith(d);
    
    for (int i = 0; i < n - 1; i++) {
        h_zith = Rcpp::as<VectorXd>(h(z.row(i).transpose(), th)); // Might be slower,
        // use R to build instead?
        h_list[i] = h_zith;
        
        H_Zth.row(i) = h_zith.transpose();
    }
    
    // Need brackets around  fraction for some reason? Maybe matrix mult takes priority?
    h_znth = -(a / (n - 1)) * H_Zth.colwise().sum();
    H_Zth.row(n - 1) = h_znth;
    
    h_list[n - 1] = h_znth.row(0);
    
    // -----------------------------
    // Lambda Calculation
    // -----------------------------
    MatrixXd lambda = compute_lambda_Rcpp(h_list, H_Zth, lam0, a, T, n, d);
    
    // -----------------------------
    // AEL Calculation
    // -----------------------------
    double accum = 0, log_AEL;
    
    for (int i = 0; i < n; i++) {
        accum += log(1.0 + (lambda.transpose() * h_list[i])[0]); // Get double vs mat hence [0]
    }
    log_AEL = -(accum + n * log(n));
    
    return log_AEL;
}

/*** R
# -----------------------------
# Main
# -----------------------------
set.seed(1)
x    <- runif(30, min = -5, max = 5)
elip <- rnorm(30, mean = 0, sd = 1)
y    <- 0.75 - x + elip
lam0 <- matrix(c(0,0), nrow = 2)
th   <- matrix(c(0.8277, -1.0050), nrow = 2)
a <- 0.001
z    <- cbind(x, y)
h    <- function(z, th) {
    xi <- z[1]
    yi <- z[2]
    h_zith <- c(yi - th[1] - th[2] * xi)
    h_zith[2] <- xi * h_zith[1]
    matrix(h_zith, nrow = 2)
}

# compute_AEL_Rcpp <- function(th, h, lam0, a, z){
#   
# }
# ans <-tester(th, h, lam0, a, z)
ans <-compute_AEL_Rcpp(th, h, lam0, a, z)


ans
*/