// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_lambda_Rcpp
Eigen::MatrixXd compute_lambda_Rcpp(std::vector<Eigen::VectorXd> h_list, Eigen::MatrixXd H_Zth, Eigen::VectorXd lam0, double a, int T, int n, int d);
RcppExport SEXP _VBel_compute_lambda_Rcpp(SEXP h_listSEXP, SEXP H_ZthSEXP, SEXP lam0SEXP, SEXP aSEXP, SEXP TSEXP, SEXP nSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<Eigen::VectorXd> >::type h_list(h_listSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type H_Zth(H_ZthSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_lambda_Rcpp(h_list, H_Zth, lam0, a, T, n, d));
    return rcpp_result_gen;
END_RCPP
}
// compute_AEL_Rcpp_inner
Rcpp::List compute_AEL_Rcpp_inner(Eigen::VectorXd th, Rcpp::Function h, Eigen::VectorXd lam0, double a, Eigen::MatrixXd z, int T);
RcppExport SEXP _VBel_compute_AEL_Rcpp_inner(SEXP thSEXP, SEXP hSEXP, SEXP lam0SEXP, SEXP aSEXP, SEXP zSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type th(thSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type h(hSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_AEL_Rcpp_inner(th, h, lam0, a, z, T));
    return rcpp_result_gen;
END_RCPP
}
// compute_AEL_Rcpp_inner_prez
Rcpp::List compute_AEL_Rcpp_inner_prez(Eigen::VectorXd th, Eigen::MatrixXd H_Zth, Eigen::VectorXd lam0, double a, Eigen::MatrixXd z, int T);
RcppExport SEXP _VBel_compute_AEL_Rcpp_inner_prez(SEXP thSEXP, SEXP H_ZthSEXP, SEXP lam0SEXP, SEXP aSEXP, SEXP zSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type th(thSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type H_Zth(H_ZthSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_AEL_Rcpp_inner_prez(th, H_Zth, lam0, a, z, T));
    return rcpp_result_gen;
END_RCPP
}
// compute_GVA_Rcpp_inner_IVtoXII
std::vector<Eigen::MatrixXd> compute_GVA_Rcpp_inner_IVtoXII(const double rho, const double elip, Eigen::VectorXd Egmu, Eigen::VectorXd Edelmu, Eigen::MatrixXd EgC, Eigen::MatrixXd EdelC, Eigen::VectorXd gmu, Eigen::VectorXd mu_t, Eigen::MatrixXd C_t, const Eigen::MatrixXd xi, Eigen::MatrixXd M, const int p, const int i);
RcppExport SEXP _VBel_compute_GVA_Rcpp_inner_IVtoXII(SEXP rhoSEXP, SEXP elipSEXP, SEXP EgmuSEXP, SEXP EdelmuSEXP, SEXP EgCSEXP, SEXP EdelCSEXP, SEXP gmuSEXP, SEXP mu_tSEXP, SEXP C_tSEXP, SEXP xiSEXP, SEXP MSEXP, SEXP pSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const double >::type elip(elipSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Egmu(EgmuSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Edelmu(EdelmuSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type EgC(EgCSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type EdelC(EdelCSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gmu(gmuSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu_t(mu_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type C_t(C_tSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_GVA_Rcpp_inner_IVtoXII(rho, elip, Egmu, Edelmu, EgC, EdelC, gmu, mu_t, C_t, xi, M, p, i));
    return rcpp_result_gen;
END_RCPP
}
// compute_GVA_Rcpp_inner_full
Rcpp::List compute_GVA_Rcpp_inner_full(Eigen::VectorXd mu, Eigen::MatrixXd C, Rcpp::Function h, Rcpp::Function delthh, Rcpp::Function delth_logpi, Eigen::MatrixXd z, Eigen::VectorXd lam0, Eigen::MatrixXd xi, double rho, double elip, double a, int T, int T2, int p, int verbosity);
RcppExport SEXP _VBel_compute_GVA_Rcpp_inner_full(SEXP muSEXP, SEXP CSEXP, SEXP hSEXP, SEXP delthhSEXP, SEXP delth_logpiSEXP, SEXP zSEXP, SEXP lam0SEXP, SEXP xiSEXP, SEXP rhoSEXP, SEXP elipSEXP, SEXP aSEXP, SEXP TSEXP, SEXP T2SEXP, SEXP pSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type delthh(delthhSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type delth_logpi(delth_logpiSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type z(zSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lam0(lam0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type elip(elipSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type T2(T2SEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_GVA_Rcpp_inner_full(mu, C, h, delthh, delth_logpi, z, lam0, xi, rho, elip, a, T, T2, p, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_hello_world
Eigen::MatrixXd rcppeigen_hello_world();
RcppExport SEXP _VBel_rcppeigen_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcppeigen_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_outerproduct
Eigen::MatrixXd rcppeigen_outerproduct(const Eigen::VectorXd& x);
RcppExport SEXP _VBel_rcppeigen_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_innerproduct
double rcppeigen_innerproduct(const Eigen::VectorXd& x);
RcppExport SEXP _VBel_rcppeigen_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_bothproducts
Rcpp::List rcppeigen_bothproducts(const Eigen::VectorXd& x);
RcppExport SEXP _VBel_rcppeigen_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VBel_compute_lambda_Rcpp", (DL_FUNC) &_VBel_compute_lambda_Rcpp, 7},
    {"_VBel_compute_AEL_Rcpp_inner", (DL_FUNC) &_VBel_compute_AEL_Rcpp_inner, 6},
    {"_VBel_compute_AEL_Rcpp_inner_prez", (DL_FUNC) &_VBel_compute_AEL_Rcpp_inner_prez, 6},
    {"_VBel_compute_GVA_Rcpp_inner_IVtoXII", (DL_FUNC) &_VBel_compute_GVA_Rcpp_inner_IVtoXII, 13},
    {"_VBel_compute_GVA_Rcpp_inner_full", (DL_FUNC) &_VBel_compute_GVA_Rcpp_inner_full, 15},
    {"_VBel_rcppeigen_hello_world", (DL_FUNC) &_VBel_rcppeigen_hello_world, 0},
    {"_VBel_rcppeigen_outerproduct", (DL_FUNC) &_VBel_rcppeigen_outerproduct, 1},
    {"_VBel_rcppeigen_innerproduct", (DL_FUNC) &_VBel_rcppeigen_innerproduct, 1},
    {"_VBel_rcppeigen_bothproducts", (DL_FUNC) &_VBel_rcppeigen_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_VBel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
