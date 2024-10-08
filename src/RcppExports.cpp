// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// colSortC
arma::mat colSortC(arma::mat X);
RcppExport SEXP _pARI_colSortC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(colSortC(X));
    return rcpp_result_gen;
END_RCPP
}
// lambdaCalibrate
NumericVector lambdaCalibrate(arma::mat X, arma::vec alpha, double delta, std::string family, double m);
RcppExport SEXP _pARI_lambdaCalibrate(SEXP XSEXP, SEXP alphaSEXP, SEXP deltaSEXP, SEXP familySEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(lambdaCalibrate(X, alpha, delta, family, m));
    return rcpp_result_gen;
END_RCPP
}
// permDiscoveries
arma::vec permDiscoveries(NumericVector ix, NumericVector cv, NumericVector praw);
RcppExport SEXP _pARI_permDiscoveries(SEXP ixSEXP, SEXP cvSEXP, SEXP prawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ix(ixSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type praw(prawSEXP);
    rcpp_result_gen = Rcpp::wrap(permDiscoveries(ix, cv, praw));
    return rcpp_result_gen;
END_RCPP
}
// permT
arma::mat permT(arma::mat X, double B, arma::vec label);
RcppExport SEXP _pARI_permT(SEXP XSEXP, SEXP BSEXP, SEXP labelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type label(labelSEXP);
    rcpp_result_gen = Rcpp::wrap(permT(X, B, label));
    return rcpp_result_gen;
END_RCPP
}
// rowSortC
arma::mat rowSortC(arma::mat X);
RcppExport SEXP _pARI_rowSortC(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSortC(X));
    return rcpp_result_gen;
END_RCPP
}
// signFlip
arma::mat signFlip(arma::mat X, double B);
RcppExport SEXP _pARI_signFlip(SEXP XSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(signFlip(X, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pARI_colSortC", (DL_FUNC) &_pARI_colSortC, 1},
    {"_pARI_lambdaCalibrate", (DL_FUNC) &_pARI_lambdaCalibrate, 5},
    {"_pARI_permDiscoveries", (DL_FUNC) &_pARI_permDiscoveries, 3},
    {"_pARI_permT", (DL_FUNC) &_pARI_permT, 3},
    {"_pARI_rowSortC", (DL_FUNC) &_pARI_rowSortC, 1},
    {"_pARI_signFlip", (DL_FUNC) &_pARI_signFlip, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_pARI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
