// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_dotL2
double cpp_dotL2(NumericVector f1, NumericVector f2, NumericVector grid);
RcppExport SEXP _densityFPCA_cpp_dotL2(SEXP f1SEXP, SEXP f2SEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dotL2(f1, f2, grid));
    return rcpp_result_gen;
END_RCPP
}
// cpp_calCdf
NumericVector cpp_calCdf(NumericVector pdf, NumericVector grid);
RcppExport SEXP _densityFPCA_cpp_calCdf(SEXP pdfSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pdf(pdfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_calCdf(pdf, grid));
    return rcpp_result_gen;
END_RCPP
}
// cpp_distWass
double cpp_distWass(NumericVector cdf1, NumericVector cdf2, NumericVector grid, double p);
RcppExport SEXP _densityFPCA_cpp_distWass(SEXP cdf1SEXP, SEXP cdf2SEXP, SEXP gridSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type cdf1(cdf1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdf2(cdf2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_distWass(cdf1, cdf2, grid, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_densityFPCA_cpp_dotL2", (DL_FUNC) &_densityFPCA_cpp_dotL2, 3},
    {"_densityFPCA_cpp_calCdf", (DL_FUNC) &_densityFPCA_cpp_calCdf, 2},
    {"_densityFPCA_cpp_distWass", (DL_FUNC) &_densityFPCA_cpp_distWass, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_densityFPCA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
