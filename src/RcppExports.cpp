// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Cpp_fwpopsim
List Cpp_fwpopsim(int G, IntegerMatrix H0, IntegerVector N0, NumericVector alpha, List mutmodel, bool SNP, IntegerVector save_gs, bool progress, bool trace);
RcppExport SEXP fwsim_Cpp_fwpopsim(SEXP GSEXP, SEXP H0SEXP, SEXP N0SEXP, SEXP alphaSEXP, SEXP mutmodelSEXP, SEXP SNPSEXP, SEXP save_gsSEXP, SEXP progressSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type G(GSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type H0(H0SEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type N0(N0SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< List >::type mutmodel(mutmodelSEXP );
        Rcpp::traits::input_parameter< bool >::type SNP(SNPSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type save_gs(save_gsSEXP );
        Rcpp::traits::input_parameter< bool >::type progress(progressSEXP );
        Rcpp::traits::input_parameter< bool >::type trace(traceSEXP );
        List __result = Cpp_fwpopsim(G, H0, N0, alpha, mutmodel, SNP, save_gs, progress, trace);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}