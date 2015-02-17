// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// count_cpp
NumericMatrix count_cpp(NumericMatrix data, bool unique = true, bool freq = true);
RcppExport SEXP equi_count_cpp(SEXP dataSEXP, SEXP uniqueSEXP, SEXP freqSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP );
        Rcpp::traits::input_parameter< bool >::type unique(uniqueSEXP );
        Rcpp::traits::input_parameter< bool >::type freq(freqSEXP );
        NumericMatrix __result = count_cpp(data, unique, freq);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
