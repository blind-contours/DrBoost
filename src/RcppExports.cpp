// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// find_best_logical_region_cpp
List find_best_logical_region_cpp(const NumericMatrix& X, const NumericVector& residuals, const IntegerVector& used_rows, double min_obs_pct);
RcppExport SEXP _DrBoost_find_best_logical_region_cpp(SEXP XSEXP, SEXP residualsSEXP, SEXP used_rowsSEXP, SEXP min_obs_pctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type residuals(residualsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type used_rows(used_rowsSEXP);
    Rcpp::traits::input_parameter< double >::type min_obs_pct(min_obs_pctSEXP);
    rcpp_result_gen = Rcpp::wrap(find_best_logical_region_cpp(X, residuals, used_rows, min_obs_pct));
    return rcpp_result_gen;
END_RCPP
}
// find_topk_candidates_beam_cpp
List find_topk_candidates_beam_cpp(const NumericMatrix& X, const NumericVector& residuals, const IntegerVector& coverageCount, const int max_overlap, const NumericMatrix& thresholds_list, double min_obs_pct, double max_obs_frac, int K1, int K2, int K3, int K4, int topK);
RcppExport SEXP _DrBoost_find_topk_candidates_beam_cpp(SEXP XSEXP, SEXP residualsSEXP, SEXP coverageCountSEXP, SEXP max_overlapSEXP, SEXP thresholds_listSEXP, SEXP min_obs_pctSEXP, SEXP max_obs_fracSEXP, SEXP K1SEXP, SEXP K2SEXP, SEXP K3SEXP, SEXP K4SEXP, SEXP topKSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type residuals(residualsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type coverageCount(coverageCountSEXP);
    Rcpp::traits::input_parameter< const int >::type max_overlap(max_overlapSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type thresholds_list(thresholds_listSEXP);
    Rcpp::traits::input_parameter< double >::type min_obs_pct(min_obs_pctSEXP);
    Rcpp::traits::input_parameter< double >::type max_obs_frac(max_obs_fracSEXP);
    Rcpp::traits::input_parameter< int >::type K1(K1SEXP);
    Rcpp::traits::input_parameter< int >::type K2(K2SEXP);
    Rcpp::traits::input_parameter< int >::type K3(K3SEXP);
    Rcpp::traits::input_parameter< int >::type K4(K4SEXP);
    Rcpp::traits::input_parameter< int >::type topK(topKSEXP);
    rcpp_result_gen = Rcpp::wrap(find_topk_candidates_beam_cpp(X, residuals, coverageCount, max_overlap, thresholds_list, min_obs_pct, max_obs_frac, K1, K2, K3, K4, topK));
    return rcpp_result_gen;
END_RCPP
}
// find_best_binary_rule_3way_topK
List find_best_binary_rule_3way_topK(const NumericMatrix& X, const NumericVector& residuals, const IntegerVector& unassigned, int K_twoWay, int K_threeWay);
RcppExport SEXP _DrBoost_find_best_binary_rule_3way_topK(SEXP XSEXP, SEXP residualsSEXP, SEXP unassignedSEXP, SEXP K_twoWaySEXP, SEXP K_threeWaySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type residuals(residualsSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type unassigned(unassignedSEXP);
    Rcpp::traits::input_parameter< int >::type K_twoWay(K_twoWaySEXP);
    Rcpp::traits::input_parameter< int >::type K_threeWay(K_threeWaySEXP);
    rcpp_result_gen = Rcpp::wrap(find_best_binary_rule_3way_topK(X, residuals, unassigned, K_twoWay, K_threeWay));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DrBoost_find_best_logical_region_cpp", (DL_FUNC) &_DrBoost_find_best_logical_region_cpp, 4},
    {"_DrBoost_find_topk_candidates_beam_cpp", (DL_FUNC) &_DrBoost_find_topk_candidates_beam_cpp, 12},
    {"_DrBoost_find_best_binary_rule_3way_topK", (DL_FUNC) &_DrBoost_find_best_binary_rule_3way_topK, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_DrBoost(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
