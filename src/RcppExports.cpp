// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pinp_num
LogicalVector pinp_num(NumericVector x, NumericVector y);
RcppExport SEXP _buffr_pinp_num(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(pinp_num(x, y));
    return rcpp_result_gen;
END_RCPP
}
// paste_int
NumericVector paste_int(IntegerVector row, IntegerVector col);
RcppExport SEXP _buffr_paste_int(SEXP rowSEXP, SEXP colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type row(rowSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type col(colSEXP);
    rcpp_result_gen = Rcpp::wrap(paste_int(row, col));
    return rcpp_result_gen;
END_RCPP
}
// dupe
LogicalVector dupe(NumericVector vec);
RcppExport SEXP _buffr_dupe(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(dupe(vec));
    return rcpp_result_gen;
END_RCPP
}
// rowSumsC
IntegerMatrix rowSumsC(IntegerMatrix mat);
RcppExport SEXP _buffr_rowSumsC(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSumsC(mat));
    return rcpp_result_gen;
END_RCPP
}
// f1
LogicalVector f1(IntegerMatrix mat, int max_row);
RcppExport SEXP _buffr_f1(SEXP matSEXP, SEXP max_rowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type max_row(max_rowSEXP);
    rcpp_result_gen = Rcpp::wrap(f1(mat, max_row));
    return rcpp_result_gen;
END_RCPP
}
// f2
LogicalVector f2(IntegerMatrix mat, int max_col);
RcppExport SEXP _buffr_f2(SEXP matSEXP, SEXP max_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type max_col(max_colSEXP);
    rcpp_result_gen = Rcpp::wrap(f2(mat, max_col));
    return rcpp_result_gen;
END_RCPP
}
// f3
IntegerVector f3(IntegerMatrix mat);
RcppExport SEXP _buffr_f3(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(f3(mat));
    return rcpp_result_gen;
END_RCPP
}
// f4
IntegerVector f4(IntegerVector vec);
RcppExport SEXP _buffr_f4(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(f4(vec));
    return rcpp_result_gen;
END_RCPP
}
// expand_gridC_2
IntegerMatrix expand_gridC_2(IntegerVector vec1_1, IntegerVector vec1_2, IntegerVector vec2_1, IntegerVector vec2_2);
RcppExport SEXP _buffr_expand_gridC_2(SEXP vec1_1SEXP, SEXP vec1_2SEXP, SEXP vec2_1SEXP, SEXP vec2_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type vec1_1(vec1_1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vec1_2(vec1_2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vec2_1(vec2_1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vec2_2(vec2_2SEXP);
    rcpp_result_gen = Rcpp::wrap(expand_gridC_2(vec1_1, vec1_2, vec2_1, vec2_2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_buffr_pinp_num", (DL_FUNC) &_buffr_pinp_num, 2},
    {"_buffr_paste_int", (DL_FUNC) &_buffr_paste_int, 2},
    {"_buffr_dupe", (DL_FUNC) &_buffr_dupe, 1},
    {"_buffr_rowSumsC", (DL_FUNC) &_buffr_rowSumsC, 1},
    {"_buffr_f1", (DL_FUNC) &_buffr_f1, 2},
    {"_buffr_f2", (DL_FUNC) &_buffr_f2, 2},
    {"_buffr_f3", (DL_FUNC) &_buffr_f3, 1},
    {"_buffr_f4", (DL_FUNC) &_buffr_f4, 1},
    {"_buffr_expand_gridC_2", (DL_FUNC) &_buffr_expand_gridC_2, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_buffr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
