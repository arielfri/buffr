#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector pinp_num(NumericVector x, NumericVector y) {
  LogicalVector out = in(x,y);
  return out;
}

// [[Rcpp::export]]
NumericVector paste_int(IntegerVector row, IntegerVector col) {
  int n = row.size();
  NumericVector out(n);
  int max_num = max(col);
  int digits = 0;
  while (max_num) {
    max_num /= 10;
    digits++;
  }
  for (int i = 0; i < n; ++i) {
    out[i] = row[i]+col[i]/pow(10, digits);
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector dupe(NumericVector vec) {
  LogicalVector out = duplicated(vec);
  return out;
}

// [[Rcpp::export]]
IntegerMatrix rowSumsC(IntegerMatrix mat) {
  int n = mat.nrow();
  IntegerMatrix out(n,2);
  for (int i = 0; i < n; ++i) {
    out(i, 0) = mat(i, 0)+mat(i, 2);
    out(i, 1) = mat(i, 1)+mat(i, 3);
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector f1(IntegerMatrix mat, int max_row, int max_col) {
  int n = mat.nrow();
  LogicalVector out(n);
  for (int i = 0; i < n; ++i) {
    if (mat(i,0)>0 & mat(i,0)<=max_row & mat(i,1)>0 & mat(i,1)<=max_col) {
      out(i) = TRUE;
    } else {
      out(i) = FALSE;
    }
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector f2(IntegerMatrix mat, int max_col) {
  int n = mat.nrow();
  LogicalVector out(n);
  for (int i = 0; i < n; ++i) {
    if (mat(i,1)>0 & mat(i,1)<=max_col) {
      out(i) = TRUE;
    } else {
      out(i) = FALSE;
    }
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector f3(IntegerMatrix mat) {
  int n = mat.nrow();
  IntegerVector out(n);
  for (int i = 0; i < n; ++i) {
    if (mat(i,0)>=1) {
      out(i) = mat(i,2);
    } else {
      out(i) = mat(i,2)-(abs(mat(i,0))+1);
    }
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector f4(IntegerVector vec) {
  int n = vec.size();
  IntegerVector out(n);
  for (int i = 0; i < n; ++i) {
    if (vec[i]>1) {
      out[i] = vec[i];
    } else {
      out(i) = 1;
    }
  }
  return out;
}

// [[Rcpp::export]]
IntegerMatrix expand_gridC_2(IntegerVector vec1_1, IntegerVector vec1_2, IntegerVector vec2_1, IntegerVector vec2_2) {
  int n1 = vec1_1.size();
  int n2 = vec2_1.size();
  IntegerVector col1 = rep(vec1_1, n2);
  IntegerVector col2 = rep(vec1_2, n2);
  IntegerVector col3 = rep_each(vec2_1, n1);
  IntegerVector col4 = rep_each(vec2_2, n1);
  IntegerMatrix out = cbind(col1, col2, col3, col4);
  return out;
}
