#include <Rcpp.h>
using namespace Rcpp;
//' GLambda (Description)
//'
//' This function calculates genomic constant lambda from observed
//' p values. Values are calcuated as median of p-values' chi-test
//' statistic, divided by 0.4549364. This result indicate deflation
//' (< 1) or inflation (> 1) of the GWAS p values. If pvalues are
//' extremely inflated, then adjustments to model are necessary.
//' @title GLambda
//' @param x A numeric vector of P values
//' @return numeric value
//' @author Franjo, Mingjing, and Xiaoxiao
//' @examples
//' GLambda(ADHDMetaP$P)
//' @export
// [[Rcpp::export]]

double GLambda(NumericVector x) {
  NumericVector y = qchisq(x, 1, false, false);
  double p = 0.4549364;
  double m = median(y);
  double gc = m/p;
  return gc;
}
