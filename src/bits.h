#ifndef BITS_H
#define BITS_H

#include <RcppArmadillo.h>

Rcpp::List customLm(Rcpp::NumericMatrix Xr, Rcpp::NumericVector yr);
Rcpp::NumericMatrix getInteractionFeature(Rcpp::NumericMatrix X, Rcpp::IntegerVector vars, Rcpp::NumericVector neg_offset);
Rcpp::NumericMatrix getAdditiveFeatures(Rcpp::NumericMatrix X, Rcpp::IntegerVector vars);
Rcpp::NumericMatrix cbind1 (Rcpp::NumericVector x, Rcpp::NumericVector y);
Rcpp::NumericMatrix matrixVectorMult(Rcpp::NumericMatrix Z, Rcpp::NumericVector y);
Rcpp::NumericMatrix matrixMatrixMult(Rcpp::NumericMatrix Z, Rcpp::NumericMatrix Y);
bool evalModel(const Rcpp::List & result, Rcpp::NumericMatrix X, Rcpp::NumericVector y,
               Rcpp::NumericMatrix Z, Rcpp::NumericVector neg_offset, bool use_Z, double gamma, Rcpp::IntegerMatrix vars, int n_vars_before,
               bool found_better_model, bool force_model);
bool evalAdditiveModels(const Rcpp::List & result, Rcpp::NumericMatrix X, Rcpp::NumericVector y,
                        Rcpp::NumericMatrix Z, Rcpp::NumericVector neg_offset, bool use_Z, double gamma, Rcpp::IntegerMatrix vars);
Rcpp::List greedyFit(Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::Nullable<Rcpp::NumericMatrix> Z_,
                     Rcpp::NumericVector neg_offset,
                     int max_vars, double gamma, bool force_model, bool adjust_shady_int,
                     bool modify_vars, bool remove_vars);
double obj_fun_rcpp(double& x, Rcpp::NumericVector& y, Rcpp::NumericVector& oldEstimates, Rcpp::NumericVector& evaluatedWeakLearners);
double getRho(Rcpp::NumericVector y, Rcpp::NumericVector oldEstimates, Rcpp::NumericVector evaluatedWeakLearners, bool y_bin);
Rcpp::IntegerMatrix rbind1 (Rcpp::IntegerMatrix x, Rcpp::IntegerMatrix y);
Rcpp::List boosting(Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::Nullable<Rcpp::NumericMatrix> Z_,
                    Rcpp::NumericVector neg_offset,
                    bool y_bin, int max_vars, double gamma,
                    int boosting_iter, double learning_rate,
                    bool adjust_shady_int);

#endif

