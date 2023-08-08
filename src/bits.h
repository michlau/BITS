#ifndef BITS_H
#define BITS_H

#include <RcppArmadillo.h>
#include <queue>
#include <memory>

struct stumpRef
{
  std::shared_ptr<Rcpp::IntegerVector> vars;
  int n_vars;
  mutable double score;
  mutable double upper_bound;
  std::shared_ptr<Rcpp::NumericVector> feature;
  mutable bool investigated;

  double sd;
  std::shared_ptr<Rcpp::NumericVector> min_sd;
  mutable double numerator;

  // ~stumpRef() { Rcpp::Rcout << "Test destroyed.\n"; } // DEBUG
};

struct stumpPointEq {
  bool operator () ( std::shared_ptr<stumpRef> const lhs, std::shared_ptr<stumpRef> const rhs ) const {
    Rcpp::IntegerVector v1 = *(lhs->vars); Rcpp::IntegerVector v2 = *(rhs->vars);
    for(int i = 0; i < v1.length(); i++) {
      if(v1[i] != v2[i]) {
        return false;
      }
    }
    return true;
  }
};

struct cmp_scoresRef {
  template <typename T>
  bool operator() (const T& a, const T& b) const {
    return a->score < b->score;
  }
};

struct ArrayHasherRef {
  template <typename T>
  std::size_t operator()(const T& a) const {
    std::size_t h = 0;

    for (auto e : *(a->vars)) {
      h ^= std::hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

void propagateSDs(Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>> set_vars, std::shared_ptr<stumpRef> model, int max_vars);
void initializeModel(Rcpp::NumericMatrix X, Rcpp::NumericVector neg_offset, int max_vars, std::shared_ptr<stumpRef> model2, int new_var, std::shared_ptr<Rcpp::NumericVector> old_feature);
SEXP initializeTerms(Rcpp::NumericMatrix X, Rcpp::NumericVector neg_offset, int max_vars);
void descendantsInvestigated(std::shared_ptr<stumpRef> model, std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>* set_vars, int p);
std::vector<double> evalModel3(Rcpp::NumericMatrix X, Rcpp::NumericVector y, const stumpRef& model, int max_vars, double gamma);
Rcpp::List completeSearch(Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::NumericVector neg_offset, int max_vars, double gamma, SEXP set_vars_R, bool reuse_terms, int max_iter, bool force_model, bool adjust_shady_int);
SEXP memoryTest();

Rcpp::List customLm(Rcpp::NumericMatrix Xr, Rcpp::NumericVector yr);
Rcpp::NumericMatrix getInteractionFeature(Rcpp::NumericMatrix X, Rcpp::IntegerVector vars, Rcpp::NumericVector neg_offset);
Rcpp::NumericMatrix getAdditiveFeatures(Rcpp::NumericMatrix X, Rcpp::IntegerVector vars);
Rcpp::NumericMatrix cbind1 (Rcpp::NumericVector x, Rcpp::NumericVector y);
Rcpp::NumericMatrix matrixVectorMult(Rcpp::NumericMatrix Z, Rcpp::NumericVector y);
Rcpp::NumericMatrix matrixMatrixMult(Rcpp::NumericMatrix Z, Rcpp::NumericMatrix Y);
bool evalModel(const Rcpp::List & result, Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::NumericVector neg_offset, double gamma, Rcpp::IntegerMatrix vars, int n_vars_before, bool found_better_model, bool force_model);
bool evalAdditiveModels(const Rcpp::List & result, Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::NumericVector neg_offset, double gamma, Rcpp::IntegerMatrix vars);
double obj_fun_rcpp(double& x, Rcpp::NumericVector& y, Rcpp::NumericVector& oldEstimates, Rcpp::NumericVector& evaluatedWeakLearners);
double getRho(Rcpp::NumericVector y, Rcpp::NumericVector oldEstimates, Rcpp::NumericVector evaluatedWeakLearners, bool y_bin);
Rcpp::IntegerMatrix rbind1 (Rcpp::IntegerMatrix x, Rcpp::IntegerMatrix y);

#endif

