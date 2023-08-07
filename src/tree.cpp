#include "bits.h"
#include <queue>
#include <memory>

struct stump
{
  Rcpp::IntegerVector vars;
  int n_vars;
  mutable double score;
  mutable double upper_bound;
  Rcpp::NumericVector feature;
  bool investigated;

  double sd;
  mutable Rcpp::NumericVector min_sd;
  mutable double numerator;

  bool operator==(const stump& b) const
  {
    bool ret = true;
    Rcpp::IntegerVector v1 = this->vars; Rcpp::IntegerVector v2 = b.vars;
    for(int i = 0; i < v1.length(); i++) {
      if(v1[i] != v2[i]) {
        ret = false; break;
      }
    }
    return ret;
  }
};

struct cmp_scores {
  template <typename T>
  bool operator() (const T& a, const T& b) const {
    return a.score < b.score;
  }
};

struct cmp_vars {
  bool operator() (const stump& a, const stump& b) const {
    bool ret = false;
    for(int i = 0; i < a.vars.length(); i++) {
      if(a.vars[i] < b.vars[i]) {
        ret = true; break;
      } else if(a.vars[i] > b.vars[i]) {
        ret = false; break;
      }
    }
    return ret;
  }
};

struct ArrayHasher {
  template <typename T>
  std::size_t operator()(const T& a) const {
    std::size_t h = 0;

    for (auto e : a.vars) {
      h ^= std::hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

Rcpp::NumericVector g_sds;
int g_max_vars, g_n_vars;
double g_gamma;
int g_optimize;
bool g_independent;

std::vector<double> evalModel2(Rcpp::NumericMatrix X, Rcpp::NumericVector y,
                               Rcpp::NumericVector old_feature, Rcpp::NumericVector & new_feature,
                               Rcpp::NumericVector neg_offset, int new_var,
                               Rcpp::IntegerVector vars) {
  std::copy( old_feature.begin(), old_feature.end(), new_feature.begin() ) ;
  for(int i = 0; i < X.nrow(); i++) {
    if(new_var > 0)
      new_feature[i] *= X( i , new_var-1 );
    else
      new_feature[i] *= neg_offset[-new_var-1] - X( i , -new_var-1 );
  }

  // ret = {score, positive score, negative score}
  std::vector<double> ret(3, 0);
  for(int i = 0; i < X.nrow(); i++) {
    if(y[i] > 0) ret[1] += new_feature[i] * y[i];
    else if(y[i] < 0) ret[2] -= new_feature[i] * y[i];
  }
  ret[0] = fabs(ret[1] - ret[2]);

  if(g_optimize) {
    ret[0] /= X.nrow();
    double upper = std::max(ret[1], ret[2])/X.nrow() - g_gamma * (g_n_vars + 1);
    ret[1] = upper; ret[2] = upper;
  } else {
    double sx = Rcpp::var(new_feature);
    sx = sqrt(sx * (X.nrow()-1) * X.nrow());
    if(fabs(sx) <= 1e-10) ret[0] = R_NegInf; else ret[0] /= sx;

    if(g_independent) {
      double sx2 = 0;
      for(int i = 0; i < X.nrow(); i++)
        sx2 += new_feature[i] * new_feature[i];
      sx2 = sqrt(sx2 * X.nrow());

      // Good
      Rcpp::NumericVector sds_sub(X.ncol() - g_n_vars);
      int ii = 0;
      for(int i = 0; i < X.ncol(); i++) {
        bool found = false;
        for(int j = 0; j < g_n_vars; j++) {
          if(vars[j] == i+1 || vars[j] == -(i+1)) {
            found = true;
            break;
          }
        }
        if(!found) sds_sub[ii++] = g_sds[i];
      }
      std::sort(sds_sub.begin(), sds_sub.end());

      double upper = R_NegInf;
      double min_sd = 1;
      for(int i = 1; i <= g_max_vars - g_n_vars; i++) {
        min_sd *= sds_sub[i-1];
        double new_denom = sx2 * min_sd;
        double tmp_upper = std::max(ret[1], ret[2]) / new_denom;
        tmp_upper -= g_gamma * (g_n_vars + i);
        if(tmp_upper > upper) upper = tmp_upper;
      }
      ret[1] = upper; ret[2] = upper;
    }

    /*double min_sd = R_PosInf;
    for(int i = 0; i < X.ncol(); i++) {
      bool found = false;
      for(int j = 0; j < g_n_vars; j++) {
        if(vars[j] == i+1 || vars[j] == -(i+1)) {
          found = true;
          break;
        }
      }
      if(!found && g_sds[i] < min_sd) min_sd = g_sds[i];
    }
    double upper = R_NegInf;
    for(int i = 1; i <= g_max_vars - g_n_vars; i++) {
      double new_denom = sx2 * pow(min_sd, i);
      double tmp_upper = std::max(ret[1], ret[2]) / new_denom;
      tmp_upper -= g_gamma * (g_n_vars + i);
      if(tmp_upper > upper) upper = tmp_upper;
    }
    ret[1] = upper; ret[2] = upper;*/

    /*double tmp_sum1 = 0; double tmp_sum2 = 0;
     for(int i = 0; i < X.nrow(); i++) {
     tmp_sum1 += new_feature[i] * new_feature[i] * y[i] * y[i];
     tmp_sum2 += new_feature[i] * new_feature[i];
     }
     double upper = sqrt(tmp_sum1 / tmp_sum2) - g_gamma * (g_n_vars + 1);
     ret[1] = upper; ret[2] = upper;*/
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::List treeFit(Rcpp::NumericMatrix X, Rcpp::NumericVector y_orig,
                   Rcpp::NumericVector neg_offset,
                   int max_vars, double gamma,
                   int optimize, bool independent, Rcpp::NumericVector sds,
                   bool force_model, bool adjust_shady_int = true,
                   int max_iter = 10000) {
  // Center y
  // Rcpp::NumericVector y = Rcpp::clone(y_orig);
  // y = y - Rcpp::mean(y);
  Rcpp::NumericVector y = y_orig;

  g_sds = sds;
  g_max_vars = max_vars;
  g_gamma = gamma;
  g_optimize = optimize;
  g_independent = independent;

  double best_score = 0; // R_NegInf;
  stump best_model;

  std::multiset<stump,cmp_scores> set_scores;
  // std::set<stump,cmp_vars> set_vars;
  std::unordered_set<stump,ArrayHasher> set_vars;

  int p = X.ncol(); int N = X.nrow();
  Rcpp::NumericMatrix Z;
  Rcpp::IntegerVector vars;
  Rcpp::List best_result;
  Rcpp::NumericVector old_feature(N, 1.0);
  Rcpp::NumericVector new_feature;
  // best_result["score"] = 0; // R_NegInf;
  stump model;
  bool found_better_model = false;

  vars = Rcpp::IntegerVector(max_vars);
  std::fill(vars.begin(), vars.end(), Rcpp::IntegerVector::get_na());
  best_model.feature = old_feature; best_model.vars = vars; best_model.upper_bound = R_PosInf;
  best_model.investigated = false; best_model.n_vars = 0; best_model.score = 0;
  // set_vars.insert(best_model); set_scores.insert(best_model);

  Rcpp::List null_mod = customLm(Rcpp::NumericMatrix(N, 1, Rcpp::NumericVector(N, 1.0).begin()), y);
  Rcpp::NumericVector null_preds = Rcpp::as<Rcpp::NumericVector>(null_mod["preds"]);
  best_result["score"] = Rcpp::as<double>(null_mod["deviance"]);
  best_result["model"] = null_mod; best_result["preds"] = null_preds;
  vars.attr("dim") = Rcpp::IntegerVector::create(1, max_vars);
  best_result["vars"] = Rcpp::IntegerMatrix(vars);
  // best_result["score"] = R_PosInf;

  // Initialization: Evaluate and add all single features
  for(int j = 0; j < 2*p; j++) {
    vars = Rcpp::IntegerVector(max_vars);
    std::fill(vars.begin(), vars.end(), Rcpp::IntegerVector::get_na());
    int ind = j / 2 + 1;
    int sign = -(2 * (j % 2) - 1);
    vars[0] = sign * ind;
    int n_vars = 1;
    new_feature = Rcpp::NumericVector(N);
    g_n_vars = n_vars;
    std::vector<double> scores = evalModel2(X, y, old_feature, new_feature, neg_offset, sign * ind, vars);
    stump model2;
    model2.feature = new_feature; model2.score = scores[0] - gamma * n_vars; model2.vars = vars; model2.n_vars = n_vars;
    // (1+n_vars) due to the upper bound being valid for descendant nodes, i.e., starting from nodes with one more variable
    //// if(scores[1] > scores[2]) model2.upper_bound = scores[1] - gamma * (1+n_vars); else model2.upper_bound = scores[2] - gamma * (1+n_vars);
    model2.upper_bound = scores[1];
    model2.investigated = false;
    set_vars.insert(model2); set_scores.insert(model2);
    if(model2.score > best_score || (!found_better_model && force_model)) {
      best_score = model2.score; best_model = model2;
      found_better_model = true;
    }
  }

  int evaluated_terms = 2*p;

  while(set_scores.size() > 0) {
    model = *set_scores.rbegin();
    set_scores.erase(std::prev(set_scores.end()));

    if(model.n_vars == max_vars) continue;

    if((optimize || independent) && model.upper_bound <= best_score) continue; ////
    // Rcpp::Rcout << "Best Score: " << best_score << " | Upper Bound: " << model.upper_bound << "\n";

    if(max_iter > 0 && evaluated_terms == max_iter) break;

    // Add all (possible) variables to the interaction term
    for(int j = 0; j < 2*p; j++) {
      if(max_iter > 0 && evaluated_terms == max_iter) break;
      int ind = j / 2 + 1;
      int sign = -(2 * (j % 2) - 1);

      bool skip = false;
      // Plausibility check
      for(int k = 0; k < model.n_vars; k++) {
        if(abs(model.vars[k]) == ind) skip = true;
      }
      if(skip) continue;

      vars = Rcpp::clone(model.vars);
      vars[model.n_vars] = sign * ind;

      // ToDo: Utilize model.feature, scale input vars to [0,1]?

      int n_vars = model.n_vars + 1;

      std::sort(&vars[0], &vars[n_vars]);

      stump model2; model2.vars = vars;
      if(set_vars.count(model2)) continue;

      new_feature = Rcpp::NumericVector(N);
      g_n_vars = n_vars;
      std::vector<double> scores = evalModel2(X, y, model.feature, new_feature, neg_offset, sign * ind, vars);
      model2.feature = new_feature; model2.score = scores[0] - gamma * n_vars; model2.n_vars = n_vars;
      // (1+n_vars) due to the upper bound being valid for descendant nodes, i.e., starting from nodes with one more variable
      //// if(scores[1] > scores[2]) model2.upper_bound = scores[1] - gamma * (1+n_vars); else model2.upper_bound = scores[2] - gamma * (1+n_vars);
      model2.upper_bound = scores[1];
      model2.investigated = false;
      set_vars.insert(model2); set_scores.insert(model2);
      if(model2.score > best_score) {
        best_score = model2.score; best_model = model2;
      }
      // Rcpp::Rcout << "Vars: " << vars << "\n";
      evaluated_terms++;
    }
  }

  vars = best_model.vars;
  vars.attr("dim") = Rcpp::IntegerVector::create(1, max_vars);
  bool status = evalModel(best_result, X, y, Z, neg_offset, false, gamma, Rcpp::IntegerMatrix(vars), 0, false, true);

  if(status && adjust_shady_int) {
    Rcpp::IntegerMatrix vars_mat = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(best_result["vars"]));
    evalAdditiveModels(best_result, X, y, Z, neg_offset, false, gamma, vars_mat);
  }

  int possible_terms = 0;
  for(int i = 0; i < max_vars; i++) {
    possible_terms += Rf_choose(2*p, i+1);
    if(i > 0) possible_terms -= p * Rf_choose(2*(p-1), i+1-2); // Remove implausible combinations such as -1 & 1 & 3
  }
  best_result["possible_terms"] = possible_terms; best_result["evaluated_terms"] = evaluated_terms;
  best_result["corr"] = best_score;

  best_result.attr("class") = "intstump";
  return best_result;
}
