#include "bits.h"
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

void propagateSDs(Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>> set_vars, std::shared_ptr<stumpRef> model, int max_vars) {
  if(model->n_vars == 1) return;
  for(int j = 0; j < model->n_vars; j++) {
    std::shared_ptr<stumpRef> model2(new stumpRef());

    std::shared_ptr<Rcpp::IntegerVector> vars(new Rcpp::IntegerVector(max_vars));
    std::copy( model->vars->begin(), model->vars->end(), vars->begin() ) ;

    for(int k = j; k < model->n_vars; k++)
      (*vars)[k] = (*vars)[k+1];
    (*vars)[model->n_vars-1] = Rcpp::IntegerVector::get_na();
    model2->vars = vars;
    auto & model2_orig = *(set_vars->find(model2));
    // Rcpp::Rcout << "IN: " << *(vars) << "\n";
    // Rcpp::Rcout << "OUT: " << *(model2_orig->vars) << "\n";
    for(int k = 0; k < model->min_sd->length(); k++) {
      if((*(model->min_sd))[k] < (*(model2_orig->min_sd))[k+1]) (*(model2_orig->min_sd))[k+1] = (*(model->min_sd))[k];
    }
    propagateSDs(set_vars, model2_orig, max_vars);
  }
}

// static void _finalizer(SEXP set_vars) {
//   if (R_ExternalPtrAddr(set_vars) == NULL)
//     return;
//   // Rprintf("finalizing\n");
//   std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>* ptr = (std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>*) R_ExternalPtrAddr(set_vars);
//   delete ptr;
//   R_ClearExternalPtr(set_vars);
// }

void initializeModel(Rcpp::NumericMatrix X, Rcpp::NumericVector neg_offset,
                     int max_vars,
                     std::shared_ptr<stumpRef> model2, int new_var,
                     std::shared_ptr<Rcpp::NumericVector> old_feature) {
  int N = X.nrow();
  std::shared_ptr<Rcpp::NumericVector> new_feature = std::shared_ptr<Rcpp::NumericVector>(new Rcpp::NumericVector(N));
  if(new_var > 0) {
    for(int i = 0; i < N; i++)
      (*new_feature)[i] = (*old_feature)[i] * X( i , new_var-1 );
  } else {
    for(int i = 0; i < N; i++)
      (*new_feature)[i] = (*old_feature)[i] * (neg_offset[-new_var-1] - X( i , -new_var-1 ));
  }

  double sx = Rcpp::var(*new_feature);
  sx = sqrt(sx * (N-1) / N);
  sx = std::max(sx, 1e-10);

  model2->sd = sx;
  model2->feature = new_feature;

  model2->min_sd = std::shared_ptr<Rcpp::NumericVector>(new Rcpp::NumericVector(max_vars-(model2->n_vars)+1, R_PosInf));
  (*(model2->min_sd))[0] = sx;
}

// [[Rcpp::export]]
SEXP initializeTerms(Rcpp::NumericMatrix X,
                     Rcpp::NumericVector neg_offset, int max_vars) {
  // std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>* set_vars = new std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>();
  Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>> set_vars( new std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>(), true );
  int p = X.ncol(); int N = X.nrow();
  std::shared_ptr<stumpRef> model(new stumpRef());
  model->n_vars = 0;

  std::shared_ptr<Rcpp::IntegerVector> vars(new Rcpp::IntegerVector(max_vars));
  std::fill(vars->begin(), vars->end(), Rcpp::IntegerVector::get_na());
  model->vars = vars;
  model->feature = std::shared_ptr<Rcpp::NumericVector>(new Rcpp::NumericVector(N, 1.0));
  model->min_sd = std::shared_ptr<Rcpp::NumericVector>(new Rcpp::NumericVector(max_vars+1, R_PosInf));
  set_vars->insert(model);

  std::queue<std::shared_ptr<stumpRef>> queue_vars;
  queue_vars.push(model);
  int order;
  int start_ind;
  std::shared_ptr<Rcpp::NumericVector> new_feature;

  std::list<std::shared_ptr<stumpRef>> list_prop;

  while(queue_vars.size() > 0) {
    model = queue_vars.front();
    queue_vars.pop();
    order = model->n_vars;

    if(order > 0) start_ind = (*(model->vars))[order-1]+1; else start_ind = -p;

    for(int j = start_ind; j <= p; j++) {
      if(j == 0) continue;
      // int ind = abs(j);
      // int sign = (j > 0)*2 - 1;

      if(j > 0) {
        bool skip = false;
        for(int k = 0; k < order; k++) {
          if(-(*(model->vars))[k] == j) skip = true;
        }
        if(skip) continue;
      }

      vars = std::shared_ptr<Rcpp::IntegerVector>(new Rcpp::IntegerVector(max_vars));
      std::copy( model->vars->begin(), model->vars->end(), vars->begin() ) ;
      // vars = Rcpp::clone(model->vars);
      (*vars)[order] = j;

      std::shared_ptr<stumpRef> model2(new stumpRef());
      model2->vars = vars;
      model2->n_vars = order+1;

      initializeModel(X, neg_offset, max_vars,
                      model2, j, model->feature);

      // new_feature = std::shared_ptr<Rcpp::NumericVector>(new Rcpp::NumericVector(N));
      // if(j > 0) {
      //   for(int i = 0; i < N; i++)
      //     (*new_feature)[i] = (*(model->feature))[i] * X( i , j-1 );
      // } else {
      //   for(int i = 0; i < N; i++)
      //     (*new_feature)[i] = (*(model->feature))[i] * (neg_offset[-j-1] - X( i , -j-1 ));
      // }
      //
      // double sx = Rcpp::var(*new_feature);
      // sx = sqrt(sx * (N-1) / N);
      // sx = std::max(sx, 1e-10);
      //
      // model2->sd = sx;
      // model2->feature = new_feature;
      //
      // model2->min_sd = std::shared_ptr<Rcpp::NumericVector>(new Rcpp::NumericVector(max_vars-order, R_PosInf));
      // (*(model2->min_sd))[0] = sx;

      if(order+1 < max_vars)
        queue_vars.push(model2);
      else
        list_prop.push_back(model2);
      set_vars->insert(model2);
    }
  }

  for (auto const& i : list_prop)
    propagateSDs(set_vars, i, max_vars);

  return set_vars;

  // SEXP ptr_set = PROTECT(R_MakeExternalPtr(set_vars, R_NilValue, R_NilValue));
  // R_RegisterCFinalizerEx(ptr_set, _finalizer, TRUE);

  // Rcpp::XPtr<std::list<stump*>> ptr_list( list_vars, true );
  // Rcpp::List L = Rcpp::List::create(ptr_set, ptr_list);
  // SEXP L = PROTECT(allocVector(VECSXP, 2));
  // SET_VECTOR_ELT(L, 0, ptr_set);
  // return L;
}

void descendantsInvestigated(std::shared_ptr<stumpRef> model,
                             std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>* set_vars,
                             int p) {
  int max_vars = model->vars->length();
  std::queue<std::shared_ptr<stumpRef>> queue_vars;
  queue_vars.push(model);
  while(queue_vars.size() > 0) {
    model = queue_vars.front();
    queue_vars.pop();

    int order = model->n_vars;
    for(int j = -p; j <= p; j++) {
      if(j == 0) continue;
      bool skip = false;
      for(int k = 0; k < order; k++) {
        if(abs((*(model->vars))[k]) == abs(j)) skip = true;
      }
      if(skip) continue;

      std::shared_ptr<Rcpp::IntegerVector> vars(new Rcpp::IntegerVector(max_vars));
      std::copy( model->vars->begin(), model->vars->end(), vars->begin() ) ;
      (*vars)[order] = j;
      std::sort(&(*vars)[0], &(*vars)[order+1]);

      std::shared_ptr<stumpRef> model2(new stumpRef());
      model2->vars = vars;
      model2->n_vars = order+1;

      auto & model2_orig = *(set_vars->find(model2));
      model2_orig->investigated = true;

      if(order+1 < max_vars) queue_vars.push(model2_orig);
    }
  }
}

std::vector<double> evalModel3(Rcpp::NumericMatrix X, Rcpp::NumericVector y,
                               const stumpRef& model,
                               int max_vars, double gamma) {
  std::shared_ptr<Rcpp::NumericVector> new_feature = model.feature;
  // ret = {score, positive score, negative score}
  std::vector<double> ret(3, 0.0);
  for(int i = 0; i < X.nrow(); i++) {
    if(y[i] > 0) ret[1] += (*new_feature)[i] * y[i];
    else if(y[i] < 0) ret[2] -= (*new_feature)[i] * y[i];
  }
  ret[0] = fabs(ret[1] - ret[2])/X.nrow();
  ret[0] /= model.sd;
  // ret[0] = ret[0] * ret[0]; //////

  ret[2] = std::max(ret[1], ret[2])/X.nrow();

  if(model.n_vars < max_vars) {
    double upper = R_NegInf;
    for(int i = 1; i <= max_vars - model.n_vars; i++) {
      // double new_denom = model.sd * (*(model.min_sd))[i];
      double new_denom = (*(model.min_sd))[i];
      double tmp_upper = ret[2] / new_denom;
      // tmp_upper = tmp_upper * tmp_upper; //////
      tmp_upper -= gamma * (model.n_vars + i);
      if(tmp_upper > upper) upper = tmp_upper;
    }
    ret[1] = upper;
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::List completeSearch(Rcpp::NumericMatrix X, Rcpp::NumericVector y,
                          Rcpp::NumericVector neg_offset,
                          int max_vars, double gamma,
                          SEXP set_vars_R, bool reuse_terms, int max_iter,
                          bool force_model, bool adjust_shady_int = true) {
  bool complete = max_iter <= 0;

  std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>* set_vars;
  // Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>> set_vars( new std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>(), true );
  if(Rf_isNull(set_vars_R))
    set_vars = new std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>();
  else
    set_vars = Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>>(set_vars_R);

  int p = X.ncol(); int N = X.nrow();
  std::shared_ptr<stumpRef> model(new stumpRef());
  model->n_vars = 0;

  std::shared_ptr<Rcpp::IntegerVector> vars(new Rcpp::IntegerVector(max_vars));
  std::fill(vars->begin(), vars->end(), Rcpp::IntegerVector::get_na());
  model->vars = vars;
  model->upper_bound = R_PosInf;
  model->score = 0.0;
  model->numerator = R_PosInf;

  model->feature = std::shared_ptr<Rcpp::NumericVector>(new Rcpp::NumericVector(N, 1.0));
  model->min_sd = std::shared_ptr<Rcpp::NumericVector>(new Rcpp::NumericVector(max_vars+1, R_PosInf));

  std::multiset<std::shared_ptr<stumpRef>,cmp_scoresRef> set_scores;
  set_scores.insert(model);
  int order;
  int start_ind;

  double best_score = 0.0;
  std::shared_ptr<stumpRef> best_model = model;
  int evaluated_terms = 0;
  bool found_better_model = false;

  /* double y_sd = 0;
  for(int i = 0; i < N; i++)
    y_sd += y[i] * y[i];
  y_sd = sqrt(y_sd/N); */
  // double y_sd = 1.0; //////

  for (const auto& elem: *set_vars)
    elem->investigated = false;

  while(set_scores.size() > 0) {
    model = *set_scores.rbegin();
    set_scores.erase(std::prev(set_scores.end()));
    order = model->n_vars;

    if(complete && model->upper_bound <= best_score) {
      // Rcpp::Rcout << "\r" << model->upper_bound << " | " << best_score;
      descendantsInvestigated(model, set_vars, p);
      continue;
    }

    if(max_iter > 0 && evaluated_terms == max_iter) break;

    if(complete && order > 0) start_ind = (*(model->vars))[order-1]+1; else start_ind = -p;

    // start_ind = -p; // DEBUG
    for(int j = start_ind; j <= p; j++) {
      if(max_iter > 0 && evaluated_terms == max_iter) break;
      if(j == 0) continue;
      // if(j > 0) {
      //   bool skip = false;
      //   for(int k = 0; k < order; k++) {
      //     if(-(*(model->vars))[k] == j) skip = true;
      //   }
      //   if(skip) continue;
      // }
      bool skip = false;
      for(int k = 0; k < order; k++) {
        if(abs((*(model->vars))[k]) == abs(j)) skip = true;
      }
      if(skip) continue;

      vars = std::shared_ptr<Rcpp::IntegerVector>(new Rcpp::IntegerVector(max_vars));
      std::copy( model->vars->begin(), model->vars->end(), vars->begin() ) ;
      // vars = Rcpp::clone(model->vars);
      (*vars)[order] = j;
      if(!complete) std::sort(&(*vars)[0], &(*vars)[order+1]);
      // std::sort(&(*vars)[0], &(*vars)[order+1]); // DEBUG

      std::shared_ptr<stumpRef> model2(new stumpRef());
      model2->vars = vars;
      model2->n_vars = order+1;
      model2->investigated = false;

      // auto & model2_orig = *(set_vars->find(model2));
      auto model2_iterator = set_vars->find(model2);
      std::shared_ptr<stumpRef> model2_orig;
      if(model2_iterator == set_vars->end()) {
        initializeModel(X, neg_offset, max_vars,
                        model2, j, model->feature);
        set_vars->insert(model2);
        model2_orig = model2;
      } else {
        model2_orig = *model2_iterator;
      }
      if(model2_orig->investigated) continue;

      // Rcpp::Rcout << model2_orig.sd << "\n";

      double tmp_upper = model->numerator/model2_orig->sd - gamma * model2_orig->n_vars;
      // double tmp_upper = pow(model->numerator/(model2_orig->sd * y_sd), 2) - gamma * model2_orig->n_vars; //////
      if(complete && tmp_upper <= best_score) {
        model2_orig->score = tmp_upper;
        model2_orig->upper_bound = model->upper_bound;
        model2_orig->numerator = model->numerator;
      } else {
        std::vector<double> res = evalModel3(X, y, *model2_orig, max_vars, gamma);
        model2_orig->score = res[0] - gamma * model2_orig->n_vars;
        model2_orig->upper_bound = res[1];
        model2_orig->numerator = res[2];
        evaluated_terms++;
      }
      model2_orig->investigated = true;

      if(model2_orig->score > best_score || (!found_better_model && force_model)) {
        best_score = model2_orig->score;
        best_model = model2_orig;
        found_better_model = true;
      }

      if(order+1 < max_vars) set_scores.insert(model2_orig);
    }
  }

  Rcpp::List best_result;

  Rcpp::List null_mod = customLm(Rcpp::NumericMatrix(N, 1, Rcpp::NumericVector(N, 1.0).begin()), y);
  Rcpp::NumericVector null_preds = Rcpp::as<Rcpp::NumericVector>(null_mod["preds"]);
  best_result["score"] = Rcpp::as<double>(null_mod["deviance"]);
  best_result["model"] = null_mod; best_result["preds"] = null_preds;
  vars = std::shared_ptr<Rcpp::IntegerVector>(new Rcpp::IntegerVector(max_vars));
  // vars = Rcpp::IntegerVector(max_vars);
  std::fill(vars->begin(), vars->end(), Rcpp::IntegerVector::get_na());
  vars->attr("dim") = Rcpp::IntegerVector::create(1, max_vars);
  best_result["vars"] = Rcpp::IntegerMatrix(*vars);

  vars = std::shared_ptr<Rcpp::IntegerVector>(new Rcpp::IntegerVector(max_vars));
  std::copy( best_model->vars->begin(), best_model->vars->end(), vars->begin() ) ;
  // vars = Rcpp::clone(best_model->vars);
  vars->attr("dim") = Rcpp::IntegerVector::create(1, max_vars);
  Rcpp::NumericMatrix Z;
  bool status = evalModel(best_result, X, y, Z, neg_offset, false, gamma, Rcpp::IntegerMatrix(*vars), 0, false, true);

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
  if(reuse_terms) {
    Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>> set_vars2 = (Rf_isNull(set_vars_R)) ? Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>>( set_vars, true ) : Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>>(set_vars_R);
    best_result["set_vars"] = set_vars2;
  } else {
    if(Rf_isNull(set_vars_R)) delete set_vars;
  }

  best_result.attr("class") = "intstump";
  return best_result;
}

// [[Rcpp::export]]
SEXP memoryTest() {
  Rcpp::XPtr<std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>> set_vars( new std::unordered_set<std::shared_ptr<stumpRef>,ArrayHasherRef,stumpPointEq>(), true );
  std::shared_ptr<Rcpp::IntegerVector> test(new Rcpp::IntegerVector(1e9));
  std::shared_ptr<stumpRef> model(new stumpRef());
  model->vars = test;
  set_vars->insert(model);
  return set_vars;
}


