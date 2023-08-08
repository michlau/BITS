// [[Rcpp::depends(RcppArmadillo)]]
#include "bits.h"

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
  bool status = evalModel(best_result, X, y, neg_offset, gamma, Rcpp::IntegerMatrix(*vars), 0, false, true);

  if(status && adjust_shady_int) {
    Rcpp::IntegerMatrix vars_mat = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(best_result["vars"]));
    evalAdditiveModels(best_result, X, y, neg_offset, gamma, vars_mat);
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



// [[Rcpp::export]]
Rcpp::List customLm(Rcpp::NumericMatrix Xr, Rcpp::NumericVector yr) {
  arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
  arma::colvec y(yr.begin(), yr.size(), false);

  Rcpp::List lm;

  if(Xr.ncol() == 1) {
    double ymean = Rcpp::mean(yr);
    lm["preds"] = Rcpp::NumericVector(Xr.nrow(), ymean);
    lm["coef"] = Rcpp::NumericVector(1, ymean);
    lm["deviance"] = (double) (Rcpp::mean(Rcpp::pow(yr - ymean, 2)));
    lm["success"] = true;
    return lm;
  } else if (Xr.ncol() == 2) {
    int N = Xr.nrow();
    double ysum = 0;
    double xsum = 0; double xysum = 0;
    double xxsum = 0;
    for(int i = 0; i < N; i++) {
      ysum += yr[i];
      xsum += Xr(i, 1);
      xysum += Xr(i, 1) * yr[i];
      xxsum += Xr(i, 1) * Xr(i, 1);
    }
    double beta1 = (xysum - xsum*ysum/N)/(xxsum - xsum*xsum/N);
    double beta0 = (ysum - beta1*xsum)/N;
    Rcpp::NumericVector preds(N, beta0);
    for(int i = 0; i < N; i++) {
      preds[i] += beta1 * Xr(i, 1);
    }
    lm["preds"] = preds;
    Rcpp::NumericVector coef = Rcpp::NumericVector::create(beta0, beta1);
    lm["coef"] = coef;
    double dev = Rcpp::mean(Rcpp::pow(yr - preds, 2));
    lm["deviance"] = dev;
    lm["success"] = Rcpp::traits::is_finite<REALSXP>(beta0) && Rcpp::traits::is_finite<REALSXP>(beta1);
    return lm;
  }

  arma::colvec coef;
  bool success = arma::solve(coef, X, y, arma::solve_opts::no_approx);

  if(success) {
    Rcpp::NumericVector coef2 = Rcpp::NumericVector(coef.begin(), coef.end());
    arma::colvec preds = X*coef;
    Rcpp::NumericVector preds2 = Rcpp::NumericVector(preds.begin(), preds.end());
    double deviance = arma::accu(arma::pow(y - preds, 2)) / yr.size();

    lm["coef"] = coef2; lm["preds"] = preds2; lm["deviance"] = deviance; lm["success"] = true;
  } else {
    lm["success"] = false;
  }

  return lm;
}

Rcpp::NumericMatrix getInteractionFeature(Rcpp::NumericMatrix X, Rcpp::IntegerVector vars, Rcpp::NumericVector neg_offset) {
  int n_vars;
  for(n_vars = 0; n_vars < vars.length(); n_vars++) {
    if(Rcpp::IntegerVector::is_na(vars[n_vars])) break;
  }
  Rcpp::NumericMatrix X_tmp(X.nrow(), n_vars);
  for(int i = 0; i < n_vars; i++) {
    // Rcpp::Rcout << "Column index : " << abs(vars[i])-1 << "\n";
    if(vars[i] > 0)
      X_tmp( Rcpp::_ , i ) = X( Rcpp::_ , vars[i]-1 );
    else
      X_tmp( Rcpp::_ , i ) = neg_offset[-vars[i]-1] - X( Rcpp::_ , -vars[i]-1 );
  }
  Rcpp::NumericVector ret(X.nrow());
  ret.fill(1.0);
  for(int i = 0; i < X.nrow(); i++) {
    for(int j = 0; j < n_vars; j++)
      ret[i] *= X_tmp( i, j );
  }
  Rcpp::NumericMatrix ret2(X.nrow(), 1, ret.begin());
  return ret2;
}

Rcpp::NumericMatrix getAdditiveFeatures(Rcpp::NumericMatrix X, Rcpp::IntegerVector vars) {
  int n_vars;
  for(n_vars = 0; n_vars < vars.length(); n_vars++) {
    if(Rcpp::IntegerVector::is_na(vars[n_vars])) break;
  }
  Rcpp::NumericMatrix X_tmp(X.nrow(), n_vars);
  for(int i = 0; i < n_vars; i++) {
    X_tmp( Rcpp::_ , i ) = X( Rcpp::_ , abs(vars[i])-1 );
  }
  return X_tmp;
}

Rcpp::NumericMatrix cbind1 (Rcpp::NumericVector x, Rcpp::NumericVector y) {
  Rcpp::NumericMatrix out(x.size(), 2);
  out( Rcpp::_ , 0 ) = x; out( Rcpp::_ , 1 ) = y;
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix matrixVectorMult(Rcpp::NumericMatrix Z, Rcpp::NumericVector y) {
  Rcpp::NumericMatrix m(Z.nrow(), Z.ncol());
  for(int i = 0; i < Z.nrow(); i++)
    m(i, Rcpp::_) = Z(i, Rcpp::_) * y[i];
  return m;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix matrixMatrixMult(Rcpp::NumericMatrix Z, Rcpp::NumericMatrix Y) {
  Rcpp::NumericMatrix m(Z.nrow(), Z.ncol() * Y.ncol());
  for(int i = 0; i < Z.nrow(); i++) {
    for(int j = 0; j < Z.ncol(); j++) {
      for(int k = 0; k < Y.ncol(); k++) {
        m(i, j*Y.ncol()+k) = Z(i, j) * Y(i, k);
      }
    }
  }
  return m;
}

bool evalModel(const Rcpp::List & result, Rcpp::NumericMatrix X, Rcpp::NumericVector y,
               Rcpp::NumericVector neg_offset, double gamma, Rcpp::IntegerMatrix vars, int n_vars_before,
               bool found_better_model, bool force_model) {
  // Rcpp::NumericVector interaction_feature;
  int n_terms = vars.nrow();
  Rcpp::NumericMatrix interaction_feature(y.size(), n_terms);
  for(int j = 0; j < n_terms; j++)
    interaction_feature(Rcpp::_, j) = getInteractionFeature(X, vars(j, Rcpp::_), neg_offset);
  Rcpp::NumericVector ones(y.size(), 1.0);
  ones.attr("dim") = Rcpp::IntegerVector::create(y.size(), 1);
  // Rcpp::NumericMatrix dm = cbind1(ones, interaction_feature);
  Rcpp::NumericMatrix dm = Rcpp::cbind(Rcpp::NumericMatrix(ones), interaction_feature);
  Rcpp::List mod = customLm(dm, y);
  Rcpp::NumericVector preds;
  double score;
  double best_score = Rcpp::as<double>(result["score"]);
  bool status = false;
  int n_vars = 0;
  for(int j = 0; j < vars.length(); j++) {
    if(!Rcpp::IntegerVector::is_na(vars[j])) n_vars++;
  }
  if(Rcpp::as<bool>(mod["success"])) {
    preds = Rcpp::as<Rcpp::NumericVector>(mod["preds"]);
    score = Rcpp::as<double>(mod["deviance"]) + gamma * n_vars;
    if(score < best_score || (n_vars_before == 0 && !found_better_model && force_model)) {
      result["score"] = score; result["vars"] = Rcpp::clone(vars); result["model"] = Rcpp::clone(mod); result["preds"] = Rcpp::clone(preds);
      status = true;
    }
  }
  return status;
}

bool evalAdditiveModels(const Rcpp::List & result, Rcpp::NumericMatrix X, Rcpp::NumericVector y,
                        Rcpp::NumericVector neg_offset, double gamma, Rcpp::IntegerMatrix vars) {
  int n_vars;
  for(n_vars = 0; n_vars < vars.length(); n_vars++) {
    if(Rcpp::IntegerVector::is_na(vars[n_vars])) break;
  }
  if(n_vars < 2)
    return false;
  Rcpp::IntegerMatrix vars_mat(n_vars, 1, vars.begin());
  bool status = false;
  if(evalModel(result, X, y, neg_offset, gamma, vars_mat, 1, true, false))
    status = true;

  if(n_vars > 2) {
    /* Split one variable off and test the rest as interaction
       This guarantees that every partition for up to 3 variables is checked
       For more than 3 variables, an algorithm should be implemented that
       generates every possible partition */
    vars_mat = Rcpp::IntegerMatrix(2, n_vars-1);
    for(int j = 0; j < n_vars; j++) {
      std::fill(vars_mat.begin(), vars_mat.end(), Rcpp::IntegerVector::get_na());
      vars_mat(0, 0) = vars[j];
      int k2 = 0;
      for(int k = 0; k < n_vars-1; k++) {
        if(j == k) k2++;
        vars_mat(1, k) = vars[k2];
        k2++;
      }
      if(evalModel(result, X, y, neg_offset, gamma, vars_mat, 1, true, false))
        status = true;
    }
  }
  return status;
}

double obj_fun_rcpp(double& x, Rcpp::NumericVector& y, Rcpp::NumericVector& oldEstimates, Rcpp::NumericVector& evaluatedWeakLearners) {
  return Rcpp::sum(Rcpp::log(1+Rcpp::exp(2*(oldEstimates + x * evaluatedWeakLearners))) - 2 * y * (oldEstimates + x*evaluatedWeakLearners));
}

double getRho(Rcpp::NumericVector y, Rcpp::NumericVector oldEstimates, Rcpp::NumericVector evaluatedWeakLearners, bool y_bin) {
  if(y_bin) {
    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];
    double init_val = 1;

    // Call the optim function from R in C++
    Rcpp::List opt_results = optim(Rcpp::_["par"]    = init_val,
                                   // Make sure this function is not exported!
                                   Rcpp::_["fn"]     = Rcpp::InternalFunction(&obj_fun_rcpp),
                                   Rcpp::_["method"] = "Brent",
                                   // Pass in the other parameters as everything
                                   // is scoped environmentally
                                   Rcpp::_["y"] = y,
                                   Rcpp::_["oldEstimates"] = oldEstimates,
                                   Rcpp::_["evaluatedWeakLearners"] = evaluatedWeakLearners);
    double rho = Rcpp::as<Rcpp::NumericVector>(opt_results[0])[0];
    return rho;
  } else {
    return Rcpp::sum((y - oldEstimates) * evaluatedWeakLearners)/Rcpp::sum(Rcpp::pow(evaluatedWeakLearners, 2));
  }
}

Rcpp::IntegerMatrix rbind1 (Rcpp::IntegerMatrix x, Rcpp::IntegerMatrix y) {
  Rcpp::IntegerMatrix out (x.nrow()+y.nrow(), x.ncol());
  for(int i = 0; i < x.nrow(); i++)
    out(i, Rcpp::_) = x(i, Rcpp::_);
  for(int i = 0; i < y.nrow(); i++)
    out(i+x.nrow(), Rcpp::_) = y(i, Rcpp::_);
  return out;
}


