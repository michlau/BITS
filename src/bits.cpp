// [[Rcpp::depends(RcppArmadillo)]]
//// #include <RcppArmadillo.h>
#include "bits.h"
// using namespace Rcpp;

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
               Rcpp::NumericMatrix Z, Rcpp::NumericVector neg_offset, bool use_Z, double gamma, Rcpp::IntegerMatrix vars, int n_vars_before,
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
  if(use_Z) {
    // Rcpp::NumericVector temp_vec = Z * interaction_feature;
    // temp_vec.attr("dim") = Rcpp::Dimension(Z.nrow(), Z.ncol());
    // dm = Rcpp::cbind(dm, Z, Rcpp::as<Rcpp::NumericMatrix>(temp_vec));
    //// dm = Rcpp::cbind(dm, Z, matrixVectorMult(Z, interaction_feature));
    dm = Rcpp::cbind(dm, Z, matrixMatrixMult(Z, interaction_feature));
  }
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
                        Rcpp::NumericMatrix Z, Rcpp::NumericVector neg_offset, bool use_Z, double gamma, Rcpp::IntegerMatrix vars) {
  int n_vars;
  for(n_vars = 0; n_vars < vars.length(); n_vars++) {
    if(Rcpp::IntegerVector::is_na(vars[n_vars])) break;
  }
  if(n_vars < 2)
    return false;
  Rcpp::IntegerMatrix vars_mat(n_vars, 1, vars.begin());
  bool status = false;
  if(evalModel(result, X, y, Z, neg_offset, use_Z, gamma, vars_mat, 1, true, false))
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
      if(evalModel(result, X, y, Z, neg_offset, use_Z, gamma, vars_mat, 1, true, false))
        status = true;
    }
  }
  return status;
}

// [[Rcpp::export]]
Rcpp::List greedyFit(Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::Nullable<Rcpp::NumericMatrix> Z_,
                     Rcpp::NumericVector neg_offset,
                     int max_vars, double gamma, bool force_model, bool adjust_shady_int = true,
                     bool modify_vars = false, bool remove_vars = false) {
  Rcpp::List best_result;
  int p = X.ncol(); int N = X.nrow();
  bool use_Z = Z_.isNotNull();
  Rcpp::NumericMatrix Z;
  if(use_Z) Z = Rcpp::NumericMatrix(Z_);
  Rcpp::IntegerMatrix vars(1, max_vars);
  std::fill(vars.begin(), vars.end(), Rcpp::IntegerVector::get_na());
  Rcpp::IntegerMatrix best_vars = Rcpp::clone(vars);
  // bool y_bin = !any(!in(y, NumericVector::create(0, 1)));
  Rcpp::List best_mod = customLm(Rcpp::NumericMatrix(N, 1, Rcpp::NumericVector(N, 1.0).begin()), y);
  // class(best.mod) <- "glm"
  Rcpp::NumericVector best_preds = Rcpp::as<Rcpp::NumericVector>(best_mod["preds"]);
  double best_score = Rcpp::as<double>(best_mod["deviance"]);

  bool found_better_model;
  Rcpp::NumericMatrix X_tmp, dm;
  Rcpp::NumericVector interaction_feature;
  Rcpp::List mod;
  Rcpp::NumericVector preds;
  Rcpp::NumericVector ones(N, 1.0);

  // Only addition of variables allowed
  // for(int i = 0; i < max_vars; i++) {
  //   found_better_model = false;
  //   vars = Rcpp::clone(best_vars);
  //   for(int j = 1; j <= 2*p; j++) {
  //     if(j <= p) k = j; else k = -(j-p);
  //     vars[i] = k;
  //     interaction_feature = getInteractionFeature(X, vars);
  //     // if(!use_Z)
  //     //   dm = cbind(ones, interaction_feature);
  //     // else
  //     //   dm = cbind(ones, interaction_feature, Z, Z * interaction_feature);
  //     dm = cbind1(ones, interaction_feature);
  //     // dm = Rcpp::NumericMatrix(N, 1, 1.0);
  //     mod = customLm(dm, y);
  //     if(Rcpp::as<bool>(mod["success"]) == false) continue;
  //     preds = Rcpp::as<Rcpp::NumericVector>(mod["preds"]);
  //     score = Rcpp::as<double>(mod["deviance"]) + gamma * (i+1);
  //     if(score <= best_score || (i == 0 && !found_better_model && force_model)) {
  //       best_score = score; best_vars = Rcpp::clone(vars); best_mod = Rcpp::clone(mod); best_preds = Rcpp::clone(preds);
  //       found_better_model = true;
  //     }
  //   }
  //
  //   if(!found_better_model) break;
  // }

  best_result["score"] = best_score; best_result["model"] = best_mod; best_result["preds"] = best_preds; best_result["vars"] = best_vars;
  std::vector<int> unused_vars(p, 0);
  std::vector<int> available_vars(2*p, 0);
  int n_free_vars;
  found_better_model = true;
  int n_vars, j, var_buffer;
  bool eval_status;
  while(found_better_model) {
    found_better_model = false;
    vars = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(best_result["vars"]));
    for(n_vars = 0; n_vars < vars.length(); n_vars++) {
      if(Rcpp::IntegerVector::is_na(vars[n_vars])) break;
    }
    std::fill(unused_vars.begin(), unused_vars.end(), 0);
    for(j = 0; j < n_vars; j++)
      unused_vars[abs(vars[j]) - 1] = 1;
    n_free_vars = 0;
    for(j = 0; j < p; j++) {
      if(!unused_vars[j]) {
        available_vars[n_free_vars] = j+1;
        if(neg_offset[j] != 0) {
          available_vars[n_free_vars+1] = -(j+1);
          n_free_vars += 2;
        } else {
          n_free_vars++;
        }
      }
    }

    // Add variable
    if (n_vars < max_vars) {
      for(j = 0; j < n_free_vars; j++) {
        vars[n_vars] = available_vars[j];
        eval_status = evalModel(best_result, X, y, Z, neg_offset, use_Z, gamma, vars, n_vars,
                                found_better_model, force_model);
        if(eval_status) found_better_model = true;
      }
      vars[n_vars] = NA_INTEGER;
    }

    // Modify variable
    if(modify_vars) {
      for(j = 0; j < n_vars; j++) {
        var_buffer = vars[j];
        for(int l = 0; l < n_free_vars+1; l++) {
          if(l == n_free_vars)
            vars[j] = -vars[j];
          else
            vars[j] = available_vars[l];
          eval_status = evalModel(best_result, X, y, Z, neg_offset, use_Z, gamma, vars, n_vars,
                                  found_better_model, force_model);
          if(eval_status) found_better_model = true;
        }
        vars[j] = var_buffer;
      }
    }

    // Remove variable
    if (n_vars > 1 && remove_vars) {
      for(j = 0; j < n_vars; j++) {
        var_buffer = vars[j];
        vars[j] = vars[n_vars-1];
        vars[n_vars-1] = NA_INTEGER;
        eval_status = evalModel(best_result, X, y, Z, neg_offset, use_Z, gamma, vars, n_vars,
                                found_better_model, force_model);
        if(eval_status) found_better_model = true;
        vars[n_vars-1] = vars[j];
        vars[j] = var_buffer;
      }
    }

    // Rcpp::Rcout << "Vars : " << vars << " Score : " << Rcpp::as<double>(best_result["score"]) << std::endl;
  }
  // i = 26, cont_no_E, gamma = 0.001

  /* NEW */
  if(adjust_shady_int) {
    vars = Rcpp::clone(Rcpp::as<Rcpp::IntegerMatrix>(best_result["vars"]));
    eval_status = evalAdditiveModels(best_result, X, y, Z, neg_offset, use_Z, gamma, vars);
  }

  best_result.attr("class") = "intstump";
  return best_result;
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

// [[Rcpp::export]]
Rcpp::List boosting(Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::Nullable<Rcpp::NumericMatrix> Z_,
                    Rcpp::NumericVector neg_offset,
                    bool y_bin, int max_vars, double gamma,
                    int boosting_iter, double learning_rate,
                    bool adjust_shady_int = true) {

  int N = y.size();
  double priorProb = Rcpp::mean(y);
  double initialModel = priorProb;
  if(y_bin)
    initialModel = 0.5 * log(priorProb/(1-priorProb));

  Rcpp::IntegerMatrix disj(0, max_vars);
  Rcpp::NumericVector currentEstimates(N, initialModel);
  Rcpp::NumericVector probs; Rcpp::NumericVector gradient;
  Rcpp::List model;
  Rcpp::IntegerMatrix vars;
  int n_vars;
  double currentRho;
  for(int i = 0; i < boosting_iter; i++) {
    if(y_bin) {
      probs = 1/(1+Rcpp::exp(-2 * currentEstimates));
      gradient = 2*probs - 2*y;
    } else {
      gradient = currentEstimates - y;
    }
    model = greedyFit(X, -gradient, Z_, neg_offset, max_vars, gamma,
                      (i == 0), adjust_shady_int);
    vars = Rcpp::as<Rcpp::IntegerMatrix>(model["vars"]);
    for(n_vars = 0; n_vars < vars.length(); n_vars++) {
      if(Rcpp::IntegerVector::is_na(vars[n_vars])) break;
    }
    if(n_vars == 0) break;

    if(max_vars - vars.ncol() > 0) {
      Rcpp::IntegerMatrix fill_mat(vars.nrow(), max_vars - vars.ncol());
      std::fill(fill_mat.begin(), fill_mat.end(), Rcpp::IntegerVector::get_na());
      vars = Rcpp::cbind(vars, fill_mat);
    }
    disj = rbind1(disj, vars);

    Rcpp::NumericVector preds = Rcpp::as<Rcpp::NumericVector>(model["preds"]);
    currentRho = getRho(y, currentEstimates, preds, y_bin);

    currentEstimates = currentEstimates + learning_rate * currentRho * preds;

    Rprintf("\rIteration %d/%d (%.0f%%)", i, boosting_iter, i/boosting_iter * 100);
  }
  Rprintf("\rBoosting done\n");
  Rcpp::List ret; ret["disj"] = disj;
  return ret;
}


