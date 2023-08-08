#' @export
#' @importFrom glmnet glmnet cv.glmnet relax.glmnet
#' @importFrom stats as.formula glm step binomial gaussian
BITS <- function(X, y,
                 max.vars = 3, gamma = NULL,
                 boosting.iter = 50, learning.rate = 0.1,
                 lambda = NULL, alpha = 1, nfolds = 0, nlambda = 100,
                 relax = TRUE, gamma2 = 0.5,
                 adjust.shady.int = TRUE, term.select = "glinternet",
                 negate = TRUE,
                 max.iter = function(p) 4*max.vars*p,
                 set_vars = NULL, reuse.terms = TRUE) {
  y_bin <- !any(!(y %in% 0:1))
  p <- ncol(X)
  N <- nrow(X)
  X <- as.matrix(X)

  neg.offset <- rep(0, p)
  if(is.logical(negate) && negate) {
    ranges <- apply(X, 2, range)
    to.negate <- ranges[1,] == 0
    neg.offset[to.negate] <- ranges[2, to.negate]
  } else if(length(negate) == p) {
    neg.offset <- negate
  }

  disj <- matrix(NA_integer_, nrow = 0, ncol = max.vars)
  rhoVector <- numeric()
  priorProb <- mean(y)

  if(y_bin)
    initialModel <- 0.5 * log(priorProb/(1-priorProb))
  else
    initialModel <- priorProb

  currentEstimates <- rep(initialModel, N)
  evaluatedWeakLearners <- array(0, c(boosting.iter, N))

  # Old parameters
  search.algo <- "complete"; optimize <- "correlation"; independent <- FALSE

  X2 <- apply(X, 2, function(x) {
    m <- max(abs(x))
    if(m > 1e-10) x/m else x
  })
  max.iter <- max.iter(2*p)
  if(max.iter <= 0 && is.null(set_vars)) set_vars <- initializeTerms(X2, neg_offset = rep(1, p), max_vars = max.vars)
  evaluated.terms <- rep(0, boosting.iter) -> possible.terms

  for(i in 1:boosting.iter) {
    if(y_bin) {
      probs <- 1/(1+exp(-2 * currentEstimates))
      gradient <- 2*probs - 2*y
    } else {
      gradient <- currentEstimates - y
    }

    grad2 <- gradient - mean(gradient)
    ###### grad.sd <- sd(gradient) * sqrt((N-1)/N)
    ###### grad2 <- (gradient - mean(gradient))/grad.sd
    model <- completeSearch(X2, -grad2,
                            neg_offset = rep(1, p),
                            max_vars = max.vars, gamma = gamma,
                            set_vars_R = set_vars, reuse_terms = reuse.terms, max_iter = max.iter,
                            force_model = (i == 1), adjust_shady_int = adjust.shady.int)
    if(reuse.terms) set_vars <- model$set_vars
    model$preds <- model$preds - mean(gradient)
    ###### model$preds <- model$preds * grad.sd - mean(gradient)

    evaluated.terms[i] <- model$evaluated_terms
    possible.terms[i] <- model$possible_terms

    if(sum(!is.na(model$vars)) == 0) break

    new.vars <- cbind(model$vars, matrix(NA_integer_, nrow=nrow(model$vars), ncol=max.vars-ncol(model$vars)))
    disj <- rbind(disj, new.vars)

    currentRho <- getRho(y, currentEstimates, model$preds, y_bin)

    rhoVector[i] <- currentRho
    currentEstimates <- currentEstimates + learning.rate * currentRho * model$preds

    cat("\r", sprintf("Iteration %d/%d (%.0f%%)", i, boosting.iter, i/boosting.iter * 100))
    flush.console()
  }

  cat("Boosting done\n")

  disj <- dontNegateSinglePredictors(disj)
  disj.unique <- unique(t(apply(disj, 1, sort.int, method = "quick", na.last = TRUE)))
  disj.unique <- disj.unique[rowMeans(is.na(disj.unique)) < 1,,drop=FALSE]
  disj.unique <- disj.unique[,colMeans(is.na(disj.unique)) < 1,drop=FALSE]

  if(nrow(disj.unique) == 0) {
    warning("An appropriate model could not be fitted.")
    return(NULL)
  }

  dm <- getDesignMatrix(X, disj.unique, neg.offset)

  if(term.select == "step") {
    n.terms <- nrow(disj.unique)
    dat <- data.frame(y = y, dm)
    if(y_bin) fam <- binomial() else fam <- gaussian()
    biggest <- as.formula(paste("y ~", paste(colnames(dm), sep="", collapse = " + ")))

    smallest <- y ~ 1
    tmp <- 0

    lin.mod2 <- step(glm(smallest, data = dat, family = fam), direction = "both",
                     scope = list(lower= smallest, upper = biggest), trace = 0)
    detected.terms <- names(lin.mod2$coefficients)[-(1:(tmp+1))]
    final.formula <- as.formula(paste("y ~", paste(detected.terms, sep="", collapse = " + ")))
    lin.mod <- glm(final.formula, data = dat, family = fam)
    s <- NULL
  } else {
    if(ncol(dm) == 1) dm <- cbind(dm, 1)
    if(y_bin) fam <- "binomial" else fam <- "gaussian"
    if(nfolds > 0) {
      lin.mod <- cv.glmnet(dm, y, family = fam, lambda = lambda, alpha = alpha, nfolds = nfolds, nlambda = nlambda,
                           relax = relax)
      s <- lin.mod$lambda.min
      if(relax) {
        s <- lin.mod$relaxed$lambda.min
        gamma2 <- lin.mod$relaxed$gamma.min
      }
    } else {
      lin.mod <- glmnet(dm, y, family = fam, lambda = lambda, alpha = alpha, nlambda = nlambda)
      if(relax) lin.mod <- relax.glmnet(lin.mod, x = dm, y = y, family = fam, lambda = lambda,
                                        alpha = alpha, nlambda = nlambda)
      s <- lambda
    }
  }

  ret <- list(lin.mod = lin.mod, disj = disj.unique, s = s,
              y_bin = y_bin, neg.offset = neg.offset,
              gamma2 = gamma2) ###
  ret$evaluated.terms <- evaluated.terms; ret$possible.terms <- possible.terms
  if(reuse.terms) {
    ret$set_vars <- set_vars
  }
  class(ret) <- "BITS"
  return(ret)
}

#' @export
#' @importFrom glmnet predict.glmnet
predict.BITS <- function(model, X) {
  dm <- getDesignMatrix(X, model$disj, model$neg.offset)
  if(inherits(model$lin.mod, "glm")) {
    return(predict(model$lin.mod, as.data.frame(dm), type = "response"))
  } else if(inherits(model$lin.mod, "glmnet") || inherits(model$lin.mod, "cv.glmnet")) {
    if(ncol(dm) == 1) dm <- cbind(dm, 1)
    gamma2 <- model$gamma2
    if(inherits(model$lin.mod, "relaxed")) {
      pred.penalized <- predict.glmnet(model$lin.mod, dm, s = model$s, type = "link")
      pred.unpenalized <- predict.glmnet(model$lin.mod$relaxed, dm, s = model$s, type = "link")
      pred <- (1-gamma2) * pred.unpenalized + gamma2 * pred.penalized
      if(model$y_bin) pred <- 1/(1+exp(-pred))
      return(pred)
    } else {
      return(predict(model$lin.mod, dm, s = model$s, type = "response", gamma = gamma2))
    }
  }
}

dontNegateSinglePredictors <- function(disj) {
  n_vars <- rowSums(!is.na(disj[,,drop=FALSE]))
  disj[n_vars == 1, 1] <- abs(disj[n_vars == 1, 1])
  disj
}

#' @export
intstump <- function(X, y, negate = TRUE, max.vars = 3, gamma = NULL, max.iter = function(p) 4*max.vars*p, adjust.shady.int = TRUE) {
  p <- ncol(X)
  max.iter <- max.iter(2*p)
  neg.offset <- rep(0, p)
  if(is.logical(negate) && negate) {
    ranges <- apply(X, 2, range)
    to.negate <- ranges[1,] == 0
    neg.offset[to.negate] <- ranges[2, to.negate]
  } else if(length(negate) == p) {
    neg.offset <- negate
  }
  set_vars <- NULL
  X.factors <- apply(X, 2, function(x) max(abs(x)))
  X.factors[X.factors <= 1e-10] <- 1
  X2 <- X/matrix(X.factors, nrow=nrow(X), ncol=ncol(X), byrow=TRUE)
  if(max.iter <= 0) set_vars <- initializeTerms(X2, neg_offset = rep(1, p), max_vars = max.vars)

  yy <- mean(y)
  y2 <- y - yy
  stump <- completeSearch(X2, y2, rep(1, p), max.vars, gamma, set_vars, FALSE, max.iter,
                          TRUE, adjust.shady.int)
  stump$preds <- stump$preds + yy
  stump$model$preds <- stump$model$preds + yy

  # Adjust intercept (because y was centered)
  stump$model$coef[1] <- stump$model$coef[1] + yy
  # Adjust coefficients (because X was scaled to [0, 1])
  vars <- stump$vars
  for(i in 1:nrow(vars)) {
    v <- abs(vars[i,]); v <- v[!is.na(v)]
    stump$model$coef[1+i] <- stump$model$coef[1+i]/prod(X.factors[v])
  }

  class(stump) <- "intstump"
  stump$neg.offset <- neg.offset
  stump
}

#' @export
predict.intstump <- function(model, X) {
  coef <- model$model$coef
  if(length(coef) == 1) return(rep(coef, nrow(X)))
  vars <- model$vars
  dm <- cbind(1, getDesignMatrix(X, vars, model$neg.offset))
  return(as.numeric(dm %*% matrix(coef, ncol=1)))
}

calcScore <- function(preds, y, y_bin) {
  if(y_bin) return(-2*sum(log(y * preds + (1-y) * (1-preds))))
  else return(mean((preds - y)^2))
}

#' @importFrom stats optim
getRho <- function(y, oldEstimates, evaluatedWeakLearners, y_bin) {
  if(y_bin) {
    calcScore <- function(x) {
      return(sum(log(1+exp(2*(oldEstimates + x * evaluatedWeakLearners))) - 2 * y * (oldEstimates + x*evaluatedWeakLearners)))
    }
    return(optim(1, calcScore, method = "Brent", lower = 0, upper = 8)$par[1])
  } else {
    calcScore <- function(x) {
      return(sum((y - oldEstimates - x * evaluatedWeakLearners)^2)/2)
    }
    return(sum((y - oldEstimates) * evaluatedWeakLearners)/sum((evaluatedWeakLearners)^2))
  }
}

#' @export
getDesignMatrix <- function(X, disj, neg.offset) {
  p.new <- nrow(disj); N <- nrow(X)
  dm <- matrix(nrow = N, ncol = p.new)
  colnames(dm) <- paste("V", 1:ncol(dm), sep="")
  for(i in 1:nrow(disj)) {
    vars <- disj[i,]
    vars <- vars[!is.na(vars)]
    abs.vars <- abs(vars)
    X.tmp <- X[,abs.vars,drop=FALSE]
    tmp.offset <- neg.offset[-vars[vars < 0]]
    tmp.offset <- matrix(tmp.offset, nrow=nrow(X.tmp), ncol=length(tmp.offset), byrow=TRUE)
    X.tmp[,vars < 0] <- tmp.offset - X.tmp[,vars < 0]
    interaction.feature <- apply(X.tmp, 1, prod)
    dm[,i] <- interaction.feature
    colnames(dm)[i] <- paste("T", i, sep="")
  }
  return(dm)
}

#' @export
get.included.vars <- function(model) {
  lin.mod <- model$lin.mod
  if(inherits(lin.mod, "glm")) {
    disj <- model$disj
  } else if(inherits(lin.mod, "glmnet") || inherits(lin.mod, "cv.glmnet")) {
    if(is.null(model$s)) stop("The regularization parameter $s has to be properly set!")
    if(inherits(lin.mod, "cv.glmnet")) lin.mod <- lin.mod$glmnet.fit
    # lambda.ind <- match(model$s, lin.mod$lambda)[1]
    zero.tol <- 1e-10
    lambda.ind <- which(abs(model$s - lin.mod$lambda) < zero.tol)[1]
    beta <- as.numeric(lin.mod$beta[,lambda.ind])
    p.new <- nrow(model$disj)
    # disj <- model$disj[beta[1:p.new] != 0,,drop=FALSE]
    disj <- model$disj[abs(beta[1:p.new]) > zero.tol,,drop=FALSE]
  }
  return(disj)
}

#' @importFrom stats sd
#' @importFrom logicDT calcAUC
#' @export
get.ideal.penalty <- function(model, X, y, choose = "min", metric = "auc") {
  lambda <- model$lin.mod$lambda
  preds <- predict(model, X)
  y_bin <- setequal(unique(y), c(0, 1))
  per.observation <- !y_bin || metric == "dev"
  if(y_bin) y <- as.integer(y)

  val.res <- data.frame()
  for(i in 1:length(lambda)) {
    if(per.observation) {
      scores <- calcScorePerObservation(as.numeric(preds[,i]), y, y_bin)
      m <- mean(scores)
      se <- sd(scores)/sqrt(length(scores))
      val.res <- rbind(val.res, data.frame(s = lambda[i], score = m, se = se, score.plus.1se = m + se))
    } else {
      m <- 1 - calcAUC(as.numeric(preds[,i]), y)
      val.res <- rbind(val.res, data.frame(s = lambda[i], score = m))
    }
  }

  min.ind <- min(which(val.res$score == min(val.res$score)))
  if(choose == "1se" && per.observation) {
    max.val <- val.res$score.plus.1se[min.ind]
    min.ind <- min(which(val.res$score <= max.val))
  }

  best.s <- lambda[min.ind]
  return(list(val.res = val.res, best.s = best.s))
}

#' @importFrom stats sd
get.ideal.penalty2 <- function(model, X, y, choose = "min") {
  if(!inherits(model$lin.mod, "relaxed")) return(get.ideal.penalty(model, X, y, Z, choose))
  lambda <- model$lin.mod$lambda
  gamma2 <- c(1, 0.75, 0.5, 0.25, 0)
  # preds <- list()
  y_bin <- setequal(unique(y), c(0, 1))
  val.res <- data.frame()

  for(g in 1:length(gamma2)) {
    model$gamma2 <- gamma2[g]
    preds <- predict(model, X)
    for(i in 1:length(lambda)) {
      scores <- calcScorePerObservation(as.numeric(preds[,i]), y, y_bin)
      m <- mean(scores)
      se <- sd(scores)/sqrt(length(scores))
      val.res <- rbind(val.res, data.frame(gamma2 = gamma2[g], s = lambda[i], score = m, se = se, score.plus.1se = m + se))
    }
  }

  min.ind <- min(which(val.res$score == min(val.res$score)))
  if(choose == "1se") {
    max.val <- val.res$score.plus.1se[min.ind]
    min.ind <- min(which(val.res$score <= max.val))
  }

  best.s <- val.res$s[min.ind]; best.gamma2 <- val.res$gamma2[min.ind]
  return(list(val.res = val.res, best.s = best.s, best.gamma2 = best.gamma2))
}

calcScorePerObservation <- function(preds, y, y_bin) {
  scores <- vector(mode = "numeric", length = length(preds))
  for(i in 1:length(preds)) scores[i] <- calcScore(preds[i], y[i], y_bin)
  scores
}

#' @export
custom.lm <- function(X, y) {
  customLm(X, y)
}

#' @export
matrix.vector.mult <- function(Z, y) {
  matrixVectorMult(Z, y)
}

addSubTerms <- function(disj) {
  if(is.null(disj) || nrow(disj) == 0) return(disj)
  disj2 <- disj
  for(i in 1:nrow(disj)) disj2 <- addSubTerm(disj2, disj[i,])
  disj2
}

#' @importFrom utils combn
addSubTerm <- function(disj, term) {
  if(is.null(disj) || nrow(disj) == 0 || is.null(term)) return(disj)
  disj2 <- disj
  if(sum(!is.na(term)) > 1) {
    vars <- term; vars <- vars[!is.na(vars)]
    # disj2 <- rbind(disj2, cbind(abs(vars), matrix(NA_integer_, nrow=length(vars), ncol=ncol(disj)-1)))
    # if(length(vars) > 2) {
    #   ind <- combn(1:length(vars), 2, simplify = TRUE)
    #   tmp <- t(apply(ind, 2, function(x) vars[x]))
    #   disj2 <- rbind(disj2, cbind(tmp, matrix(NA_integer_, nrow=length(vars), ncol=ncol(disj)-2)))
    # }
    n <- length(vars)
    for(i in 1:(n-1)) {
      ind <- combn(vars, i, simplify = TRUE)
      tmp <- t(ind)
      disj2 <- rbind(disj2, cbind(tmp, matrix(NA_integer_, nrow=nrow(tmp), ncol=ncol(disj)-i)))
    }
  }
  disj2 <- dontNegateSinglePredictors(disj2)
  unique(disj2)
}

translateDisj <- function(disj, X) {
  if(is.null(colnames(X))) return(disj)
  translated <- t(apply(disj, 1, function(row) ifelse(is.na(row), NA, paste(ifelse(row < 0, "-", ""), colnames(X)[abs(row)], sep=""))))
  if(ncol(disj) == 1)
    translated <- t(translated)
  return(translated)
}

getPredictorNames <- function(real_disj, sort_conj = FALSE) {
  n_conj <- sum(rowSums(!is.na(real_disj)) > 0)
  disj2 <- split(real_disj[1:n_conj,,drop=FALSE], 1:n_conj)
  if(sort_conj)
    disj2 <- lapply(disj2, sort)
  disj2 <- lapply(disj2, function(x) x[!is.na(x)])
  disj2 <- lapply(disj2, paste, collapse=":")
  disj2 <- unlist(disj2, use.names = FALSE)
  return(disj2)
}

#' @export
gammaPath <- function(X, y,
                      gmin = function(gmax) 0.001 * gmax,
                      steps = 50) {
  y_bin <- !any(!(y %in% 0:1))
  p <- ncol(X)
  N <- nrow(X)
  X <- as.matrix(X)

  priorProb <- mean(y)
  if(y_bin) {
    initialModel <- 0.5 * log(priorProb/(1-priorProb))
    probs <- 1/(1+exp(-2 * initialModel))
    gradient <- 2*probs - 2*y
  } else {
    initialModel <- priorProb
    gradient <- initialModel - y
  }

  # Old parameters
  search.algo <- "complete"; optimize <- "correlation"

  X2 <- apply(X, 2, function(x) {
    m <- max(abs(x))
    if(m > 1e-10) x/m else x
  })
  grad2 <- gradient - mean(gradient)
  model <- completeSearch(X2, -grad2,
                          neg_offset = rep(0, p),
                          max_vars = 1, gamma = 0,
                          set_vars_R = NULL, reuse_terms = FALSE, max_iter = -1,
                          force_model = FALSE, adjust_shady_int = FALSE)
  # Get gmax such that 0 = corr.1 - gmax, since null model has zero correlation
  gmax <- model$corr

  if(gmax <= 0) gmax <- 0.1

  gmin <- gmin(gmax)
  exp(seq(log(gmax), log(gmin), length.out = steps))
}

compareNA <- function(v1, v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

combine.custom <- function(a, b) {
  depth <- function(this,thisdepth=0){
    if(!is.list(this)){
      return(thisdepth)
    }else{
      return(max(unlist(lapply(this,depth,thisdepth=thisdepth+1))))
    }
  }

  if(depth(a) > depth(b)) {
    res <- c(a, list(b))
  } else {
    res <- c(list(a), list(b))
  }
  return(res)
}

#' @importFrom foreach `%dopar%` `%do%`
ParMethod <- function(x) if(x) {`%dopar%`} else {`%do%`}

#' @importFrom foreach getDoParWorkers foreach
#' @export
BITS.complete <- function(X, y,
                          max.vars = 3, gamma = NULL,
                          boosting.iter = 50, learning.rate = 0.1,
                          lambda = NULL, alpha = 1, nfolds = 0, nlambda = 100,
                          relax = TRUE, gamma2 = 0.5,
                          adjust.shady.int = TRUE, term.select = "glinternet",
                          negate = TRUE,
                          max.iter = function(p) -1,
                          reuse.terms = TRUE, parallel = TRUE,
                          gmin = function(gmax) 0.001 * gmax, gsteps = 50) {
  X <- as.matrix(X)
  X2 <- apply(X, 2, function(x) {
    m <- max(abs(x))
    if(m > 1e-10) x/m else x
  })
  if(is.null(gamma)) gamma <- gammaPath(X, y, gmin = gmin, steps = gsteps)
  gamma <- sort(gamma, decreasing = TRUE)

  max.iter.FUN <- max.iter
  max.iter <- max.iter.FUN(2*ncol(X))
  if(max.iter <= 0)
    set_vars <- initializeTerms(X2, neg_offset = rep(1, ncol(X)), max_vars = max.vars)
  else
    set_vars <- NULL

  null.model <- BITS(X, y,
                     max.vars = max.vars, gamma = 0,
                     boosting.iter = boosting.iter, learning.rate = learning.rate,
                     lambda = lambda, alpha = alpha, nfolds = nfolds, nlambda = nlambda,
                     relax = relax, gamma2 = gamma2,
                     adjust.shady.int = adjust.shady.int, term.select = term.select,
                     negate = negate, max.iter = max.iter.FUN,
                     set_vars = set_vars, reuse.terms = reuse.terms)
  if(reuse.terms) set_vars <- null.model$set_vars
  null.model$gamma <- 0

  if(reuse.terms || max.iter <= 0 || getDoParWorkers() == 1) parallel <- FALSE
  `%op%` <- ParMethod(parallel)
  should.break <- length(gamma) + 1
  models <- foreach(i=1:length(gamma), .combine=combine.custom, .export = c()) %op%
  {
    if(i < should.break) {
      model <- BITS(X, y,
                    max.vars = max.vars, gamma = gamma[i],
                    boosting.iter = boosting.iter, learning.rate = learning.rate,
                    lambda = lambda, alpha = alpha, nfolds = nfolds, nlambda = nlambda,
                    relax = relax, gamma2 = gamma2,
                    adjust.shady.int = adjust.shady.int, term.select = term.select,
                    negate = negate, max.iter = max.iter.FUN,
                    set_vars = set_vars, reuse.terms = reuse.terms)
      if(reuse.terms) set_vars <- model$set_vars
      model$gamma <- gamma[i]

      if(i > 1 && all(dim(model$disj) == dim(null.model$disj)) && all(compareNA(model$disj, null.model$disj))) should.break <- i
      model
    }
  }
  models <- models[lengths(models) != 0]
  models <- models[order(unlist(lapply(models, function(x) x$gamma)), decreasing = TRUE)]
  models[[length(models) + 1]] <- null.model

  class(models) <- "BITS.list"
  models
}

#' @importFrom stats sd
#' @importFrom logicDT calcAUC
#' @export
get.ideal.model <- function(models, X, y, choose = "min", metric = "auc") {
  if(!inherits(models, "BITS.list")) stop("The first argument needs to be an object of class 'BITS.list' (generated using BITS.complete)!")
  y_bin <- setequal(unique(y), c(0, 1))
  per.observation <- !y_bin || metric == "dev"
  if(y_bin) y <- as.integer(y)

  val.res <- data.frame()

  for(i in 1:length(models)) {
    model <- models[[i]]
    lambda <- model$lin.mod$lambda
    preds <- predict(model, X)

    for(j in 1:length(lambda)) {
      if(per.observation) {
        scores <- calcScorePerObservation(as.numeric(preds[,j]), y, y_bin)
        m <- mean(scores)
        se <- sd(scores)/sqrt(length(scores))
        val.res <- rbind(val.res, data.frame(g = model$gamma, s = lambda[j], score = m, se = se, score.plus.1se = m + se))
      } else {
        m <- 1 - calcAUC(as.numeric(preds[,j]), y)
        val.res <- rbind(val.res, data.frame(g = model$gamma, s = lambda[j], score = m))
      }
    }
  }

  min.ind <- which.min(val.res$score)
  if(choose == "1se" && per.observation) {
    max.val <- val.res$score.plus.1se[min.ind]
    min.ind <- min(which(val.res$score <= max.val))
  }

  best.g <- val.res$g[min.ind]; best.s <- val.res$s[min.ind]
  return(list(val.res = val.res, best.g = best.g, best.s = best.s))
}

#' @importFrom stats sd
#' @importFrom logicDT calcAUC
#' @export
cv.BITS <- function(X, y, nfolds = 5, parallel = TRUE,
                    gmin = function(gmax) 0.001 * gmax, gsteps = 50,
                    choose = "min", metric = "auc", ...) {
  N <- nrow(X)
  y_bin <- setequal(unique(y), c(0, 1))
  if(y_bin) y <- as.integer(y)
  folds <- cut(seq(1, N), breaks=nfolds, labels=FALSE)
  shuffle <- sample(N)
  X_shuffle <- X[shuffle,]
  y_shuffle <- y[shuffle]
  per.observation <- !y_bin || metric == "dev"

  gamma <- gammaPath(X, y, gmin = gmin, steps = gsteps)

  cv.res <- data.frame()

  if(per.observation)
    scores <- matrix(Inf, nrow = N, ncol = length(gamma)+1)
  else
    scores <- matrix(Inf, nrow = nfolds, ncol = length(gamma)+1)
  lambdas <- matrix(Inf, nrow = nfolds, ncol = length(gamma)+1)

  for(i in 1:nfolds) {
    test_ind <- which(folds == i, arr.ind = TRUE)
    train_ind <- -test_ind
    X_train <- X_shuffle[train_ind,,drop=FALSE]
    y_train <- y_shuffle[train_ind]
    X_val <- X_shuffle[test_ind,,drop=FALSE]
    y_val <- y_shuffle[test_ind]

    models <- BITS.complete(X_train, y_train, gamma = gamma, parallel = parallel, ...)

    for(j in 1:length(models)) {
      model <- models[[j]]

      model$s <- get.ideal.penalty(model, X_val, y_val, choose = choose, metric = metric)$best.s
      preds <- predict(model, X_val)
      j2 <- which(c(0, gamma) == model$gamma)
      lambdas[i, j2] <- model$s
      if(per.observation)
        scores[test_ind, j2] <- calcScorePerObservation(preds, y_val, y_bin)
      else
        scores[i, j2] <- 1 - calcAUC(preds, y_val)
    }
  }

  gamma <- c(0, gamma)

  for(j in 1:length(gamma)) {
    m <- mean(scores[,j])
    s <- exp(mean(log(lambdas[,j]))) # Geometric mean
    if(per.observation) {
      se <- sd(scores[,j])/sqrt(length(scores[,j]))
      cv.res <- rbind(cv.res, data.frame(g = gamma[j], s = s, score = m, se = se, score.plus.1se = m + se))
    } else {
      cv.res <- rbind(cv.res, data.frame(g = gamma[j], s = s, score = m))
    }
  }

  min.ind <- which.min(cv.res$score)
  if(choose == "1se" && per.observation) {
    max.val <- cv.res$score.plus.1se[min.ind]
    min.ind <- min(which(cv.res$score <= max.val))
  }

  best.gamma <- gamma[min.ind]; best.s <- cv.res$s[min.ind]

  model <- BITS(X, y, nfolds = 0, gamma = best.gamma, lambda = best.s, ...)
  cv <- list(best.gamma = best.gamma, best.s = best.s, cv.res = cv.res)
  model$cv <- cv
  return(model)
}

#' @export plot.val.performance
plot.val.performance <- function(val.res) {
  library(ggplot2)
  p <- ggplot(val.res, aes(x=log(s), y=score, color = log(g), group = g)) +
    geom_line() +
    theme_bw() +
    xlab(expression(log(lambda))) + ylab("Score") +
    scale_colour_gradient(name = expression(log(gamma)), low="#00B4EF", high="#FF6C91")
  p
}

#' @export
calcNoPossibleTerms <- function(p, max.vars, negate = TRUE) {
  if(!negate) return(sum(choose(p, 1:max.vars)))
  tmp <- sum(choose(2*p, 1:max.vars))
  tmp - sum(p*choose(2*(p-1), 0:(max.vars-2)))
}

memory.test <- function() {
  memoryTest()
}

#' @useDynLib BITS
#' @importFrom Rcpp sourceCpp
NULL


