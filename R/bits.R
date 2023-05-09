#' @export
#' @importFrom glmnet glmnet cv.glmnet relax.glmnet
#' @importFrom stats as.formula glm step binomial gaussian
BITS <- function(X, y, Z = NULL,
                 max.vars = 3, gamma = 0.01,
                 modify.vars = FALSE, remove.vars = FALSE,
                 boosting.iter = 50, learning.rate = 0.1,
                 lambda = NULL, alpha = 1, nfolds = 0, nlambda = 100,
                 relax = TRUE, gamma2 = 0.5,
                 adjust.shady.int = TRUE, term.select = "glinternet",
                 negate = TRUE) {
  y_bin <- !any(!(y %in% 0:1))
  p <- ncol(X)
  N <- nrow(X)
  X <- as.matrix(X)
  use.Z <- !is.null(Z)

  if(use.Z) Z <- as.matrix(Z)

  neg.offset <- rep(0, p)
  if(is.logical(negate) && negate) {
    ranges <- apply(X, 2, range)
    to.negate <- ranges[1,] == 0
    neg.offset[to.negate] <- ranges[2, to.negate]
  } else if(length(negate) == p) {
    neg.offset <- negate
  }

  models <- list()
  disj <- matrix(NA_integer_, nrow = 0, ncol = max.vars)
  rhoVector <- numeric()
  priorProb <- mean(y)

  if(y_bin)
    initialModel <- 0.5 * log(priorProb/(1-priorProb))
  else
    initialModel <- priorProb

  currentEstimates <- rep(initialModel, N)
  evaluatedWeakLearners <- array(0, c(boosting.iter, N))

  for(i in 1:boosting.iter) {
    if(y_bin) {
      probs <- 1/(1+exp(-2 * currentEstimates))
      gradient <- 2*probs - 2*y
    } else {
      gradient <- currentEstimates - y
    }

    model <- greedyFit(X, -gradient, Z = Z, neg_offset = neg.offset,
                       max_vars = max.vars, gamma = gamma,
                       force_model = (i == 1), adjust_shady_int = adjust.shady.int,
                       modify_vars = modify.vars, remove_vars = remove.vars)

    if(sum(!is.na(model$vars)) == 0) break

    models[[i]] <- model
    new.vars <- cbind(model$vars, matrix(NA_integer_, nrow=nrow(model$vars), ncol=max.vars-ncol(model$vars)))
    disj <- rbind(disj, new.vars)

    currentRho <- getRho(y, currentEstimates, model$preds, y_bin)

    rhoVector[i] <- currentRho
    currentEstimates <- currentEstimates + learning.rate * currentRho * model$preds

    cat("\r", sprintf("Iteration %d/%d (%.0f%%)", i, boosting.iter, i/boosting.iter * 100))
    flush.console()
  }

  cat("Boosting done\n")

  # disj <- boosting(X, y, Z = Z, y_bin = y_bin, max_vars = max.vars, gamma = gamma,
  #                  boosting_iter = boosting.iter, learning_rate = learning.rate,
  #                  adjust_shady_int = adjust.shady.int)$disj

  disj <- dontNegateSinglePredictors(disj)
  disj.unique <- unique(t(apply(disj, 1, sort.int, method = "quick", na.last = TRUE)))
  disj.unique <- disj.unique[rowMeans(is.na(disj.unique)) < 1,,drop=FALSE]
  disj.unique <- disj.unique[,colMeans(is.na(disj.unique)) < 1,drop=FALSE]

  if(nrow(disj.unique) == 0) {
    warning("An appropriate model could not be fitted.")
    return(NULL)
  }

  dm <- getDesignMatrix(X, Z, disj.unique, neg.offset)

  if(term.select == "step") {
    n.terms <- nrow(disj.unique)
    dat <- data.frame(y = y, dm)
    if(y_bin) fam <- binomial() else fam <- gaussian()
    biggest <- as.formula(paste("y ~", paste(colnames(dm), sep="", collapse = " + ")))
    if(use.Z) {
      p.Z <- ncol(Z)
      smallest <- as.formula(paste("y ~", paste(colnames(dm)[1:(n.terms+p.Z)], sep="", collapse = " + ")))
      lin.mod1 <- step(glm(smallest, data = dat, family = fam), direction = "both",
                       scope = list(lower= smallest, upper = biggest), trace = 0)
      detected.ints <- names(lin.mod1$coefficients)[-(1:(n.terms+p.Z+1))]
      Z.names <- colnames(dm)[(n.terms+1):(n.terms+p.Z)]
      smallest <- as.formula(paste("y ~", paste(c(Z.names, detected.ints),
                                                sep="", collapse = " + ")))
      biggest <- as.formula(paste("y ~", paste(c(colnames(dm)[1:(n.terms+p.Z)], detected.ints), sep="", collapse = " + ")))
      tmp <- p.Z + length(detected.ints)
    } else {
      detected.ints <- character() -> Z.names
      smallest <- y ~ 1
      tmp <- 0
    }
    lin.mod2 <- step(glm(smallest, data = dat, family = fam), direction = "both",
                     scope = list(lower= smallest, upper = biggest), trace = 0)
    detected.terms <- names(lin.mod2$coefficients)[-(1:(tmp+1))]
    final.formula <- as.formula(paste("y ~", paste(c(detected.terms, Z.names, detected.ints), sep="", collapse = " + ")))
    lin.mod <- glm(final.formula, data = dat, family = fam)
    s <- NULL
    main.disj <- disj.unique[match(detected.terms, colnames(dm)),,drop=FALSE]
    Z.disj <- disj.unique[match(detected.ints, colnames(dm)[(n.terms+length(Z.names)+1):ncol(dm)]),,drop=FALSE]
  } else if(!use.Z || term.select == "elnet") {
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
    main.disj <- NULL -> Z.disj
  } else {
    n.terms <- nrow(disj.unique); p.Z <- ncol(Z)
    tmp <- getDesignMatrixGL(X, Z, disj.unique, neg.offset)
    # cat("gglasso fitting...\n")
    dm.gl <- tmp$dm; groups <- tmp$groups; pf <- tmp$factors
    # pf <- rep(1, length(unique(groups)))
    y2 <- y
    if(y_bin) y2[y2 == 0] <- -1
    if(nfolds > 0) {
      lin.mod <- gglasso::cv.gglasso(dm.gl, y2, group = groups, pf = pf,
                                     loss = ifelse(y_bin, "logit", "ls"),
                                     lambda = lambda, nfolds = nfolds, nlambda = nlambda)
      s <- lin.mod$lambda.min
    } else {
      lin.mod <- gglasso::gglasso(dm.gl, y2, group = groups, pf = pf,
                                  loss = ifelse(y_bin, "logit", "ls"),
                                  lambda = lambda, nlambda = nlambda)
      s <- lambda
    }
    main.disj <- NULL -> Z.disj
  }

  ret <- list(lin.mod = lin.mod, disj = disj.unique, s = s, main.disj = main.disj, Z.disj = Z.disj,
              y_bin = y_bin, neg.offset = neg.offset,
              gamma2 = gamma2) ###
  class(ret) <- "BITS"
  return(ret)
}

#' @export
predict.intstump <- function(model, X, Z = NULL) {
  coef <- model$model$coef
  if(length(coef) == 1) return(rep(coef, nrow(X)))
  vars <- model$vars
  dm <- cbind(1, getDesignMatrix(X, Z, matrix(vars, nrow=1)))
  return(as.numeric(dm %*% matrix(coef, ncol=1)))
}

#' @export
#' @importFrom glmnet predict.glmnet
#' @importFrom gglasso gglasso
predict.BITS <- function(model, X, Z = NULL) {
  if(!is.null(Z)) Z <- as.matrix(Z)
  dm <- getDesignMatrix(X, Z, model$disj, model$neg.offset)
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
  } else {
    # glinternet
    dm.gl <- getDesignMatrixGL(X, Z, model$disj, model$neg.offset)$dm
    pred <- predict(model$lin.mod, dm.gl, type = "link", s = model$s)
    if(model$y_bin) pred <- 1/(1+exp(-pred))
    return(pred)
  }
}

dontNegateSinglePredictors <- function(disj) {
  n_vars <- rowSums(!is.na(disj[,,drop=FALSE]))
  disj[n_vars == 1, 1] <- abs(disj[n_vars == 1, 1])
  disj
}

#' @export
greedy.fit <- function(X, y, Z = NULL, neg.offset, max_vars = 3, gamma = 0.001, force_model = FALSE) {
  greedyFit(X, y, Z, neg.offset, max_vars, gamma, force_model)
}

#' @export
#' @importFrom stats lm.fit predict.glm gaussian binomial
greedy.fit.old <- function(X, y, Z = NULL, max.vars = 3, gamma = 0.01, force.model = FALSE) {
  p <- ncol(X); N <- nrow(X)
  use.Z <- !is.null(Z)
  vars <- integer()
  y_bin <- !any(!(y %in% 0:1))
  if(y_bin) fam <- gaussian() else fam <- binomial(link = "logit")

  best.vars <- vars
  best.mod <- lm.fit(matrix(1, nrow = N), y)
  class(best.mod) <- "glm"
  best.preds <- predict(best.mod, type = "response")
  best.score <- best.mod$deviance
  best.score <- calcScore(best.preds, y, y_bin)

  for(i in 1:max.vars) {
    found.better.model <- FALSE
    vars <- best.vars
    for(j in 1:(2*p)) {
      k <- ifelse(j <= p, j, -j)
      vars[i] <- k
      X.tmp <- X[,vars,drop=FALSE]
      X.tmp[,vars < 0,drop=FALSE] <- 1 - X.tmp[,vars < 0,drop=FALSE]
      interaction.feature <- apply(X.tmp, 1, prod)
      if(!use.Z)
        dm <- cbind(1, interaction.feature)
      else
        dm <- cbind(1, interaction.feature, Z, Z * interaction.feature)
      mod <- lm.fit(dm, y)
      class(mod) <- "glm"
      preds <- predict(mod, type = "response")
      score <- mod$deviance + gamma * i
      score <- calcScore(preds, y, y_bin) + gamma * i
      if(score <= best.score || (i == 1 && !found.better.model && force.model)) {
        best.score <- score; best.vars <- vars; best.mod <- mod; best.preds <- preds
        found.better.model <- TRUE
      }
    }
    if(!found.better.model) break
  }

  return(list(score = best.score, model = best.mod, preds = best.preds, vars = best.vars))
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
getDesignMatrix <- function(X, Z = NULL, disj, neg.offset, Z.interactions = TRUE) {
  use.Z <- !is.null(Z)
  p.new <- nrow(disj); N <- nrow(X)
  if(use.Z) p.Z <- ncol(Z) else p.Z <- 0
  if(!use.Z) dm <- matrix(nrow = N, ncol = p.new)
  else if(use.Z & Z.interactions) dm <- matrix(nrow = N, ncol = (p.Z+1) * p.new + p.Z)
  else dm <- matrix(nrow = N, ncol = p.new + p.Z)
  colnames(dm) <- paste("V", 1:ncol(dm), sep="")
  if(use.Z) {
    dm[,(p.new+1):(p.new+p.Z)] <- Z
    colnames(dm)[(p.new+1):(p.new+p.Z)] <- paste("Z", 1:p.Z, sep="")
  }
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
    col.start <- p.new + p.Z + (i-1)*p.Z + 1; col.end <- col.start + p.Z - 1
    if(use.Z & Z.interactions) {
      dm[,col.start:col.end] <- Z * interaction.feature
      colnames(dm)[col.start:col.end] <- paste(paste("T", i, sep=""), "Z", 1:p.Z, sep="")
    }
  }
  return(dm)
}

# Custom normalized Frobenius norm such that
# the norm is longer dependent on the sample size
frobenius.norm <- function(X) sqrt(sum(X^2)/nrow(X))

#' @export
getDesignMatrixGL <- function(X, Z, disj, neg.offset) {
  use.Z <- !is.null(Z)
  if(!use.Z) stop("Z must be supplied!")
  p.new <- nrow(disj); N <- nrow(X)
  p.Z <- ncol(Z)
  dm <- matrix(nrow = N, ncol = p.new + p.Z + p.new*p.Z*3)
  colnames(dm) <- paste("V", 1:ncol(dm), sep="")
  dm[,(p.new+1):(p.new+p.Z)] <- Z
  colnames(dm)[(p.new+1):(p.new+p.Z)] <- paste("Z", 1:p.Z, sep="")
  groups <- 1:(p.new+p.Z)
  factors <- rep(1, p.new+p.Z)
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

    col.start <- p.new + p.Z + (i-1)*p.Z*3 + 1
    for(j in 1:p.Z) {
      dm[,col.start] <- interaction.feature;
      dm[,col.start+1] <- Z[,j]; dm[,col.start+2] <- Z[,j] * interaction.feature
      colnames(dm)[col.start:(col.start+2)] <- paste(c("T", "Z", "I"), paste(i,j,sep=":"), sep="")
      groups <- c(groups, rep(groups[length(groups)]+1, 3))
      frob <- frobenius.norm(dm[,col.start:(col.start+2)])
      if(frob < 1e-30) frob <- 1
      dm[,col.start:(col.start+2)] <- dm[,col.start:(col.start+2)] / frob ###
      factors <- c(factors, frob)
    }
  }
  factors[1:(p.new+p.Z)] <- apply(dm[,1:(p.new+p.Z)], 2, function(x) frobenius.norm(matrix(x, ncol=1)))
  factors[factors < 1e-30] <- 1
  for(i in 1:(p.new+p.Z)) dm[,i] <- dm[,i] / factors[i] ###
  factors <- rep(1, length(factors)) ###
  return(list(dm=dm, groups=groups, factors=factors))
}

#' @export
get.included.vars <- function(model) {
  lin.mod <- model$lin.mod
  if(inherits(lin.mod, "glm")) {
    main.disj <- model$main.disj; Z.disj <- model$Z.disj
  } else if(inherits(lin.mod, "glmnet") || inherits(lin.mod, "cv.glmnet")) {
    if(is.null(model$s)) stop("The regularization parameter $s has to be properly set!")
    if(inherits(lin.mod, "cv.glmnet")) lin.mod <- lin.mod$glmnet.fit
    # lambda.ind <- match(model$s, lin.mod$lambda)[1]
    zero.tol <- 1e-10
    lambda.ind <- which(abs(model$s - lin.mod$lambda) < zero.tol)[1]
    beta <- as.numeric(lin.mod$beta[,lambda.ind])
    p.new <- nrow(model$disj)
    p.Z <- (length(beta) - p.new) / (p.new + 1)
    # main.disj <- model$disj[beta[1:p.new] != 0,,drop=FALSE]
    main.disj <- model$disj[abs(beta[1:p.new]) > zero.tol,,drop=FALSE]
    Z.disj <- NULL
    if(length(beta) > p.new + p.Z) {
      for(i in 1:p.Z) {
        col.start <- p.new + p.Z + (i-1)*p.new + 1; col.end <- col.start + p.new - 1
        # Z.disj <- rbind(Z.disj, model$disj[beta[col.start:col.end] != 0,,drop=FALSE])
        Z.disj <- rbind(Z.disj, model$disj[abs(beta[col.start:col.end]) > zero.tol,,drop=FALSE])
      }
      Z.disj <- unique(Z.disj)
    }
  } else {
    # glinternet
    if(is.null(model$s)) stop("The regularization parameter $s has to be properly set!")
    if("cv.gglasso" %in% class(lin.mod)) lin.mod <- lin.mod$gglasso.fit
    zero.tol <- 1e-10
    lambda.ind <- which(abs(model$s - lin.mod$lambda) < zero.tol)[1]
    beta <- lin.mod$beta[,lambda.ind]
    p.new <- nrow(model$disj)
    p.Z <- (length(beta) - p.new) / (p.new*3 + 1)
    chosen.beta <- abs(beta) > zero.tol
    main.disj <- matrix(nrow = 0, ncol = ncol(model$disj))
    Z.disj <- NULL
    for(i in 1:p.new) {
      ind <- which(startsWith(names(beta), paste("T", i, ":", sep="")))
      if(any(chosen.beta[c(i, ind)])) {
        main.disj <- rbind(main.disj, model$disj[i,,drop=FALSE])
      }
      ind <- which(startsWith(names(beta), paste("I", i, ":", sep="")))
      if(any(chosen.beta[ind])) {
        Z.disj <- rbind(Z.disj, model$disj[i,,drop=FALSE])
      }
    }
  }

  return(list(main.disj = main.disj, Z.disj = Z.disj))
}

#' @export
get.ideal.penalty <- function(model, X, y, Z = NULL, choose = "min") {
  lambda <- model$lin.mod$lambda
  preds <- predict(model, X, Z = Z)
  y_bin <- setequal(unique(y), c(0, 1))

  val.res <- data.frame()
  for(i in 1:length(lambda)) {
    scores <- calcScorePerObservation(as.numeric(preds[,i]), y, y_bin)
    m <- mean(scores)
    se <- sd(scores)/sqrt(length(scores))
    val.res <- rbind(val.res, data.frame(s = lambda[i], score = m, se = se, score.plus.1se = m + se))
  }

  min.ind <- min(which(val.res$score == min(val.res$score)))
  if(choose == "1se") {
    max.val <- val.res$score.plus.1se[min.ind]
    min.ind <- min(which(val.res$score <= max.val))
  }

  best.s <- lambda[min.ind]
  return(list(val.res = val.res, best.s = best.s))
}

#' @export
get.ideal.penalty2 <- function(model, X, y, Z = NULL, choose = "min") {
  if(!inherits(model$lin.mod, "relaxed")) return(get.ideal.penalty(model, X, y, Z, choose))
  lambda <- model$lin.mod$lambda
  gamma2 <- c(1, 0.75, 0.5, 0.25, 0)
  # preds <- list()
  y_bin <- setequal(unique(y), c(0, 1))
  val.res <- data.frame()

  for(g in 1:length(gamma2)) {
    model$gamma2 <- gamma2[g]
    preds <- predict(model, X, Z = Z)
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

#' @export
linearize <- function(model, X, y, Z = NULL) {
  disj <- get.included.vars(model)
  linearize.restricted(X, y, Z, disj$main.disj, disj$Z.disj)
}

#' @export
linearize.restricted <- function(X, y, Z = NULL, main.disj, Z.disj = NULL) {
  if(!is.null(Z)) Z <- as.matrix(Z)
  y_bin <- !any(!(y %in% 0:1))
  if(y_bin) fam <- binomial() else fam <- gaussian()
  dm <- getDesignMatrix2(X, Z, main.disj, Z.disj)

  df <- data.frame(y = y, dm)
  # colnames(df)[-1] <- paste("D", 1:ncol(dm), sep="")
  lin.mod <- glm(y ~ ., data = df, family = fam)
  # class(lin.mod) <- "glm"

  ret <- list(lin.mod = lin.mod, main.disj = main.disj, Z.disj = Z.disj)
  class(ret) <- "linBITS"
  return(ret)
}

#' @export
predict.linBITS <- function(model, X, Z = NULL) {
  if(!is.null(Z)) Z <- as.matrix(Z)
  dm <- getDesignMatrix2(X, Z, model$main.disj, model$Z.disj)
  df <- data.frame(dm)
  # if(ncol(df) > 0) colnames(df) <- paste("D", 1:ncol(dm), sep="")
  return(predict(model$lin.mod, df, type = "response"))
}

#' @export
getDesignMatrix2 <- function(X, Z = NULL, disj, Z.disj) {
  use.Z <- !is.null(Z)
  p.new <- nrow(disj); N <- nrow(X)
  if(use.Z) {
    p.Z <- ncol(Z)
    if(!is.null(Z.disj)) p.new.Z <- nrow(Z.disj) else p.new.Z <- 0
  } else {
    p.Z <- 0
    p.new.Z <- 0
  }
  if(!use.Z) dm <- matrix(nrow = N, ncol = p.new)
  else dm <- matrix(nrow = N, ncol = p.new + p.Z + p.new.Z*p.Z)
  if(ncol(dm) > 0) colnames(dm) <- 1:ncol(dm)
  if(use.Z) {
    dm[,(p.new+1):(p.new+p.Z)] <- Z
    colnames(dm)[(p.new+1):(p.new+p.Z)] <- paste("Z", 1:p.Z, sep="")
  }
  if(p.new > 0) {
    for(i in 1:nrow(disj)) {
      vars <- disj[i,]
      vars <- vars[!is.na(vars)]
      abs.vars <- abs(vars)
      X.tmp <- X[,abs.vars,drop=FALSE]
      X.tmp[,vars < 0] <- 1 - X.tmp[,vars < 0]
      interaction.feature <- apply(X.tmp, 1, prod)
      dm[,i] <- interaction.feature
      colnames(dm)[i] <- paste(vars, collapse=":")
    }
  }
  if(use.Z && p.new.Z > 0) {
    for(i in 1:p.new.Z) {
      vars <- Z.disj[i,]
      vars <- vars[!is.na(vars)]
      abs.vars <- abs(vars)
      X.tmp <- X[,abs.vars,drop=FALSE]
      X.tmp[,vars < 0] <- 1 - X.tmp[,vars < 0]
      interaction.feature <- apply(X.tmp, 1, prod)
      col.start <- p.new + p.Z + (i-1)*p.Z + 1; col.end <- col.start + p.Z - 1
      dm[,col.start:col.end] <- Z * interaction.feature
      colnames(dm)[col.start:col.end] <- paste(paste(vars, collapse=":"), ":Z", 1:p.Z, sep="")
    }
  }
  # dm <- cbind(1, dm)
  return(dm)
}

#' Gene-environment (GxE) interaction test based on boosted linear models
#'
#' This function takes a fitted \code{BITS} model and independent test
#' data as input for testing if there is a general GxE interaction.
#' This hypothesis test is based on a likelihood-ratio test.
#'
#' In detail, the null hypothesis
#' \deqn{H_0: \delta_1 = \ldots = \delta_B = 0}
#' using the supplied linear model
#' \deqn{g(E[Y]) = \beta_0 + \sum_{i=1}^B \beta_i \times 1[C_i] + \delta_0
#' \times E + \sum_{i=1}^B \delta_i \times 1[C_i] \cdot E}
#' is tested.
#'
#' @param model A fitted \code{BITS} model (i.e., a model created via
#'   \code{\link{BITS}})
#' @param X Matrix or data frame of input data.
#'   This object should correspond to the matrix for fitting the model.
#' @param y Response vector. 0-1 coding for binary outcomes.
#' @param Z Quantitative covariable supplied as a matrix or data frame
#' @return A list containing
#'   \item{\code{Deviance}}{The deviance used for performing the
#'     likelihood-ratio test}
#'   \item{\code{p.value}}{The p-value of the test}
#'
#' @importFrom stats anova
#' @export
gxe.test.boosting <- function(model, X, y, Z) {
  if(class(model) != "BITS") stop("The supplied model has to be of class 'BITS'!")
  if(is.null(Z)) return(NULL)

  ### GOOD
  mod.complete <- linearize(model, X, y, Z)$lin.mod
  included.vars <- get.included.vars(model)
  mod.reduced <- linearize.restricted(X, y, Z, included.vars$main.disj, matrix(nrow=0, ncol=ncol(included.vars$main.disj)))$lin.mod

  ###
  # included.vars <- get.included.vars(model)
  # mod.complete <- linearize.restricted(X, y, Z, model$disj, included.vars$Z.disj)$lin.mod
  # mod.reduced <- linearize.restricted(X, y, Z, model$disj, matrix(nrow=0, ncol=ncol(included.vars$main.disj)))$lin.mod

  ### GOOD
  res <- anova(mod.reduced, mod.complete, test="LRT")
  ret <- list(Deviance = res$Deviance[2], p.value = res$`Pr(>Chi)`[2])
  if(is.na(ret$p.value)) ret$p.value <- 1

  ###
  # Variance component test
  # included.vars <- get.included.vars(model)
  # dm <- getDesignMatrix2(X, NULL, included.vars$main.disj, NULL)
  # y_bin <- !any(!(y %in% 0:1))
  # library(penalized)
  # hm <- iSKAT::GESAT(dm, matrix(y, ncol=1), as.matrix(Z), out_type = ifelse(y_bin, "D", "C"),
  #                    type="liu")
  # p.val <- hm$pvalue[1]
  # if(is.null(p.val)) p.val <- 1
  # ret <- list(p.value = p.val, Deviance = 0)

  return(ret)
}

# perm.test <- function(diff, alternative, n.perm = 10000, mu = 0, stat = "t") {
#   if(stat == "t")
#     erg.orig <- as.numeric(t.test(diff - mu, alternative = alternative)$statistic)
#   else
#     erg.orig <- mean(diff - mu)
#   n <- length(diff)
#   null.stats <- vector(mode="numeric", length=n.perm)
#   for(i in 1:n.perm) {
#     perm <- sample(c(-1,1), n, replace = TRUE)
#     perm.diff <- perm * (diff - mu)
#     if(stat == "t")
#       null.stats[i] <- as.numeric(t.test(perm.diff, alternative = alternative)$statistic)
#     else
#       null.stats[i] <- mean(perm.diff)
#   }
#   null.stats <- sort(null.stats)
#   if(alternative == "greater") {
#     p.value <- mean(null.stats >= erg.orig)
#   } else if(alternative == "less") {
#     p.value <- mean(null.stats <= erg.orig)
#   } else {
#     p.value <- mean(abs(null.stats) >= abs(erg.orig))
#   }
#   return(p.value)
# }

#' Term importance test based on boosted linear models
#'
#' This function takes a fitted \code{BITS} model and independent test
#' data as input for testing if the included terms are influential with respect
#' to the outcome.
#' This hypothesis test is based on a likelihood-ratio test.
#'
#' In detail, the null hypotheses
#' \deqn{H_0: \beta_j = \delta_j = 0}
#' using the linear model
#' \deqn{g(E[Y]) = \beta_0 + \sum_{i=1}^B \beta_i \times 1[C_i] + \delta_0
#' \times E + \sum_{i=1}^B \delta_i \times 1[C_i] \cdot E}
#' are tested for each \eqn{j \in \lbrace 1,\ldots,B \rbrace}
#' if a quantitative covariable \code{Z} was included.
#' Otherwise, the null hypotheses
#' \deqn{H_0: \beta_j = 0}
#' using the linear model
#' \deqn{g(E[Y]) = \beta_0 + \sum_{i=1}^B \beta_i \times 1[C_i] +
#' \delta_0 \times E}
#' are tested.
#'
#' @param model A fitted \code{BITS} model (i.e., a model created via
#'   \code{\link{BITS}})
#' @param X Matrix or data frame of input data.
#'   This object should correspond to the matrix for fitting the model.
#' @param y Response vector. 0-1 coding for binary outcomes.
#' @param Z Optional quantitative covariables supplied as a matrix or
#'   data frame. Only used (and required) if the model was fitted using them.
#' @return A data frame consisting of three columns,
#'   \item{\code{var}}{The tested term,}
#'   \item{\code{vim}}{The associated variable importance, and}
#'   \item{\code{p.value}}{The corresponding p-value for testing if the term
#'     is influential.}
#'
#' @importFrom stats anova
#' @export
importance.test.boosting <- function(model, X, y, Z = NULL) {
  if(class(model) != "BITS") stop("The supplied model has to be of class 'BITS'!")

  tmp <- get.included.vars(model)
  main.disj <- tmp$main.disj; Z.disj <- tmp$Z.disj
  disj <- unique(rbind(main.disj, Z.disj))
  n_terms <- nrow(disj)

  # mod.complete <- linearize(model, X, y, Z = Z)$lin.mod
  # main.disj <- addSubTerms(main.disj); Z.disj <- addSubTerms(Z.disj) ###
  # mod.complete <- linearize.restricted(X, y, Z, main.disj, Z.disj)$lin.mod ###
  main.disj.orig <- main.disj; Z.disj.orig <- Z.disj ###

  vims <- data.frame(matrix(nrow = n_terms, ncol = 3))
  colnames(vims) <- c("var", "vim", "p.value")
  if(n_terms == 0) return(vims)
  vims$var <- getPredictorNames(translateDisj(disj, X))
  vims$vim <- 0 -> vims$p.value

  for(i in 1:n_terms) {
    main.disj <- addSubTerm(main.disj.orig, term = disj[i,]); Z.disj <- addSubTerm(Z.disj.orig, term = disj[i,]) ###
    mod.complete <- linearize.restricted(X, y, Z, main.disj, Z.disj)$lin.mod ###

    main.rem <- which(apply(main.disj, 1, function(x) all.equal(x, disj[i,])) == "TRUE")
    if(length(main.rem) > 0)
      main.disj.rem <- main.disj[-main.rem,,drop=FALSE]
    else
      main.disj.rem <- main.disj

    if(!is.null(Z.disj)) {
      Z.rem <- which(apply(Z.disj, 1, function(x) all.equal(x, disj[i,])) == "TRUE")
      if(length(Z.rem) > 0)
        Z.disj.rem <- Z.disj[-Z.rem,,drop=FALSE]
      else
        Z.disj.rem <- Z.disj
    } else {
      Z.disj.rem <- NULL
    }

    mod.reduced <- linearize.restricted(X, y, Z, main.disj.rem, Z.disj.rem)$lin.mod

    res <- anova(mod.reduced, mod.complete, test="LRT")
    vims$vim[i] <- res$Deviance[2]
    vims$p.value[i] <- res$`Pr(>Chi)`[2]
  }

  vims <- vims[order(vims$vim, decreasing = TRUE),]; rownames(vims) <- 1:nrow(vims)
  return(vims)
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
gammaPath <- function(X, y, Z = NULL,
                      gmin = function(gmax) 0.0001 * gmax,
                      steps = 20) {
  y_bin <- !any(!(y %in% 0:1))
  p <- ncol(X)
  N <- nrow(X)
  X <- as.matrix(X)
  use.Z <- !is.null(Z)

  if(use.Z) Z <- as.matrix(Z)

  priorProb <- mean(y)
  if(y_bin) {
    initialModel <- 0.5 * log(priorProb/(1-priorProb))
    probs <- 1/(1+exp(-2 * initialModel))
    gradient <- 2*probs - 2*y
  } else {
    initialModel <- priorProb
    gradient <- initialModel - y
  }
  mse.0 <- mean((gradient - mean(gradient))^2)

  model <- greedyFit(X, -gradient, Z = Z, neg_offset = rep(0, p),
                     max_vars = 1, gamma = 0, force_model = FALSE)
  mse.1 <- model$score

  # Get gmax such that
  # mse.0 = mse.1 + gmax
  gmax <- mse.0 - mse.1
  if(gmax <= 0) gmax <- 0.1

  gmin <- gmin(gmax)
  exp(seq(log(gmin), log(gmax), length.out = steps))
}

#' @export
cv.BITS <- function(X, y, Z = NULL,
                    nfolds = 5, choose = "min",
                    gmin = function(gmax) 0.0001 * gmax, steps = 20,
                    ...) {
  N <- nrow(X)
  y_bin <- setequal(unique(y), c(0, 1))
  use_Z <- !is.null(Z)
  folds <- cut(seq(1, N), breaks=nfolds, labels=FALSE)
  shuffle <- sample(N)
  X_shuffle <- X[shuffle,]
  y_shuffle <- y[shuffle]
  if(use_Z) Z_shuffle <- Z[shuffle,,drop=FALSE]

  gamma <- gammaPath(X, y, Z, gmin = gmin, steps = steps)
  cv.res <- data.frame()
  scores <- matrix(nrow = N, ncol = length(gamma))
  lambdas <- matrix(nrow = nfolds, ncol = length(gamma))

  for(i in 1:nfolds) {
    test_ind <- which(folds == i, arr.ind = TRUE)
    train_ind <- -test_ind
    X_train <- X_shuffle[train_ind,,drop=FALSE]
    y_train <- y_shuffle[train_ind]
    X_val <- X_shuffle[test_ind,,drop=FALSE]
    y_val <- y_shuffle[test_ind]
    Z_train <- NULL -> Z_val
    if(use_Z) {
      Z_train <- Z_shuffle[train_ind,,drop=FALSE]
      Z_val <- Z_shuffle[test_ind,,drop=FALSE]
    }

    for(j in 1:length(gamma)) {
      model <- BITS(X_train, y_train, Z = Z_train, nfolds = 0,
                    gamma = gamma[j], ...)
      model$s <- get.ideal.penalty(model, X_val, y_val, Z = Z_val, choose = "min")$best.s
      preds <- predict(model, X_val, Z = Z_val)
      scores[test_ind, j] <- calcScorePerObservation(preds, y_val, y_bin)
      lambdas[i, j] <- model$s
    }
  }

  for(j in 1:length(gamma)) {
    m <- mean(scores[,j])
    se <- sd(scores[,j])/sqrt(length(scores[,j]))
    s <- exp(mean(log(lambdas[,j]))) # Geometric mean
    cv.res <- rbind(cv.res, data.frame(gamma = gamma[j], s = s, score = m, se = se, score.plus.1se = m + se))
  }

  min.ind <- which.min(cv.res$score)
  if(choose == "1se") {
    max.val <- cv.res$score.plus.1se[min.ind]
    min.ind <- min(which(cv.res$score <= max.val))
  }

  best.gamma <- gamma[min.ind]; best.s <- cv.res$s[min.ind]

  model <- BITS(X, y, Z = Z, nfolds = 0,
                gamma = best.gamma, lambda = best.s, ...)
  cv <- list(best.gamma = best.gamma, best.s = best.s, cv.res = cv.res)
  model$cv <- cv
  return(model)
}

#' @useDynLib BITS
#' @importFrom Rcpp sourceCpp
NULL


