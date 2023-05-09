#' @export
BITS.bagging <- function(X, y, Z = NULL,
                         train.frac = 0.632,
                         bagging.iter = 50,
                         cv = FALSE, ...) {
  ret <- list()
  ret$models <- list()
  ret$bags <- list()
  ret$oob <- list()
  ret$X <- X
  ret$y <- y
  ret$Z <- Z
  N <- length(y)
  ret$preds <- rep(0.0, N) -> oob_n_models
  for(i in 1:bagging.iter) {
    bag <- sample(1:N, floor(train.frac * N), replace=FALSE)
    oob <- setdiff(1:N, bag)
    ret$bags[[i]] <- bag; ret$oob[[i]] <- oob
    if(is.null(Z)) {
      Z_temp <- NULL
      Z_oob <- NULL
    } else {
      Z_temp <- Z[bag,,drop=FALSE]
      Z_oob <- Z[oob,,drop=FALSE]
    }
    if(cv) {
      ret$models[[i]] <- cv.BITS(X[bag,,drop=FALSE], y[bag], Z = Z_temp, ...)
    } else {
      ret$models[[i]] <- BITS(X[bag,,drop=FALSE], y[bag], Z = Z_temp, ...)
    }
    ret$preds[oob] <- ret$preds[oob] + predict(ret$models[[i]], X[oob,,drop=FALSE], Z = Z_oob)
    oob_n_models[oob] <- oob_n_models[oob] + 1

    cat("\r", sprintf("Iteration %d/%d (%.0f%%)", i, bagging.iter, i/bagging.iter * 100))
    flush.console()
  }
  cat("\n")
  ret$preds <- ret$preds/oob_n_models
  class(ret) <- "BITS.bagged"
  return(ret)
}

get.bagging.terms <- function(model, min.term.fraction = 0.5) {
  disj.all <- lapply(model$models, function(x) get.included.vars(x))
  main.disj.all <- lapply(disj.all, function(x) x$main.disj)
  Z.disj.all <- lapply(disj.all, function(x) x$Z.disj)

  max.vars <- lapply(disj.all, function(x) lapply(x, ncol))
  max.vars <- max(unlist(max.vars))
  main.disj.all <- lapply(main.disj.all, function(x) cbind(x, matrix(nrow=nrow(x), ncol=max.vars-ncol(x))))
  Z.disj.all <- lapply(Z.disj.all, function(x) cbind(x, matrix(nrow=nrow(x), ncol=max.vars-ncol(x))))

  do.call(rbind, main.disj.all)

}

#' @export
gxe.test.bagging <- function(model, X, y, Z) {
  if(class(model) != "BITS.bagged") stop("The supplied model has to be of class 'BITS.bagged'!")
  if(is.null(Z)) return(NULL)

  bagging.iter <- length(model$models)
  p.vals <- numeric() -> devs
  for(i in 1:bagging.iter) {
    test <- gxe.test.boosting(model$models[[i]], X[model$oob[[i]],,drop=FALSE], y[model$oob[[i]]], Z[model$oob[[i]],,drop=FALSE])
    p.vals[i] <- test$p.value; devs[i] <- test$Deviance
  }

  gammin <- 0.05
  qs <- numeric()
  for(q in seq(gammin, 1, 0.01)) {
    qs <- c(qs, min(1, quantile(p.vals/q, q, na.rm = TRUE)))
  }
  p.value <- min(1, (1-log(gammin)) * min(qs))

  ### median instead of optimized quantile
  ### p.value <- min(1, median(p.vals)*2)

  Deviance <- mean(devs)

  return(list(Deviance = Deviance, p.value = p.value))
}

#' @export
importance.test.bagging <- function(model, X, y, Z = NULL) {
  if(class(model) != "BITS.bagged") stop("The supplied model has to be of class 'BITS.bagged'!")

  bagging.iter <- length(model$models)
  vims <- list()
  for(i in 1:bagging.iter) {
    Z.oob <- NULL
    if(!is.null(Z)) Z.oob <- Z[model$oob[[i]],,drop=FALSE]
    vims[[i]] <- importance.test.boosting(model$models[[i]], X[model$oob[[i]],,drop=FALSE], y[model$oob[[i]]], Z.oob)
  }

  vims <- lapply(vims, function(x) {
    x$p.value <- pmin(1, x$p.value * nrow(x))
    x
  })

  for(i in 1:bagging.iter) {
    if(i == 1) {
      a <- vims[[1]]
      a$var <- as.character(a$var); a <- a[order(a$var),]
    }
    else {
      b <- vims[[i]]
      b$var <- as.character(b$var); b <- b[order(b$var),]
      match.ind <- b$var %in% a$var
      a <- rbind(a, b[!match.ind,])
      a <- a[order(a$var),]
    }
  }

  vim.list <- list() -> p.list
  for(v in a$var) {
    vim.list[[v]] <- rep(0, bagging.iter)
    p.list[[v]] <- rep(1, bagging.iter)
    for(i in 1:bagging.iter) {
      b <- vims[[i]]
      if(v %in% b$var) {
        vim.list[[v]][i] <- b$vim[b$var == v]
        p.list[[v]][i] <- b$p.value[b$var == v]
      }
    }
    p.list[[v]][is.na(p.list[[v]])] <- 1
  }

  gammin <- 0.05
  ps <- list()
  for(v in a$var) {
    qs <- numeric()
    for(q in seq(gammin, 1, 0.01)) {
      qs <- c(qs, min(1, quantile(p.list[[v]]/q, q)))
    }
    ps[[v]] <- min(1, (1-log(gammin)) * min(qs))
  }

  ### median instead of optimized quantile
  # ps <- list()
  # for(v in a$var)
  #   ps[[v]] <- min(1, median(p.list[[v]])*2)

  vim.mean <- lapply(vim.list, mean)
  a$vim <- unname(unlist(vim.mean))
  a$p.value <- unname(unlist(ps))

  # ToDo: Test for next time!
  # FDR
  # Or: Hierarchy in terms or main effects in testing GLM
  # if(FDR && !is.null(a$vim)) {
  #   lp <- length(a$p.value)
  #   i <- lp:1L
  #   o <- order(a$p.value, decreasing = TRUE)
  #   ro <- order(o)
  #   q <- sum(1/(1L:lp))
  #   a$p.value <- pmin(1, cummin(q / i * a$p.value[o]))[ro]
  # }

  if(is.null(a$vim)) {
    a$vim <- numeric(); a$p.value <- numeric()
  } else {
    a <- a[order(a$vim, decreasing = TRUE),]; rownames(a) <- 1:nrow(a)
  }
  return(a)
}

