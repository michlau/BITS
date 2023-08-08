#' @importFrom fastglm fastglm
#' @importFrom stats glm anova binomial gaussian pchisq
#' @export
MOI <- function(X, y, fast = FALSE) {
  y_bin <- !any(!(y %in% 0:1))
  p <- ncol(X); N <- nrow(X)
  modes <- rep("A", p)
  if(y_bin) fam <- binomial() else fam <- gaussian()
  form.full <- y ~ A; form.red <- y ~ 1
  if(fast) {
    mod.red <- fastglm(matrix(1,ncol=1,nrow=N), y, family = fam)
    for(i in 1:p) {
      A <- X[,i]
      mod.A <- fastglm(cbind(1, A), y, family = fam)
      llr.stat <- 2*as.numeric(logLik(mod.A) - logLik(mod.red))
      p.A <- 1-pchisq(llr.stat, 1)

      A <- as.numeric(X[,i] > 0)
      mod.D <- fastglm(cbind(1, A), y, family = fam)
      llr.stat <- 2*as.numeric(logLik(mod.D) - logLik(mod.red))
      p.D <- 1-pchisq(llr.stat, 1)

      A <- as.numeric(X[,i] == 2)
      mod.R <- fastglm(cbind(1, A), y, family = fam)
      llr.stat <- 2*as.numeric(logLik(mod.R) - logLik(mod.red))
      p.R <- 1-pchisq(llr.stat, 1)

      min.moi.ind <- which.min(c(p.A, p.D, p.R))
      if(length(min.moi.ind) == 0) min.moi.ind <- 1
      modes[i] <- c("A", "D", "R")[min.moi.ind]
    }
  } else {
    for(i in 1:p) {
      A <- X[,i]
      mod.A <- glm(form.full, family = fam); p.A <- 1
      mod.A.red <- glm(form.red, family = fam)
      p.A <- anova(mod.A.red, mod.A, test = "LRT")$`Pr(>Chi)`[2]

      A <- as.numeric(X[,i] > 0)
      mod.D <- glm(form.full, family = fam); p.D <- 1
      mod.D.red <- glm(form.red, family = fam)
      p.D <- anova(mod.D.red, mod.D, test = "LRT")$`Pr(>Chi)`[2]

      A <- as.numeric(X[,i] == 2)
      mod.R <- glm(form.full, family = fam); p.R <- 1
      mod.R.red <- glm(form.red, family = fam)
      p.R <- anova(mod.R.red, mod.R, test = "LRT")$`Pr(>Chi)`[2]

      min.moi.ind <- which.min(c(p.A, p.D, p.R))
      if(length(min.moi.ind) == 0) min.moi.ind <- 1
      modes[i] <- c("A", "D", "R")[min.moi.ind]
    }
  }
  modes
}

#' @export
applyMOI <- function(X, modes) {
  X[,modes == "D"] <- X[,modes == "D"] > 0
  X[,modes == "R"] <- X[,modes == "R"] == 2
  apply(X, 2, as.numeric)
}

