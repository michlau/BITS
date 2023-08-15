#' Modes of inheritance
#'
#' Determining the most plausible modes of inheritance (MOIs).
#'
#' For each considered SNP (single nucleotide polymorphism) and each considered
#' MOI (additive, dominant, and recessive), likelihood-ratio tests are carried
#' out testing the association of the SNP in the considered MOI with the
#' outcome. For each SNP, the MOI that yields the highest association with the
#' outcome, i.e., the MOI that yields the lowest p-value, is returned.
#'
#' @param X Matrix or data frame of \eqn{p} SNPs coded as 0,1,2
#' @param y Numeric vector of a binary or continuous outcome
#' @param fast Shall fastglm instead of glm be employed for model fitting?
#'   Default is \code{FALSE}. However, for high-dimensional problems, setting
#'   \code{fast = TRUE} is recommended.
#' @return A character vector of length \eqn{p} containing the identified MOI,
#'   where \code{"A"} is the additive MOI, \code{"D"} is the dominant MOI, and
#'   \code{"R"} is the recessive MOI
#' @references
#' \itemize{
#'   \item Lau, M., Schikowski, T. & Schwender, H. (2023).
#'   Boosting interaction tree stumps. To be submitted.
#'   \item Scherer, N., Sekula, P., Pfaffelhuber, P. & Schlosser, P. (2021).
#'   pgainsim: an R-package to assess the mode of inheritance for quantitative
#'   trait loci in GWAS. Bioinformatics, 37(18):3061–3063.
#'   \doi{https://doi.org/10.1093/bioinformatics/btab150}
#'   \item Petersen, A. K., Krumsiek, J., Wägele, B., Theis, F. J., Wichmann,
#'   H.-E., Gieger, C. & Suhre, K. (2012). On the hypothesis-free testing of
#'   metabolite ratios in genome-wide and metabolome-wide association studies.
#'   BMC Bioinformatics, 13:120, 2012.
#'   \doi{https://doi.org/10.1186/1471-2105-13-120}
#' }
#'
#' @importFrom fastglm fastglm
#' @importFrom stats glm anova binomial gaussian pchisq logLik
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

#' Applying modes of inheritance
#'
#' Transform SNPs according to the supplied modes of inheritance (MOIs).
#'
#' @param X Matrix or data frame of SNPs coded as 0,1,2
#' @param modes Character vector of MOIs that can be obtained using
#'   \code{\link{MOI}}.
#' @return A matrix or data frame of SNPs using the given MOIs
#'
#' @export
applyMOI <- function(X, modes) {
  X[,modes == "D"] <- X[,modes == "D"] > 0
  X[,modes == "R"] <- X[,modes == "R"] == 2
  apply(X, 2, as.numeric)
}

