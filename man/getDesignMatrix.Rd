% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bits.R
\name{getDesignMatrix}
\alias{getDesignMatrix}
\title{Design matrix for the set of terms}
\usage{
getDesignMatrix(X, disj, neg.offset)
}
\arguments{
\item{X}{The original predictor matrix or data frame of \eqn{p} input
variables}

\item{disj}{The term matrix which can, e.g., be extracted from a
fitted \code{BITS} model via $disj}

\item{neg.offset}{Numeric vector of length \eqn{p} that specifies the
negation offsets for all input variables, i.e., that specifies \eqn{m_j}
for the negated variable \eqn{X_j^c = m_j-X_j}}
}
\value{
The transformed design matrix
}
\description{
Transform the original predictor matrix X into the term design matrix
which contains for each term a corresponding column.
}
