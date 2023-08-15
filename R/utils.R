#' Calculate the NMSE
#'
#' Computation of the normalized mean squared error.
#'
#' @param preds Numeric vector of predictions
#' @param y True outcomes
#' @return The NMSE
#'
#' @export
calcNMSE <- function(preds, y) {
  mse <- mean((preds - y)^2)
  v <- mean((mean(y) - y)^2)
  return(mse/v)
}
