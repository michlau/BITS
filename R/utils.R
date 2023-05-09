#' @export
calcNMSE <- function(preds, y) {
  mse <- mean((preds - y)^2)
  v <- mean((mean(y) - y)^2)
  return(mse/v)
}
