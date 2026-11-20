#' @useDynLib attitudedyn, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim pchisq qnorm pnorm na.omit na.exclude na.pass
#'   model.frame model.matrix model.response terms complete.cases setNames
#'   coef vcov logLik nobs AIC BIC confint predict simulate fitted formula
#' @importFrom broom tidy glance
NULL
