

#' Calculate intraclass correlation
#'
#' A simple function for calculating intraclass correlation using [`lme4::lmer()`].
#' The formula follows the description provided on [Wikipedia](https://en.wikipedia.org/wiki/Intraclass_correlation).
#'
#' @param x Response variable.
#' @param group Grouping variable.
#'
#' @return An estimate of intraclass correlation.
#'
#' @export

icc <- function(x, group) {

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("lme4 is needed for icc to work. Please install it.", call. = FALSE)
  }
  mod <- lme4::lmer(x ~ 1 + (1 | group))
  var <- as.data.frame(lme4::VarCorr(mod))
  grpvar <- var$vcov[var$grp == "group"]
  resvar <- var$vcov[var$grp == "Residual"]
  grpvar / (grpvar + resvar)

}
