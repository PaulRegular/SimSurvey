

#' Calculate intraclass correlation
#'
#' This is a simple function for calculating intraclass correlation. It uses
#' \code{\link[lme4]{lmer}} to run the formula described here:
#' https://en.wikipedia.org/wiki/Intraclass_correlation
#'
#' @param x      Response variable
#' @param group  Group
#'
#' @return Returns estimate of intraclass correlation.
#'
#' @export
#'

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
