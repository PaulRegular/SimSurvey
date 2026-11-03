#' Configure a generic fishery capture model
#'
#' Returns a function for the `catch_model` argument of [sim_fishery()]. Capture
#' proposals use a bounded encounter probability and an optional gear capacity.
#'
#' @param q Nonnegative capture coefficient, in grid-area units per gear-time
#' unit. It can exceed one; it is not a capture probability.
#' @param capacity_per_gear Maximum fish per gear unit in one operation, or `Inf`.
#'
#' @details
#' For gear amount G, duration T, selectivity s, and cell area A, the probability
#' is `1 - exp(-q * G * T * s / A)`. Binomial proposals are sampled from the
#' available age-length counts. If necessary, proposals are sampled without
#' replacement to respect `floor(capacity_per_gear * G)`.
#' This is an idealized encounter-and-capacity model, not a mechanistic model of
#' every gear's soak-time response. Alternative capture functions may be supplied.
#' The symbol q is conventional in fisheries catch/effort models; see
#' <https://www.fao.org/4/x5685e/x5685e04.htm>. It differs from the dimensionless
#' age/length catchability probabilities supplied as q to [sim_survey()].
#'
#' @return A function with arguments `available` and `operation`. `available` is a
#' table with `age`, `length`, `N`, `selectivity`, and `cell_area`; `operation`
#' includes `gear_amount`, `duration`, and `effort`, plus identifiers/location.
#' The function returns one nonnegative whole-fish proposal per available row.
#'
#' @examples
#' capture <- sim_capture(q = 0.1, capacity_per_gear = 1)
#' available <- data.frame(age = 1:2, length = c(10, 20), N = c(100, 50),
#'                         selectivity = c(0.5, 1), cell_area = 25)
#' capture(available, list(gear_amount = 10, duration = 4, effort = 40))
#' @export
sim_capture <- function(q, capacity_per_gear = Inf) {
  if (!is.numeric(q) || length(q) != 1 || !is.finite(q) || q < 0) {
    stop("q must be one finite nonnegative catchability.", call. = FALSE)
  }
  if (!is.numeric(capacity_per_gear) || length(capacity_per_gear) != 1 ||
      is.na(capacity_per_gear) || capacity_per_gear < 0) {
    stop("capacity_per_gear must be nonnegative or Inf.", call. = FALSE)
  }
  force(q)
  force(capacity_per_gear)
  capture <- function(available, operation) {
    prob <- numeric(nrow(available))
    selected <- available$selectivity > 0
    prob[selected] <- -expm1(-q * operation$effort * available$selectivity[selected] /
                              available$cell_area[selected])
    caught <- stats::rbinom(length(prob), available$N, prob)
    capacity <- floor(capacity_per_gear * operation$gear_amount)
    .thin_counts(caught, min(sum(caught), capacity))
  }
  attr(capture, "sim_capture") <- list(q = q, capacity_per_gear = capacity_per_gear)
  capture
}
