#' Configure weight given length
#'
#' A deterministic relationship shared by the population and all fisheries.
#' @param alpha Nonnegative coefficient in `weight = alpha * length^beta`.
#' @param beta Positive exponent.
#' @param plot Logical. Plot when the returned function is used?
#' @details Defaults use centimetres and return kilograms per fish. They are the
#' reference relationship in Table 10 of Varkey et al. (2024), Canadian Data Report
#' of Fisheries and Aquatic Sciences 1381:
#' <https://publications.gc.ca/collections/collection_2024/mpo-dfo/Fs97-13-1381-eng.pdf>.
#' They are an illustrative relationship used for 3Ps cod, not a calibration of
#' every aspect of the population. No individual weight variation is imposed.
#' Custom functions supplied to [sim_abundance()] return expected nonnegative
#' kilograms for each input length, or one constant value.
#' @return A function of `length` returning kilograms per fish.
#' @examples
#' weight <- sim_weight()
#' weight(c(20, 40, 60))
#' @export
sim_weight <- function(alpha = 10^-5.2106, beta = 3.0879, plot = FALSE) {
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha < 0 ||
      !is.numeric(beta) || length(beta) != 1 || !is.finite(beta) || beta <= 0) {
    stop("alpha must be finite and nonnegative; beta must be finite and positive.")
  }
  force(alpha)
  force(beta)
  function(length = NULL) {
    if (!is.numeric(length) || any(!is.finite(length)) || any(length < 0)) {
      stop("Lengths must be finite nonnegative numbers.")
    }
    weight <- alpha * length^beta
    if (any(!is.finite(weight))) stop("The weight relationship produced nonfinite kilograms.")
    if (plot) graphics::plot(length, weight, type = "l", xlab = "Length (cm)", ylab = "Weight (kg)")
    weight
  }
}

#' Configure maturity probability
#'
#' Uses the same logistic calculation as [sim_logistic()], with defaults for a
#' historical maturity-at-length relationship rather than survey catchability.
#' @param k Nonnegative steepness, in inverse units of the input trait.
#' @param x50 Trait value at 50 percent maturity: length or age according to
#' `maturity_by` in [sim_abundance()].
#' @param plot Logical. Plot when the returned function is used?
#' @details Defaults approximate the female 3Ps cod relationship for 1996 in
#' Table 3 of Brattey and Morgan (1996), DFO Atlantic Fisheries Research Document
#' 96/92: <https://waves-vagues.dfo-mpo.gc.ca/Library/198884.pdf>.
#' Length is in centimetres and `k` is approximately 0.28 per centimetre. Maturity
#' varies by sex and period; this historical example is not a contemporary stock
#' calibration. SimSurvey applies one curve to all fish, without a female fraction.
#' Maturity is a probability, not a knife-edge classification at `x50`. Supply
#' age-appropriate `k` and `x50` when using age-based maturity.
#' @return A function of `x` returning probabilities in `[0, 1]`. Interchangeable
#' with [sim_logistic()] when parameters and trait basis agree.
#' @examples
#' maturity <- sim_maturity()
#' maturity(c(30, 47.81, 60))
#' @export
sim_maturity <- function(k = 0.28, x50 = 47.81, plot = FALSE) {
  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) || k < 0) {
    stop("Maturity steepness k must be finite and nonnegative.")
  }
  sim_logistic(k = k, x50 = x50, plot = plot)
}

