
#' Simulate length using von Bertalanffy notation
#'
#' @description This function returns a function to use inside \code{\link{sim_distribution}}.
#' Given parameters, it generates an approximate age-length key and inverse
#' age-length key by sampling from a von Bertalanffy relationship. Note that length is
#' discretized by rounding.
#'
#' @param Linf     Mean asymptotic length
#' @param L0       Length at birth
#' @param K        Growth rate parameter
#' @param log_sd   Standard deviation of the relationship in log scale
#' @param plot     Produce a simple plot of the simulated values?
#'
#' @export
#'

sim_growth <- function(Linf = 120, L0 = 5, K = 0.1, log_sd = 0.1, plot = FALSE) {
  function(ages = NULL) {

    ## Draw many samples from the vonB relationship, round and calculate proportions
    ## (This is a bruit force approach; there's probably a better way)
    age <- rep(ages, 100000)
    pred_length <- Linf - (Linf - L0) * exp(-K * age)
    log_length <- rnorm(length(age), log(pred_length), sd = log_sd)
    length <- round(exp(log_length))

    tab <- table(age, length)
    inv_alk <- prop.table(tab, margin = 1)
    alk <- prop.table(tab, margin = 2)

    if (plot) {
      image(x = as.numeric(rownames(alk)), y = as.numeric(colnames(alk)), z = inv_alk,
            xlab = "Age", ylab = "Length", main = "Inverse ALK",
            col = colorRampPalette(c("white", "black"))(100))
      box()
    }

    list(alk = alk, inv_alk = inv_alk)

  }
}

