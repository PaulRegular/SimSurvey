
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



sim <- sim_distribution(pop = sim_abundance(years = 1:10),
                        grid = sim_grid(res = c(3.5, 3.5)),
                        space_covar = sim_sp_covar(range = 50, sd = 0.1),
                        ay_covar = sim_ay_covar(sd = 10,
                                                phi_age = 0.5,
                                                phi_year = 0.5),
                        depth_par = sim_parabola(alpha = 0, sigma = 50))

## Generate age-length keys
growth <- sim_growth()
g <- growth(ages = sim$ages)
lengths <- as.numeric(colnames(g$inv_alk))

## Expand abundance at age to abundance at age and length
## (todo: clean up this code...it's clumsy)
R <- t(replicate(ncol(g$inv_alk), sim$R)) * replicate(length(sim$years), g$inv_alk[as.character(sim$ages[1]), ])
dimnames(R) <- list(length = lengths, year = sim$years)
N0 <- g$inv_alk * replicate(ncol(g$inv_alk), sim$N0)
N <- sapply(as.character(sim$years),
            function(y) g$inv_alk * replicate(ncol(g$inv_alk), sim$N[, y]),
            simplify = "array")
dimnames(N) <- list(age = sim$ages, length = colnames(g$inv_alk), year = sim$years)
# all.equal(sim$N, apply(N, FUN = sum, MARGIN = c(1, 3)))

rowSums(head(sim$sp_N$N, 20) * g$inv_alk[head(as.character(sim$sp_N$age), 20), 1:10])

long_inv_alk <- as.data.frame.table(g$inv_alk, responseName = "prop")
sp_N <- merge(sim$sp_N, long_inv_alk, by = "age", all = TRUE)
## Not going to work well...computational drain.
## Simplify and minimize object size by simulating length of 'caught' fish


