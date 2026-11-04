small_grid <- function() {
  make_grid(x_range = c(-10, 10), y_range = c(-10, 10), res = c(5, 5),
            shelf_width = 4, shelf_depth = 50, depth_range = c(0, 100),
            strat_breaks = c(0, 101), strat_splits = 1, method = "linear")
}

small_abundance <- function() {
  sim_abundance(
    ages = 1:3, years = 2001:2003,
    R = sim_R(log_mean = log(2000), log_sd = 0, random_walk = FALSE),
    Z = sim_Z(log_mean = log(0.4), log_sd = 0, phi_age = 0, phi_year = 0),
    growth = sim_vonB(Linf = 80, L0 = 10, K = 0.2, log_sd = 0.1,
                      length_group = 3)
  )
}

zero_covariance <- function(x, ages, years, cells) {
  array(0, dim = c(length(ages), length(years), length(cells)),
        dimnames = list(age = ages, year = years, cell = cells))
}

small_population <- function() {
  sim_distribution(small_abundance(), grid = small_grid(),
                   ays_covar = zero_covariance,
                   depth_par = function(data) rep(0, nrow(data)))
}

small_survey <- function(...) {
  sim_survey(small_population(), q = function(x) rep(1, length(x)),
             trawl_dim = c(1, 5), set_den = 0, min_sets = 2, ...)
}

local_worker_library <- function(envir = parent.frame()) {
  withr::local_envvar(c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)),
                     .local_envir = envir)
  withr::defer(foreach::registerDoSEQ(), envir = envir)
}
