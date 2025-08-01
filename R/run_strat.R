
#' Calculate common error statistics
#'
#' @param error A numeric vector of errors.
#'
#' @return A named vector of error statistics including:
#'
#' - **ME**: Mean error
#' - **MAE**: Mean absolute error
#' - **MSE**: Mean squared error
#' - **RMSE**: Root mean squared error
#'
#' @export

error_stats <- function(error) {
  c(ME = mean(error),
    MAE = mean(abs(error)),
    MSE = mean(error ^ 2),
    RMSE = sqrt(mean(error ^ 2)))
}


#' Prepare simulated data for stratified analysis
#'
#' Generate set details (`setdet`), length-frequency (`lf`), and age-frequency (`af`)
#' data for stratified analysis.
#'
#' @param sim Simulation object returned by [`sim_survey()`].
#' @param length_group Size of the length frequency bins.
#' @param alk_scale Spatial scale at which to construct and apply age-length keys:
#' `"division"`, `"strat"`, or `"set"`.
#'
#' @return A list containing:
#' - `setdet`: Set details
#' - `lf`: Length-frequency data
#' - `af`: Age-frequency data
#'
#' @export

strat_data <- function(sim, length_group = 3, alk_scale = "division") {

  n_measured <- n_aged <- n <- n_lengths <- ratio <- length_freq <- age_tot <- age_freq <-
    age_prop <- age <- NULL

  ## Extract setdet and samp objects, and add sim and year to samp data
  setdet <- sim$setdet
  samp <- sim$samp
  cols <- c("set", "sim", "year")
  if (alk_scale != "set") {
    cols <- c(cols, alk_scale)
  }
  samp <- merge(setdet[, cols, with = FALSE], samp, by = "set")

  ## Produce warnings if age sampling scale does not match scale requested for age-length key
  ag_check <- setdet[n_measured > 0, list(n_aged = sum(n_aged)), by = c("sim", "year", "strat")]
  if (alk_scale == "strat" && any(ag_check$n_aged == 0)) {
    warning("Age-length keys cannot be constructed at the strat level because ages were not sampled at every strat. This discrepancy may introduce bias. Consider running sim_survey with age_space_group = 'strat'")
  }
  ag_check <- setdet[n_measured > 0, list(n_aged = sum(n_aged)), by = c("sim", "year", "set")]
  if (alk_scale == "set" && any(ag_check$n_aged == 0)) {
    warning("Age-length keys cannot be constructed at the set level because ages were not sampled at every set. This discrepancy may introduce bias. Consider running sim_survey with age_space_group = 'set'")
  }

  ## Construct length-frequency table
  lf <- samp
  lf[, n := .N, by = "set"]
  lf <- lf[lf$measured, ]
  lf[, n_lengths := .N, by = "set"]
  lf$length <- group_lengths(lf$length, length_group)
  lf[, ratio := n_lengths / n] # calculate ratio measured
  lf <- lf[, list(length_freq = .N),
           by = c("set", "ratio", "length")]
  lf[, length_freq := length_freq / ratio] # scale frequencies
  lf[, ratio := NULL] # discard ratio column (no longer needed)
  setkeyv(lf, c("set", "length"))

  ## Add zeros to length-frequency table
  cj <- CJ(set = setdet$set, length = sort(unique(lf$length)), unique = TRUE)
  setkeyv(cj, c("set", "length"))
  cj <- merge(setdet[, cols, with = FALSE], cj, by = "set", all = TRUE)
  lf <- merge(cj, lf, by = c("set", "length"), all = TRUE)
  lf$length_freq[is.na(lf$length_freq)] <- 0 # replace NA's with 0's

  ## Construct age-length key by sim and year. Split by division or strat if requested.
  ag <- samp
  ag <- ag[ag$aged, ] # discard records without corresponding length and age values
  ag$length <- group_lengths(ag$length, length_group)
  alk <- ag[, list(age_freq = .N), by = c("sim", "year", "length", "age", alk_scale)] # combined alk
  alk[, age_tot := sum(age_freq), by = c("sim", "year", "length", alk_scale)]
  alk[, age_prop := age_freq / age_tot, by = c("sim", "year", "length", alk_scale)]
  setkeyv(alk, c("sim", "year", "length", alk_scale))

  ## Merge lf and alk objects
  alf <- merge(lf, alk, by = c("sim", "year", "length", alk_scale), all = TRUE, allow.cartesian = TRUE)
  alf <- alf[!is.na(length_freq)] # drop matches lacking length frequencies (further ensure bins not observed in the sets are not represented in the results)
  alf[, age_freq := ifelse(is.na(age), length_freq, length_freq * age_prop)] # pass 100% len freq if age is unknown for that length bin, otherwise apply proportion
  af <- alf[, list(age_freq = sum(age_freq, na.rm = TRUE)), by = c("set", "age")]
  setkeyv(af, c("set", "age"))

  ## Add zeros to age-frequency table (this step may be unnecessary)
  ages <- c(NA, seq(min(ag$age, na.rm = TRUE), max(ag$age, na.rm = TRUE)))
  af <- af[CJ(set = setdet$set, age = ages, unique = TRUE), ] # expand to missing groups in each set
  af$age_freq[is.na(af$age_freq)] <- 0 # replace NA's with 0's

  ## Simplify/order and return data
  lf <- lf[, c("set", "length", "length_freq"), with = FALSE]
  setkeyv(lf, c("set", "length"))
  setkeyv(af, c("set", "age"))
  setnames(lf, "length_freq", "n")
  setnames(af, "age_freq", "n")
  list(setdet = setdet, lf = lf, af = af)

}



#' Calculate stratified means, variances, and confidence intervals across groups
#'
#' This function is primarily designed for use within [`run_strat()`]. It first calculates
#' statistics at the stratum level and then computes broader summaries like total abundance.
#'
#' @param data A `data.table` with all grouping variables in stacked format. Must include
#' `strat_area` and `tow_area` for scaling values.
#' @param metric Name of the variable in `data.table` to summarize (e.g., `"number"`, `"mass"`).
#' @param strat_groups Grouping variables for fine-scale stratum-level means.
#' Must include `"strat"` and `"strat_area"`.
#' Example: `c("year", "species", "shiptrip", "NAFOdiv", "strat", "strat_area", "age")`
#' @param survey_groups Grouping variables for large-scale summary calculations
#' Example: `c("year", "species")`
#' @param confidence Confidence limit percentage (e.g., 95 for 95% CI).
#'
#' @return A `data.table` containing stratified estimates of abundance.
#'
#' @export

strat_means <- function(data = NULL, metric = NULL, strat_groups = NULL,
                        survey_groups = NULL, confidence = 95) {

  Nh <- strat_area <- tow_area <- Wh <- total <- sumYh <- nh <- gh <- meanYh <- varYh <-
    meanYst_lcl <- meanYst <- varYst <- df <- meanYst_ucl <- sumYst <- N <- sumYst_lcl <-
    sumYst_ucl <- NULL

  ## set-up
  lc <- (100 - confidence) / 200
  uc <- (100 - confidence) / 200 + (confidence / 100)
  d <- copy(data) # make a copy of the provided data.table
  d <- d[,  c(strat_groups,  metric),  with = FALSE]
  setnames(d, names(d), c(strat_groups, "metric"))
  setkeyv(d,  strat_groups)

  ## strat.tab includes strat-level means,  variances,  etc.
  strat_tab <- d[, list(sumYh = sum(metric), meanYh = mean(metric),
                        varYh = stats::var(metric), nh = .N), by = strat_groups]
  strat_tab[, Nh := strat_area / tow_area] # number of sample units in each strat
  strat_tab[, Wh := Nh / sum(Nh), by = survey_groups]
  strat_tab[, total := Nh * sumYh / nh]
  strat_tab[, gh := Nh * (Nh - nh) / nh]

  ## survey.tab includes large-scale means,  such as mean weight per tow,  and totals,  such as total biomass + associated confidence intervals
  survey_tab <- strat_tab[, list(n = sum(nh), N = sum(Nh), meanYst = sum(Wh * meanYh),
                                 varYst = (1 / ((sum(Nh)) ^ 2)) * sum(gh * varYh),
                                 df = ((sum(gh * varYh)) ^ 2) / (sum((gh ^ 2 * varYh ^ 2) / (nh - 1)))),
                          by = survey_groups]
  survey_tab[, meanYst_lcl := (meanYst - (sqrt(varYst)) * abs(stats::qt(lc,  df)))]
  survey_tab[, meanYst_ucl := (meanYst + (sqrt(varYst)) * abs(stats::qt(lc,  df)))]
  survey_tab[, sumYst := N * meanYst]
  survey_tab[, sumYst_lcl := (sumYst - abs(stats::qt(lc, df)) * N * sqrt(varYst))]
  survey_tab[, sumYst_ucl := (sumYst + abs(stats::qt(lc, df)) * N * sqrt(varYst))]
  survey_tab[sapply(survey_tab,  is.nan)] <- NA

  ## Rename cols
  survey_tab <- survey_tab[, c(survey_groups,  "n", "N", "df", "varYst",
                               "meanYst",  "meanYst_lcl", "meanYst_ucl",
                               "sumYst", "sumYst_lcl", "sumYst_ucl"),  with = FALSE]
  survey_tab$varYst <- sqrt(survey_tab$varYst) # convert to sd
  setnames(survey_tab, names(survey_tab), c(survey_groups, "sets", "sampling_units", "df", "sd",
                                            "mean", "mean_lcl", "mean_ucl",
                                            "total", "total_lcl", "total_ucl"))

  ## return results (here I only return the survey_tab to minimize details and object size)
  survey_tab

}


#' Run stratified analysis on simulated data
#'
#' @param sim Simulation object from [`sim_survey()`].
#' @param length_group Size of the length frequency bins used for both abundance-at-length
#' calculations and age-length-key construction. By default, this is inherited from the
#' value defined in [`sim_abundance()`] via the closure supplied to `sim_length` (`"inherit"`).
#' You may also supply a numeric value; however, mismatches in length groupings may cause
#' issues with [`strat_error()`] if true vs. estimated groupings are not aligned.
#' @param alk_scale Spatial scale at which to construct and apply age-length keys:
#' `"division"` or `"strat"`.
#' @param strat_data_fun Function used to prepare data for stratified analysis (e.g., [`strat_data()`]).
#' @param strat_means_fun Function used to calculate stratified means (e.g., [`strat_means()`]).
#'
#' @details The `strat_data_fun` and `strat_means_fun` arguments allow you to use custom
#' [`strat_data()`] and [`strat_means()`] functions.
#'
#' @return Adds stratified analysis results to the `sim` list:
#' - `"total_strat"`: Results for the total population
#' - `"length_strat"`: Results aggregated by length group
#' - `"age_strat"`: Results aggregated by age
#'
#' @examples
#' \donttest{
#' sim <- sim_abundance(ages = 1:5, years = 1:5,
#'                      R = sim_R(log_mean = log(1e+7)),
#'                      growth = sim_vonB(length_group = 1)) |>
#'   sim_distribution(grid = make_grid(res = c(20, 20)),
#'                    ays_covar = sim_ays_covar(sd = 1)) |>
#'   sim_survey(n_sims = 1, q = sim_logistic(k = 2, x0 = 3)) |>
#'   run_strat()
#' }
#'
#' @export

run_strat <- function(sim,
                      length_group = "inherit",
                      alk_scale = "division",
                      strat_data_fun = strat_data,
                      strat_means_fun = strat_means) {

  sim_length_group <- get("length_group", envir = environment(sim$sim_length))
  if (is.character(length_group) && length_group == "inherit") {
    length_group <- sim_length_group
  } else {
    if (length_group != sim_length_group) {
      warning(paste0("length_group value should be set to ", sim_length_group,
                     " to match the length group defined inside sim_abundance using sim_length",
                     "; a mismatch in length groupings will cause issues with strat_error",
                     " as true vs. estimated length groupings will be mismatched."))
    }
  }

  ## Prep data for strat_means function
  data <- strat_data_fun(sim, length_group = length_group, alk_scale = alk_scale)
  data$setdet <- data$setdet[, c("sim", "year", "division", "strat", "strat_area", "tow_area", "set", "n"), with = FALSE]
  data$lf <- merge(data$setdet[, setdiff(names(data$setdet), "n"), with = FALSE], data$lf, by = "set", all = TRUE)
  data$af <- merge(data$setdet[, setdiff(names(data$setdet), "n"), with = FALSE], data$af, by = "set", all = TRUE)

  ## Run strat calculations
  strat_args <- list(data = data$setdet,
                     metric = "n",
                     strat_groups = c("sim", "year", "division", "strat", "strat_area", "tow_area"),
                     survey_groups = c("sim", "year"))
  sim$total_strat <- do.call(strat_means_fun, strat_args)

  strat_args$data <- data$lf
  strat_args$strat_groups <- c(strat_args$strat_groups, "length")
  strat_args$survey_groups <- c(strat_args$survey_groups, "length")
  sim$length_strat <- do.call(strat_means_fun, strat_args)

  strat_args$data <- data$af
  strat_args$strat_groups[strat_args$strat_groups == "length"] <- "age"
  strat_args$survey_groups[strat_args$survey_groups == "length"] <- "age"
  sim$age_strat <- do.call(strat_means_fun, strat_args)

  ## Return results
  sim

}



#' Calculate error of stratified estimates
#'
#' @param sim Object returned by [`run_strat()`], which includes the simulated population,
#' survey results, and stratified analysis outputs.
#'
#' @return Adds error details and summary statistics to the `sim` list, ending with
#' `"*_strat_error"` and `"*_strat_error_stats"`. Error statistics include:
#'
#' - **MAE**: Mean absolute error
#' - **MSE**: Mean squared error
#' - **RMSE**: Root mean squared error
#'
#' @examples
#' \donttest{
#' sim <- sim_abundance(ages = 1:5, years = 1:5,
#'                      R = sim_R(log_mean = log(1e+7)),
#'                      growth = sim_vonB(length_group = 1)) |>
#'   sim_distribution(grid = make_grid(res = c(20, 20)),
#'                    ays_covar = sim_ays_covar(sd = 1)) |>
#'   sim_survey(n_sims = 1, q = sim_logistic(k = 2, x0 = 3)) |>
#'   run_strat() |>
#'   strat_error()
#' }
#'
#' @export

strat_error <- function(sim) {

  total <- age <- NULL

  ## total_strat
  I_hat <- sim$total_strat[, list(sim, year, total)]
  names(I_hat) <- c("sim", "year", "I_hat")
  I <- data.frame(year = sim$years, I = colSums(sim$I))
  comp <- merge(I_hat, I, by = "year")
  comp$error <- comp$I_hat - comp$I
  means <- error_stats(comp$error)
  sim$total_strat_error <- comp
  sim$total_strat_error_stats <- means

  ## length_strat
  I_hat <- sim$length_strat[, list(sim, year, length, total)]
  names(I_hat) <- c("sim", "year", "length", "I_hat")
  sly <- expand.grid(sim = seq(max(sim$total_strat$sim)),
                     year = sim$years, length = sim$lengths)
  I_hat <- merge(sly, I_hat, by = c("sim", "year", "length"), all = TRUE) # expand to all lengths
  I_hat$I_hat[is.na(I_hat$I_hat)] <- 0                                    # fill missing lengths with zero
  I <- as.data.frame.table(sim$I_at_length, responseName = "I")
  I$year <- as.numeric(as.character(I$year))
  I$length <- as.numeric(as.character(I$length))
  comp <- merge(data.table(I_hat), data.table(I), by = c("year", "length"))
  comp$error <- comp$I_hat - comp$I
  means <- error_stats(comp$error)
  sim$length_strat_error <- comp
  sim$length_strat_error_stats <- means

  ## age_strat
  I_hat <- sim$age_strat[, list(sim, year, age, total)]
  names(I_hat) <- c("sim", "year", "age", "I_hat")
  say <- expand.grid(sim = seq(max(sim$total_strat$sim)),
                     year = sim$years, age = sim$ages)
  I_hat <- merge(say, I_hat, by = c("sim", "year", "age"), all = TRUE) # expand to all ages
  I_hat$I_hat[is.na(I_hat$I_hat)] <- 0                                 # fill missing ages with zero
  I <- as.data.frame.table(sim$I, responseName = "I")
  I$year <- as.numeric(as.character(I$year))
  I$age <- as.numeric(as.character(I$age))
  comp <- merge(data.table(I_hat), data.table(I), by = c("year", "age"))
  comp$error <- comp$I_hat - comp$I
  means <- error_stats(comp$error)
  sim$age_strat_error <- comp
  sim$age_strat_error_stats <- means

  sim

}





