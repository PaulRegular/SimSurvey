
#' Calculate common error statistics
#'
#' @param error Vector of errors
#'
#' @return Returns a named vector of error statistics including
#'         mean error (\code{"ME"}),
#'         mean absolute error (\code{"MAE"}),
#'         mean squared error (\code{"MSE"}) and
#'         root mean squared error (\code{"RMSE"})
#'
#' @export
#'

error_stats <- function(error) {
  c(ME = mean(error),
    MAE = mean(abs(error)),
    MSE = mean(error ^ 2),
    RMSE = sqrt(mean(error ^ 2)))
}


#' Prepare simulated data for stratified analysis
#'
#' @description Generate set details (setdet), length-frequency (lf)
#' and age-frequency (af) data for stratified analysis
#'
#' @param sim           Simulation from \code{\link{sim_survey}}
#' @param length_group  Size of the length frequency bins
#' @param alk_scale     Spatial scale at which to construct and apply age-length-keys:
#'                      "division", "strat" or "set".
#'
#' @export
#'

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



#' Calculate stratified means, variances and confidence intervals across groups
#'
#' @details Function was mainly created for use in the \code{\link{run_strat}} function.
#' It first calculates strat-level statistics and then the larger-scale statistics like total abundance
#'
#' @param data            Expects data.table with all grouping variables in stacked format (must include
#'                        strat_area and tow_area for scaling values)
#' @param metric          Variable in specificed data.table. e.g. "number", "mass"
#' @param strat_groups    Grouping variables for calculations of the fine-scale strat-level
#'                        means (strat and strat_area are required). e.g. c("year", "species",
#'                        "shiptrip", "NAFOdiv", "strat", "strat_area","age")
#' @param survey_groups   Grouping variables for large-scale summary calculations. e.g. ("year","species")
#' @param confidence      Percent for confidence limits
#'
#' @export
#'

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
#' @param sim               Simulation from \code{\link{sim_survey}}
#' @param length_group      Size of the length frequency bins for both abundance at length calculations
#'                          and age-length-key construction. By default this value is inherited from
#'                          the value defined in \code{\link{sim_abundance}} from the closure supplied to
#'                          \code{sim_length} ("inherit"). A numeric value can also be supplied, however,
#'                          a mismatch in length groupings will cause issues with \code{\link{strat_error}}
#'                          as true vs. estimated length groupings will be mismatched.
#' @param alk_scale         Spatial scale at which to construct and apply age-length-keys:
#'                          "division" or "strat".
#' @param strat_data_fun    Function for preparing data for stratified analysis (e.g. \code{\link{strat_data}})
#' @param strat_means_fun   Function for calculating stratified means (e.g. \code{\link{strat_means}})
#'
#' @details The \code{"strat_data_fun"} and \code{"strat_means_fun"} allow the use of custom
#'          \code{\link{strat_data}} and  \code{\link{strat_means}} functions.
#'
#' @return Adds stratified analysis results for the total population (\code{"total_strat"})
#'         and the population aggregated by length group and age (\code{"length_strat"} and
#'         \code{"age_strat"}, respectively) to the \code{sim} list.
#'
#' @examples
#'
#' sim <- sim_abundance(ages = 1:10, years = 1:5,
#'                      R = sim_R(log_mean = log(1e+7)),
#'                      growth = sim_vonB(length_group = 1)) %>%
#'            sim_distribution(grid = make_grid(res = c(10, 10)),
#'                             ays_covar = sim_ays_covar(sd = 1)) %>%
#'            sim_survey(n_sims = 10, q = sim_logistic(k = 2, x0 = 3)) %>%
#'            run_strat()
#'
#' @export
#'

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
#' @param sim   Object from \code{\link{run_strat}} (includes simulated population and
#'              survey along with stratified analysis results)
#'
#' @return Adds details and summary stats of stratified estimate error to the
#'         \code{sim} list, ending with \code{"_strat_error"} or
#'         \code{"_strat_error_stats"}. Error statistics includes mean absolute
#'         error (\code{"MAE"}), mean squared error (\code{"MSE"}), and
#'         root mean squared error (\code{"RMSE"})
#'
#' @examples
#'
#' sim <- sim_abundance(ages = 1:10, years = 1:5,
#'                      R = sim_R(log_mean = log(1e+7)),
#'                      growth = sim_vonB(length_group = 1)) %>%
#'            sim_distribution(grid = make_grid(res = c(10, 10)),
#'                             ays_covar = sim_ays_covar(sd = 1)) %>%
#'            sim_survey(n_sims = 10, q = sim_logistic(k = 2, x0 = 3)) %>%
#'            run_strat() %>%
#'            strat_error()
#'
#' @export
#'

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





