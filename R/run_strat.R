
#' Calculate common error statistics
#'
#' @param error Vector of errors
#'
#' @return Returns a named vector of error statistics including mean absolute
#'         error (\code{"MAE"}), mean squared error (\code{"MSE"}), and
#'         root mean squared error (\code{"RMSE"})
#'
#' @export
#'

error_stats <- function(error) {
  c(MAE = mean(abs(error)),
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
#'
#' @export
#'

strat_data <- function(sim, length_group = 3) {

  ## Extract setdet and samp objects, and add sim and year to samp data
  setdet <- sim$setdet
  samp <- sim$samp
  samp <- merge(setdet[, list(set, sim, year)], samp, by = "set")

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
  cj <- merge(setdet[, list(set, sim, year)], cj, by = "set", all = TRUE)
  lf <- merge(cj, lf, by = c("set", "length"), all = TRUE)
  lf$length_freq[is.na(lf$length_freq)] <- 0 # replace NA's with 0's

  ## Construct age-length key by sim and year
  ag <- samp
  ag <- ag[ag$aged, ] # discard records without corresponding length and age values
  ag$length <- group_lengths(ag$length, length_group)
  alk <- ag[, list(age_freq = .N), by = c("sim", "year", "length", "age")] # combined alk
  alk[, age_tot := sum(age_freq), by = c("sim", "year", "length")]
  alk[, age_prop := age_freq / age_tot, by = c("sim", "year", "length")]
  setkeyv(alk, c("sim", "year", "length"))

  ## Merge lf and alk objects
  alf <- merge(lf, alk, by = c("sim", "year", "length"), all = TRUE, allow.cartesian = TRUE)
  alf <- alf[!is.na(length_freq)] # drop matches lacking length frequencies (further ensure bins not observed in the sets are not represented in the results)
  alf[, age_freq := ifelse(is.na(age), length_freq, length_freq * age_prop)] # pass 100% len.freq if age is unknown for that length bin, otherwise apply proportion
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
#' It first calculates strat-level statistics and then the larger-scale stistics like total abundance
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

  ## set-up
  lc <- (100 - confidence) / 200
  uc <- (100 - confidence) / 200 + (confidence / 100)
  d <- copy(data) # make a copy of the provided data.table
  d <- d[,  c(strat_groups,  metric),  with = FALSE]
  setnames(d, names(d), c(strat_groups, "metric"))
  setkeyv(d,  strat_groups)

  ## strat.tab includes strat-level means,  variances,  etc.
  strat_tab <- d[, list(sumYh = sum(metric), meanYh = mean(metric),
                        varYh = var(metric), nh = .N), by = strat_groups]
  strat_tab[, Nh := strat_area / tow_area] # number of sample units in each strat
  strat_tab[, Wh := Nh / sum(Nh), by = survey_groups]
  strat_tab[, total := Nh * sumYh / nh]
  strat_tab[, gh := Nh * (Nh - nh) / nh]

  ## survey.tab includes large-scale means,  such as mean weight per tow,  and totals,  such as total biomass + associated confidence intervals
  survey_tab <- strat_tab[, list(n = sum(nh), N = sum(Nh), meanYst = sum(Wh * meanYh),
                                 varYst = (1 / ((sum(Nh)) ^ 2)) * sum(gh * varYh),
                                 df = ((sum(gh * varYh)) ^ 2) / (sum((gh ^ 2 * varYh ^ 2) / (nh - 1)))),
                          by = survey_groups]
  survey_tab[, meanYst_lcl := (meanYst - (sqrt(varYst)) * abs(qt(lc,  df)))]
  survey_tab[, meanYst_ucl := (meanYst + (sqrt(varYst)) * abs(qt(lc,  df)))]
  survey_tab[, sumYst := N * meanYst]
  survey_tab[, sumYst_lcl := (sumYst - abs(qt(lc, df)) * N * sqrt(varYst))]
  survey_tab[, sumYst_ucl := (sumYst + abs(qt(lc, df)) * N * sqrt(varYst))]
  survey_tab[sapply(survey_tab,  is.nan)] <- NA

  ## Rename cols
  survey_tab <- survey_tab[, c(survey_groups,  "n", "N", "meanYst",  "meanYst_lcl", "meanYst_ucl",
                               "sumYst", "sumYst_lcl", "sumYst_ucl"),  with = FALSE]
  setnames(survey_tab, names(survey_tab), c(survey_groups, "sets", "sampling_units",
                                            "mean", "mean_lcl", "mean_ucl",
                                            "total", "total_lcl", "total_ucl"))

  ## return results (here I only return the survey_tab to minimize details and object size)
  survey_tab

}


#' Run stratified analysis on simulated data
#'
#' @param sim    Simulation from \code{\link{sim_survey}}
#' @param length_group  Size of the length frequency bins
#'
#' @return Adds stratified analysis results for the total population (\code{"total_strat"})
#'         and the population aggregated by length group and age (\code{"length_strat"} and
#'         \code{"age_strat"}, respectively) to the \code{sim} list.
#'
#' @export
#'

run_strat <- function(sim, length_group = 3) {

  ## Prep data for strat_means function
  data <- strat_data(sim, length_group = length_group)
  data$setdet <- data$setdet[, c("sim", "year", "division", "strat", "strat_area", "tow_area", "set", "n"), with = FALSE]
  data$lf <- merge(data$setdet[, setdiff(names(data$setdet), "n"), with = FALSE], data$lf, by = "set", all = TRUE)
  data$af <- merge(data$setdet[, setdiff(names(data$setdet), "n"), with = FALSE], data$af, by = "set", all = TRUE)

  ## Run strat calculations
  strat_args <- list(data = data$setdet,
                     metric = "n",
                     strat_groups = c("sim", "year", "division", "strat", "strat_area", "tow_area"),
                     survey_groups = c("sim", "year"))
  sim$total_strat <- do.call(strat_means, strat_args)

  strat_args$data <- data$lf
  strat_args$strat_groups <- c(strat_args$strat_groups, "length")
  strat_args$survey_groups <- c(strat_args$survey_groups, "length")
  sim$length_strat <- do.call(strat_means, strat_args)

  strat_args$data <- data$af
  strat_args$strat_groups[strat_args$strat_groups == "length"] <- "age"
  strat_args$survey_groups[strat_args$survey_groups == "length"] <- "age"
  sim$age_strat <- do.call(strat_means, strat_args)

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
#' @export
#'

strat_error <- function(sim) {

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





