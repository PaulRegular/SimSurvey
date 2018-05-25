

#' Convert cut labels to mid points
#'
#' @description Helper function for converting \code{\link{base::cut}} labels
#' into length mid points provided length group (Note: this isn't a general function;
#' output alligns with DFO specific method/labeling)
#'
#' @param cut_labs       Labels from \code{\link{base::cut}}
#' @param length_group   Length group used to cut the length data
#'
#' @export
#'

mid_length <- function(cut_labs, length_group) {
  lr <- strsplit(gsub("[\\[\\)]","",cut_labs),",") # split left and right breaks
  lr <- do.call(rbind, lr)
  l <- as.numeric(lr[, 1]) # grab left breaks
  r <- as.numeric(lr[, 2]) # grab right breaks
  if (length_group == 0.5 | length_group == 1) { m <- l }
  if (length_group > 1) { m <- l + (0.5 * (length_group - 1)) }
  m
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
  setdet <- copy(sim$setdet)
  samp <- copy(sim$samp)
  samp <- merge(setdet[, list(set, sim, year)], samp, by = "set")

  ## Construct length-frequency table
  lf <- copy(samp)
  lf[, n := .N, by = "set"]
  lf <- lf[lf$measured, ]
  lf[, n_lengths := .N, by = "set"]
  length_breaks <- seq(0, max(lf$length, na.rm = TRUE) * 2, length_group)
  lf[, length_group := cut(length, length_breaks, right = FALSE)]
  lf[, ratio := n_lengths / n] # calculate ratio measured
  lf <- lf[, list(length_freq = .N),
           by = c("set", "ratio", "length_group")]
  lf[, length_freq := length_freq / ratio] # scale frequencies
  lf[, ratio := NULL] # discard ratio column (no longer needed)
  setkeyv(lf, c("set", "length_group"))

  ## Add zeros to length-frequency table
  cj <- CJ(set = setdet$set, length_group = levels(lf$length_group), unique = TRUE)
  cj <- merge(setdet[, list(set, sim, year)], cj, by = "set")
  lf <- merge(cj, lf, by = c("set", "length_group"), all = TRUE)
  lf$length_freq[is.na(lf$length_freq)] <- 0 # replace NA's with 0's
  lf$length <- mid_length(lf$length_group, length_group)

  ## Construct age-length key by sim and year
  ag <- copy(samp)
  ag <- ag[ag$aged, ] # discard records without corresponding length and age values
  ag[, length_group := cut(length, length_breaks, right = FALSE)]
  alk <- ag[, list(age_freq = .N), by = c("sim", "year", "length_group", "age")] # combined alk
  alk[, age_tot := sum(age_freq), by = c("sim", "year", "length_group")]
  alk[, age_prop := age_freq / age_tot, by = c("sim", "year", "length_group")]

  ## Merge lf and alk objects
  alf <- merge(lf, alk, by = c("sim", "year", "length_group"), all = TRUE, allow.cartesian = TRUE)
  alf <- alf[!is.na(length_freq)] # drop matches lacking length frequencies (further ensure bins not observed in the sets are not represented in the results)
  alf[, age_freq := ifelse(is.na(age), length_freq, length_freq * age_prop)] # pass 100% len.freq if age is unknown for that length bin, otherwise apply proportion
  af <- alf[, list(age_freq = sum(age_freq, na.rm = TRUE)), by = c("set", "age")]
  setkeyv(af, c("set", "age"))

  ## Add zeros to age-frequency table (this step may be unnecessary)
  ages <- c(NA, seq(min(ag$age, na.rm = TRUE), max(ag$age, na.rm = TRUE)))
  af <- af[CJ(set = setdet$set, age = ages, unique = TRUE), ] # expand to missing groups in each set
  af$age_freq[is.na(af$age_freq)] <- 0 # replace NA's with 0's

  ## Simplify/order and return data
  lf[, length_group := NULL] # no longer needed with midpoint column
  lf <- lf[, c("set", "length", "length_freq"), with = FALSE]
  setkeyv(lf, c("set", "length"))
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
  strat_tab[, totals := Nh * sumYh / nh]
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
  survey_tab <- survey_tab[, c(survey_groups,  "meanYst",  "meanYst_lcl", "meanYst_ucl",
                               "sumYst", "sumYst_lcl", "sumYst_ucl"),  with = FALSE]
  setnames(survey_tab, names(survey_tab), c(survey_groups, "mean", "mean_lcl", "mean_ucl",
                                            "total", "total_lcl", "total_ucl"))

  ## return results (here I only return the survey_tab to minimize details and object size)
  survey_tab

}




#' Run stratified analysis on simulated data
#'
#' @param sim    Simulation from \code{\link{sim_survey}}
#'
#' @export
#'

run_strat <- function(sim) {

  ## Prep data for strat_means function
  data <- strat_data(sim)
  data$setdet <- data$setdet[, c("sim", "year", "division", "strat", "strat_area", "tow_area", "set", "n"), with = FALSE]
  data$lf <- merge(data$setdet[, setdiff(names(data$setdet), "n"), with = FALSE], data$lf, by = "set", all = TRUE)
  data$af <- merge(data$setdet[, setdiff(names(data$setdet), "n"), with = FALSE], data$af, by = "set", all = TRUE)

  ## Run strat calculations
  strat_args <- list(data = data$setdet,
                     metric = "n",
                     strat_groups = c("sim", "year", "division", "strat", "strat_area", "tow_area"),
                     survey_groups = c("sim", "year"))
  strat2 <- do.call(strat_means, strat_args)

  strat1 <- list()
  strat_args$data <- data$lf
  strat_args$strat_groups <- c(strat_args$strat_groups, "length")
  strat_args$survey_groups <- c(strat_args$survey_groups, "length")
  strat1[["length"]] <- do.call(strat_means, strat_args)

  strat_args$data <- data$af
  strat_args$strat_groups[strat_args$strat_groups == "length"] <- "age"
  strat_args$survey_groups[strat_args$survey_groups == "length"] <- "age"
  strat1[["age"]] <- do.call(strat_means, strat_args)

  ## Return results
  list(strat2 = strat2, strat1 = strat1)

}




