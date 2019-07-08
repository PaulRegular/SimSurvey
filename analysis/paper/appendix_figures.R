
library(raster)
library(SimSurvey)
library(Rstrap)
library(plotly)
library(dplyr)

source("analysis/paper/helper_functions.R")

## REAL DATA ----------

## Load survey data compiled for Rstrap
load("analysis/rv_data/converted_set_details_2018-02-26.Rdata")
load("analysis/rv_data/converted_length-frequency_data_2018-02-26.Rdata")
load("analysis/rv_data/age-growth_data_2018-02-26.Rdata")

## Subset data to cod
con.setdet <- con.setdet[(con.setdet$rec == 5 | (con.setdet$rec == 6 & con.setdet$spec == 438)), ]
con.lf <- con.lf[con.lf$spec == 438, ]
ag <- ag[ag$spec == 438, ]
rv_data <- list(setdet = con.setdet, lf = con.lf, ag = ag)

## Save Rstrap output
index.strata <- c(293:300, 306:326, 705:708, 711:716, 779:783)
out <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = rv_data$ag,
                 data.series = "Campelen", program = "strat2 & strat1",
                 species = 438, survey.year = 1995:2017, season = "spring",
                 NAFOdiv = "3P", strat = c(293:300, 306:326, 705:708, 711:716, 779:783),
                 sex = c("male","female","unsexed"), length.group = 3,
                 group.by = "length & age",
                 export = NULL, plot.results = FALSE)

## Convert lat and lon to UTM
setdet <- data.table(out$raw.data$set.details)
xy <- cbind(-setdet$long.start, setdet$lat.start) %>%
  rgdal::project(., proj = sp::proj4string(survey_grid)) %>%
  data.frame(.)
names(xy) <- c("easting", "northing")
setdet <- cbind(setdet, xy)

## Make a unique identifier
setdet$id <- with(setdet, paste(survey.year, vessel, trip, set))
ag <- out$raw.data$age.growth
ag$id <- with(ag, paste(survey.year, vessel, trip, set))

## Melt age frequency data
af <- data.table::melt(setdet,
                       id.vars = c("survey.year", "vessel", "trip", "set", "easting", "northing"),
                       measure.vars = names(setdet)[grepl("^af", names(setdet))],
                       variable.name = "age", value.name = "freq")
af <- af[af$age != "afNA", ]
af$age <- as.integer(gsub("af", "", af$age))

## Melt length frequency data
idvar <- c("survey.year", "vessel", "trip", "set")
sub_lf <- merge(out$raw.data$set.details[, idvar], con.lf, by = idvar)
sub_lf$id <- with(sub_lf, paste(survey.year, vessel, trip, set))
ind <- grep("^bin|^id$", names(sub_lf)) # get length frequencies and expand to samples
lf <- sub_lf[, ind]
len_samp <- data.table::melt(lf, id.vars = "id")
len_samp <- len_samp[len_samp$value > 0, ]
len_samp$value <- round(len_samp$value)
len_samp$length <- as.numeric(gsub("bin", "", len_samp$variable))
len_samp <- len_samp[rep(row.names(len_samp), len_samp$value), c("id", "length")]
len_samp <- data.table(len_samp)


## SIMULATED DATA ----------

set.seed(438)
pop <- sim_abundance(ages = 1:20,
                     years = 1:20,
                     R = sim_R(mean = 30000000,
                               log_sd = 0.5,
                               random_walk = TRUE),
                     Z = sim_Z(mean = 0.5,
                               log_sd = 0.2,
                               phi_age = 0.9,
                               phi_year = 0.5),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0)) %>%
  sim_distribution(grid = make_grid(x_range = c(-140, 140),
                                   y_range = c(-140, 140),
                                   res = c(3.5, 3.5),
                                   shelf_depth = 200,
                                   shelf_width = 100,
                                   depth_range = c(0, 1000),
                                   n_div = 1,
                                   strat_breaks = seq(0, 1000, by = 40),
                                   strat_splits = 2,
                                   method = "spline"),
                   ays_covar = sim_ays_covar(sd = 2.8,
                                             range = 300,
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 5:20),
                   depth_par = sim_parabola(mu = 200,
                                            sigma = 70))

survey <- sim_survey(pop,
                     n_sims = 1,
                     light = FALSE,
                     set_den = 3 / 1000,
                     lengths_cap = 400,
                     ages_cap = 10,
                     q = sim_logistic(k = 2, x0 = 3))


## COMPS ----------

set.seed(12)

ax <- list(
  zeroline = FALSE,
  showline = TRUE,
  showgrid = FALSE,
  mirror = "ticks",
  showticklabels = FALSE,
  title = ""
)

no_ax <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)

yax <- ax
yax$title <- ""

real <- setdet %>%
  filter(survey.year == sample(c(1997:2005, 2007:2017), 1)) %>%
  mutate(zero = number == 0) %>%
  arrange(-number) %>%
  mutate(rank = seq(1, nrow(.)))

sub_real <- sample(real$id[real$number > 50], 3)

sel <- as.numeric(factor(real$id, levels = sub_real))
sel[is.na(sel)] <- 0
real$sel <- sel

real_hist <- real %>%
  plot_ly(x = ~number) %>%
  add_histogram(color = I("#999999")) %>%
  layout(yaxis = list(title = ""),
         xaxis = list(title = ""))

real_dist <- real %>%
  plot_ly(x = ~easting, y = ~northing, size = ~number,
          sizes = c(10, 500), symbol = ~zero,
          symbols = c(16, 4), color = ~factor(sel),
          colors = c("#999999", viridis::viridis(3))) %>%
  add_markers() %>%
  hide_legend() %>%
  layout(xaxis = ax, yaxis = yax)

simulated <- survey$setdet %>%
  filter(year == sample(1:20, 1)) %>%
  mutate(zero = n == 0) %>%
  arrange(-n) %>%
  mutate(rank = seq(1, nrow(.)))

sub_sim <- sample(simulated$set[simulated$n > 50], 3)

sel <- as.numeric(factor(simulated$set, levels = sub_sim))
sel[is.na(sel)] <- 0
simulated$sel <- sel

sim_hist <- simulated %>%
  plot_ly(x = ~n) %>%
  add_histogram(color = I("#999999")) %>%
  layout(yaxis = list(title = "Frequency"),
         xaxis = list(title = ""))

sim_dist <- simulated %>%
  plot_ly(x = ~x, y = ~y, size = ~n,
          sizes = c(10, 500), symbol = ~zero,
          symbols = c(16, 4), color = ~factor(sel),
          colors = c("#999999", viridis::viridis(3))) %>%
  add_markers() %>%
  hide_legend() %>%
  layout(xaxis = ax, yaxis = ax)


empty_plot <- plot_ly(x = 1, y = 1) %>%
  add_text(text = "Numbers caught", size = I(13)) %>%
  layout(xaxis = no_ax, yaxis = no_ax)

empty_row <- subplot(empty_plot, empty_plot) # hack to add x axis labels (margins were mis-behaving using the standard approach)
hist_row <- subplot(real_hist, sim_hist, shareY = TRUE, titleX = TRUE)
dist_row <- subplot(real_dist, sim_dist, titleY = TRUE)


real_samps <- len_samp %>%
  filter(id %in% sub_real) %>%
  mutate(sel = factor(id, levels = sub_real))

real_ag <- ag %>%
  filter(id %in% sub_real) %>%
  mutate(sel = factor(id, levels = sub_real))

real_len <- list()
cols <- viridis::viridis(3)
for (i in 1:3) {

  l <- real_samps %>%
    filter(as.numeric(sel) == i) %>%
    group_by(length) %>%
    summarise(n = n())

  a <- real_ag %>%
    filter(as.numeric(sel) == i) %>%
    group_by(length) %>%
    summarise(n = n()) %>%
    rename(n_aged = n)

  s <- left_join(l, a, by = "length")

  real_len[[i]] <- s %>%
    plot_ly(x = ~length, color = I(cols[i])) %>%
    add_bars(y = ~n_aged, alpha = 1) %>%
    add_bars(y = ~n, alpha = 0.5) %>%
    layout(barmode = "stack")

}

real_len <- subplot(real_len, nrows = 3, shareX = TRUE, margin = 0.05) %>%
  hide_legend() %>%
  layout(yaxis2 = list(title = ""),
         xaxis = list(title = "Length (cm)"))


sim_samps <- survey$samp %>%
  filter(set %in% sub_sim) %>%
  mutate(sel = factor(set, levels = sub_sim))


sim_len <- list()
cols <- viridis::viridis(3)
for (i in 1:3) {
  sim_len[[i]] <- sim_samps %>%
    filter(as.numeric(sel) == i) %>%
    group_by(length) %>%
    summarise(n = sum(measured), n_aged = sum(aged)) %>%
    plot_ly(x = ~length, color = I(cols[i])) %>%
    add_bars(y = ~n_aged, alpha = 1) %>%
    add_bars(y = ~n, alpha = 0.5) %>%
    layout(barmode = "stack")
}

sim_len <- subplot(sim_len, nrows = 3, shareX = TRUE, margin = 0.05) %>%
  hide_legend() %>%
  layout(yaxis2 = list(title = ""),
         xaxis = list(title = "Length (cm)"))


len_row <- subplot(real_len, sim_len, shareX = TRUE, titleY = TRUE, margin = 0.05) %>%
  layout(xaxis = list(range = range(c(real_samps$length, sim_samps$length))),
         xaxis2 = list(range = range(c(real_samps$length, sim_samps$length))))


p <- subplot(hist_row,
             empty_row,
             dist_row,
             len_row,
             nrows = 4,
             heights = c(0.14, 0.06, 0.4, 0.4),
             titleY = TRUE,
             titleX = TRUE) %>%
  add_labels("a) real data", x = 0, y = 1.04) %>%
  add_labels("b) simulated data", x = 0.8, y = 1.04)

p

export_plot(p, file = "analysis/paper/figures/real_sim_comp.png", width = 500, height = 650)

