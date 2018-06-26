
library(plotly)
load("analysis/cod_sim_exports/2018-06-06_test/test_output.RData")

d <- merge(sim$surveys, sim$age_strat_error_stats, by = "survey")
d <- d[d$set_den >= 0.001, ]
split_d <- split(d, d$set_den)
split_d <- lapply(split_d, function(.) {
  xtabs(RMSE ~ ages_cap + lengths_cap, data = ., subset = NULL)
})
x <- sort(unique(d$lengths_cap))
y <- sort(unique(d$ages_cap))

p <- plot_ly(x = ~x, y = ~y, cmin = min(d$RMSE), cmax = max(d$RMSE))
visible <- showscale <- rep(TRUE, length(split_d))
steps <- list()
showscale <- c(TRUE, rep(FALSE, length(split_d) - 1))
for (i in seq_along(split_d)) {
  vis <- rep(FALSE, length(split_d))
  vis[i] <- TRUE
  p <- p %>% add_surface(z = split_d[[i]],
                         visible = i == 1,
                         showscale = vis,
                         name = names(split_d)[i],
                         colorbar = list(title = "RMSE"))
  steps[[i]] <- list(args = list(list(visible = vis,
                                      showscale = vis)),
                     method = "update",
                     label = names(split_d)[i])
  if (i == length(split_d)) {
    steps[[i + 1]] <- list(args = list(list(visible = rep(TRUE, length(split_d)),
                                            showscale = showscale)),
                           method = "update",
                           label = "all")
  }
}

p %>%
  layout(sliders = list(list(
    currentvalue = list(prefix = "Set density: "),
    steps = steps
  )),
  scene = list(
    xaxis = list(title = "max(lengths)"),
    yaxis = list(title = "max(ages)"),
    zaxis = list(title = "RMSE",
                 range = range(d$RMSE)),
    camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
  ))








d <- merge(sim$surveys, sim$age_strat_error, by = "survey")
sub_d <- d[d$set_den == 5e-04 & d$lengths_cap == 1000 & d$ages_cap == 60, ]
hist(sub_d$error, breaks = 500)
hist(sub_d$error, breaks = 500, xlim = c(-1e+08, 1e+08))

d <- merge(sim$surveys, sim$age_strat_error, by = "survey")
d %>%
  filter(age == 3 & set_den == 5e-04 & lengths_cap == 1000 & ages_cap == 60) %>%
  group_by(sim) %>%
  plot_ly() %>%
  add_lines(x = ~year, y = ~I_hat, size = I(0.5), alpha = 0.5,
            name = "estimate", color = I("steelblue")) %>%
  add_lines(x = ~unique(year), y = ~unique(I), name = "true", color = I("black"))

d %>%
  filter(year == 14 & set_den == 5e-04 & lengths_cap == 1000 & ages_cap == 60) %>%
  group_by(sim) %>%
  plot_ly() %>%
  add_lines(x = ~age, y = ~I_hat, size = I(0.5), alpha = 0.5,
            name = "estimate", color = I("steelblue")) %>%
  add_lines(x = ~unique(age), y = ~unique(I), name = "true", color = I("black"))

## larger error in the core of the index??

d %>%
  filter(age == 3 & set_den == 5e-04 & lengths_cap == 1000 & ages_cap == 60 &
           error > 500000000)


## Simulate one survey
survey <- sim_survey(pop,
                     n_sims = 1,
                     light = FALSE,
                     set_den = 3 / 1000,
                     lengths_cap = 400,
                     ages_cap = 10,
                     q = sim_logistic(k = 2, x0 = 3),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0)) %>%
  run_strat() %>% strat_error()




survey <- sim_survey(pop,
                     n_sims = 10,
                     light = FALSE,
                     set_den = 10 / 1000,
                     lengths_cap = Inf,
                     ages_cap = Inf,
                     q = sim_logistic(k = 2, x0 = 3),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0),
                     binom_error = TRUE) %>%
  run_strat() %>% strat_error()

hist(survey$age_strat_error[age == 12, ]$error, breaks = 50)


survey <- sim_survey(pop,
                     n_sims = 20,
                     light = FALSE,
                     set_den = 9 / 1000,
                     lengths_cap = 1000,
                     ages_cap = 60,
                     q = sim_logistic(k = 2, x0 = 3),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0),
                     binom_error = TRUE) %>%
  run_strat() %>% strat_error()


## Quick look at distribution
sp_N <- data.frame(merge(survey$sp_N, survey$grid_xy, by = "cell"))
for (j in rev(survey$ages)) {
  z <- xtabs(N ~ x + y, subset = age == j & year == 13, data = sp_N)
  image(z = z, axes = FALSE, col = viridis::viridis(100), main = paste("age", j))
}
for (i in rev(survey$years)) {
  z <- xtabs(N ~ x + y, subset = age == 10 & year == i, data = sp_N)
  image(z = z, axes = FALSE, col = viridis::viridis(100), main = paste("year", i))
}

sim_af <- survey$full_setdet
sim_af %>%
  filter(age == 3 & sim == 1) %>%
  group_by(year) %>%
  plot_ly(x = ~x, y = ~y, size = ~n, frame = ~year,
          text = ~set, sizes = c(5, 500), showlegend = FALSE) %>%
  add_markers() %>%
  animation_opts(frame = 5)


one_samp <- survey$samp[set == 1383, ]
hist(one_samp$length[sample.int(nrow(one_samp), 7000)], breaks = 100)




plot_sets <- function(sim, sim_num = 1) {
  sim$setdet %>%
    filter(sim == sim_num) %>%
    plot_ly(x = ~x, y = ~y, size = ~n+1, frame = ~year,
            text = ~paste("n:", n)) %>%
    add_markers()
}
## The lack of percision at low set densities and high sampling efforts
## may be related to the unbalanced effort in narrow strata??

## May also be related to mismatch between 1 cm length sampling bin and 3 cm strat binning

## Also size up the composite age distribution sampled by the survey at low
## set densities. Composite may be poor at low densities and increased length
## sampling may increase variance because of that




survey <- sim_survey(pop,
                     n_sims = 1,
                     light = FALSE,
                     set_den = 0.5 / 1000,
                     lengths_cap = 20,
                     ages_cap = 10,
                     q = sim_logistic(k = 2, x0 = 3),
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0),
                     binom_error = TRUE) %>%
  run_strat() %>% strat_error()

d <- merge(survey$setdet, survey$samp, by = "set")
d <- d[d$year == 14, ]
true_dist <- obs_dist <- samp_dist <- survey$I[, "14"]
obs_dist[] <- samp_dist[] <- 0
temp <- table(d$age)
obs_dist[names(temp)] <- temp
temp <- table(d$age[d$measured])
samp_dist[names(temp)] <- temp
b <- barplot(true_dist / sum(true_dist), ylim = c(0, 1))
lines(x = b[, 1], y = obs_dist / sum(obs_dist), col = "red")
lines(x = b[, 1], y = samp_dist / sum(samp_dist), col = "blue")

## low set density results != random sample

