
library(SimSurvey)

## Simulate cod-like population
## See "imitate_cod_data.R" file for details on parameter choices

set.seed(438)
pop <- sim_abundance(ages = 1:20,
                     years = 1:20,
                     R = sim_R(mean = 100000000,
                               log_sd = 0.5,
                               random_walk = TRUE),
                     Z = sim_Z(mean = 0.5,
                               log_sd = 0.2,
                               phi_age = 0.9,
                               phi_year = 0.5)) %>%
  sim_distribution(grid = sim_grid(x_range = c(-140, 140),
                                   y_range = c(-140, 140),
                                   res = c(3.5, 3.5),
                                   shelf_depth = 200,
                                   shelf_width = 100,
                                   depth_range = c(0, 1000),
                                   n_div = 1,
                                   strat_breaks = seq(0, 1000, by = 20),
                                   strat_splits = 2),
                   space_covar = sim_sp_covar(range = 40,
                                              sd = 0.1),
                   ay_covar = sim_ay_covar(sd = 10,
                                           phi_age = 0.1,
                                           phi_year = 0.8,
                                           group_ages = 5:20),
                   depth_par = sim_parabola(mu = 250,
                                            sigma = 50))

## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index
setMKLthreads(1) # turn off MKL hyperthreading
res <- test_surveys(pop,
                    # surveys = expand_surveys(set_den = c(0.3, 0.5, 0.8,
                    #                                      1, 2, 3, 6, 9) / 1000,
                    #                          lengths_cap = c(2, 3, 5, 8, 10, 20,
                    #                                          30, 60, 90, 100, 200,
                    #                                          400, 600, 1000),
                    #                          ages_cap = c(2, 3, 5, 8, 10,
                    #                                       20, 30, 60)),
                    surveys = expand_surveys(set_den = c(0.5, 1, 3, 9) / 1000,
                                             lengths_cap = c(3, 8, 20, 60, 100,
                                                             400, 1000),
                                             ages_cap = c(3, 8, 20, 60)),
                    n_sims = 5,
                    n_loops = 50,
                    cores = 6,
                    q = sim_logistic(k = 2, x0 = 3),
                    growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1, digits = 0),
                    export = "analysis/cod_sim_exports")
# res <- resume_test(dir = "analysis/2018-05-31_test")
setMKLthreads() # turn hyperthreading on again


library(plotly)
sim <- res
d <- merge(sim$surveys, sim$age_strat_error_stats, by = "survey")
z1 <- xtabs(RMSE ~ ages_cap + lengths_cap, data = d, subset = d$set_den == 3e-03)
z2 <- xtabs(RMSE ~ ages_cap + lengths_cap, data = d, subset = d$set_den == 9e-03)
z3 <- xtabs(RMSE ~ ages_cap + lengths_cap, data = d, subset = d$set_den == 1e-03)
z4 <- xtabs(RMSE ~ ages_cap + lengths_cap, data = d, subset = d$set_den == 5e-04)
x <- as.numeric(colnames(z1))
y <- as.numeric(rownames(z1))
plot_ly(x = ~x, y = ~y, showscale = FALSE) %>%
  add_surface(z = ~z1) %>%
  add_surface(z = ~z2) %>%
  add_surface(z = ~z3) %>%
  add_surface(z = ~z4)

d <- merge(sim$surveys, sim$age_strat_error, by = "survey")
x1 <- d$error[d$lengths_cap == 20 & d$ages_cap == 8 & d$age == 3]
x2 <- d$error[d$lengths_cap == 1000 & d$ages_cap == 60 & d$age == 3]
hist(x1, breaks = 100, xlim = range(c(x1, x2)))
hist(x2, breaks = 100, xlim = range(c(x1, x2)))

d <- merge(sim$surveys, sim$age_strat_error, by = "survey")
d <- d[, list(RMSE = SimSurvey::error_stats(error)["RMSE"]),
       by = c("survey", "set_den", "lengths_cap", "ages_cap", "age")]
plot_ly(data = d[d$set_den == 5e-04 & d$age == 3, ], x = ~ages_cap, y = ~lengths_cap, z = ~RMSE) %>%
  add_heatmap()


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
                     lengths_cap = 1000,
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

