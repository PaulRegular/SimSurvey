
library(SimSurvey)
library(data.table)

## Simulate cod-like population
## See "imitate_cod_data.R" file for details on parameter choices

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
                     growth = sim_vonB(Linf = 120, L0 = 5, K = 0.1,
                                       log_sd = 0.1, length_group = 3,
                                       digits = 0)) %>%
  sim_distribution(grid = make_grid(x_range = c(-140, 140),
                                    y_range = c(-140, 140),
                                    res = c(3.5, 3.5),
                                    shelf_depth = 200,
                                    shelf_width = 100,
                                    depth_range = c(0, 1000),
                                    n_div = 1,
                                    strat_breaks = seq(0, 1000, by = 40),
                                    strat_splits = 2),
                   ays_covar = sim_ays_covar(sd = 2.8,
                                             range = 300,
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 5:20),
                   depth_par = sim_parabola(mu = 200,
                                            sigma = 70))



## ~ current sampling protocol
base_case <- sim_survey_parallel(pop, n_sims = 5, n_loops = 5, cores = 6,
                                 set_den = 2e-03, lengths_cap = 500, ages_cap = 10,
                                 q = sim_logistic(k = 2, x0 = 3),
                                 quiet = FALSE)
setdet <- base_case$setdet[, list(year, x, y, depth, strat, n, sim, set)]
samp <- base_case$samp[base_case$samp$measured, ]
samp$age[!samp$aged] <- NA
samp$measured <- samp$aged <- NULL
samp <- merge(setdet[, list(year, sim, set)], samp, by = "set")
# fwrite(setdet, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/01_setdet.csv")
# fwrite(samp, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/01_samp.csv")
# rm(base_case)
# gc()

## One example of a protocol that results in bias
bias_case <- sim_survey_parallel(pop, n_sims = 1, n_loops = 25, cores = 6,
                                 set_den = 1e-02, lengths_cap = 10, ages_cap = 10,
                                 q = sim_logistic(k = 2, x0 = 3),
                                 quiet = FALSE)
setdet <- bias_case$setdet[, list(year, x, y, depth, strat, n, sim, set)]
samp <- bias_case$samp[bias_case$samp$measured, ]
samp$age[!samp$aged] <- NA
samp$measured <- samp$aged <- NULL
samp <- merge(setdet[, list(year, sim, set)], samp, by = "set")
# fwrite(setdet, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/02_setdet.csv")
# fwrite(samp, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/02_samp.csv")
# rm(bias_case)
# gc()


## 2018-12-11 - comps to Geoff's output ----------------------------------------


"
For each set there are three lines:
  the set number, position (x,y,depth), number in the set and number measured
the age composition with spatially varying age-length key
age composition with spatially uniform key

year 4 sim 1
645 [103.25 103.25 226.  ] 516 500
[  1 199 314   1   0   0   0   0   0   0   0   0]
[  8 233 274   2   0   0   0   0   0   0   0   0]

683 [-117.25 -138.25  147.  ] 550 500
[  0   0   0 101 210 120  64  33  12   5   3   0]
[  1   0   0  30 134 123 107  88  38  17  10   0]

666 [-138.25 -127.75   15.  ] 562 500
[  0   0   0   0 291 144  71  36  10   7   3   0]
[  1   0   0   4 192 141 103  69  25  15   8   2]

year 4 sim 3
8270 [-131.25 -113.75   68.  ] 636 500
[  0   0   0   0 309 143  82  50  35  10   6   0]
[  2   0   0   4 187 143 124  67  74  19  11   1]

8190 [ 85.75  96.25 141.  ] 831 500
[  8 246 565  11   0   0   0   0   0   0   0   0]
[ 30 318 469  14   0   0   0   0   0   0   0   0]

8186 [ 89.25  96.25 148.  ] 1040 500
[ 17 253 743  28   0   0   0   0   0   0   0   0]
[ 36 390 597  17   0   0   0   0   0   0   0   0]

year 7 sim 1
1214 [103.25  33.25 226.  ] 479 479
[  0  14 464   0   0   0   0   0   0   0   0   0]
[  0  36 442   0   0   0   0   0   0   0   0   0]

1276 [-106.75 -106.75  184.  ] 802 500
[  0   1   0  36 323 255  95  38  20  20   7   2]
[  2   4   0  18 226 215 120  91  44  50  12  12]

1260 [-113.75 -138.25  161.  ] 1014 500
[  0   2   0 209 266 244 141 104  20  17   8   1]
[  1   6   0  87 221 218 154 172  55  56  26  10]

year 7 sim 3
8819 [110.25  64.75 303.  ] 731 500
[  0  46 685   0   0   0   0   0   0   0   0   0]
[  0  38 692   0   0   0   0   0   0   0   0   0]

8817 [106.75  64.75 261.  ] 837 500
[  0   4 832   0   0   0   0   0   0   0   0   0]
[  0  44 793   0   0   0   0   0   0   0   0   0]

8818 [110.25  68.25 303.  ] 1056 500
[  0  66 990   0   0   0   0   0   0   0   0   0]
[   0   56 1000    0    0    0    0    0    0    0    0    0]
"

## Function for comparing numbers
comp_fun <- function(x_y_year = c(103.25, 103.25, 4),
                     mod_N = NULL,
                     sp_N = NULL) {

  I <- sp_N[x == x_y_year[1] & y == x_y_year[2] & year == x_y_year[3], list(age, I)]
  I$prop <- I$I / sum(I$I)
  true <- I$prop[I$age %in% 1:12] * 100
  true_mat <- t(replicate(2, true))
  est <- fread(mod_N)
  est <- unname(as.matrix(est))
  est_prop <- prop.table(est, margin = 1) * 100
  dif <- true_mat - est_prop

  blank <- rep(NA, ncol(est))
  comp_tab <- rbind(est,
                    blank,
                    est_prop,
                    blank,
                    true,
                    blank,
                    dif)
  rownames(comp_tab) <- c("Est", "", "", "Est prop", "", "", "True prop", "", "Diff", "")
  colnames(comp_tab) <- 1:12

  options(knitr.kable.NA = "")
  cat("x = ", x_y_year[1], " | y = ", x_y_year[2], " | year = ", x_y_year[3])
  knitr::kable(comp_tab, digits = 0)

}


sp_N <- merge(base_case$grid_xy, base_case$sp_N, by = "cell")
q_fun <- sim_logistic(k = 2, x0 = 3)
sp_N$I <- sp_N$N * q_fun(sp_N$age)


comp_fun(x_y_year = c(103.25, 103.25, 4),
         sp_N = sp_N,
         mod_N = "1 199 314   1   0   0   0   0   0   0   0   0
                  8 233 274   2   0   0   0   0   0   0   0   0")

comp_fun(x_y_year = c(-117.25, -138.25, 4),
         sp_N = sp_N,
         mod_N = "0   0   0 101 210 120  64  33  12   5   3   0
                  1   0   0  30 134 123 107  88  38  17  10   0")

comp_fun(x_y_year = c(-138.25, -127.75, 4),
         sp_N = sp_N,
         mod_N = "0   0   0   0 291 144  71  36  10   7   3   0
                  1   0   0   4 192 141 103  69  25  15   8   2")

comp_fun(x_y_year = c(-138.25, -113.75, 4),
         sp_N = sp_N,
         mod_N = "0   0   0   0 309 143  82  50  35  10   6   0
                  2   0   0   4 187 143 124  67  74  19  11   1")

comp_fun(x_y_year = c(85.75, 96.25, 4),
         sp_N = sp_N,
         mod_N = " 8 246 565  11   0   0   0   0   0   0   0   0
                  30 318 469  14   0   0   0   0   0   0   0   0")

comp_fun(x_y_year = c(89.25, 96.25, 4),
         sp_N = sp_N,
         mod_N = "17 253 743  28   0   0   0   0   0   0   0   0
                  36 390 597  17   0   0   0   0   0   0   0   0")

comp_fun(x_y_year = c(103.25, 33.25, 7),
         sp_N = sp_N,
         mod_N = "0  14 464   0   0   0   0   0   0   0   0   0
                  0  36 442   0   0   0   0   0   0   0   0   0")

comp_fun(x_y_year = c(-106.75, -106.75, 7),
         sp_N = sp_N,
         mod_N = "0   1   0  36 323 255  95  38  20  20   7   2
                  2   4   0  18 226 215 120  91  44  50  12  12")

comp_fun(x_y_year = c(-113.75, -138.25, 7),
         sp_N = sp_N,
         mod_N = "0   2   0 209 266 244 141 104  20  17   8   1
                  1   6   0  87 221 218 154 172  55  56  26  10")

comp_fun(x_y_year = c(110.25, 64.75, 7),
         sp_N = sp_N,
         mod_N = "0  46 685   0   0   0   0   0   0   0   0   0
                  0  38 692   0   0   0   0   0   0   0   0   0")

comp_fun(x_y_year = c(106.75, 64.75, 7),
         sp_N = sp_N,
         mod_N = "0   4 832   0   0   0   0   0   0   0   0   0
                  0  44 793   0   0   0   0   0   0   0   0   0")

comp_fun(x_y_year = c(110.25, 68.25, 7),
         sp_N = sp_N,
         mod_N = "0   66  990    0    0    0    0    0    0    0    0    0
                  0   56 1000    0    0    0    0    0    0    0    0    0")



## 2018-12-17 - real 3Ps data --------------------------------------------------

library(Rstrap)

## Load survey data compiled for Rstrap
load("analysis/rv_data/converted_set_details_2018-02-26.Rdata")
load("analysis/rv_data/converted_length-frequency_data_2018-02-26.Rdata")
load("analysis/rv_data/age-growth_data_2018-02-26.Rdata")

## Subset data to cod
con.setdet <- con.setdet[(con.setdet$rec == 5 |
                            (con.setdet$rec == 6 & con.setdet$spec == 438)) &
                           con.setdet$NAFOdiv == "3P", ]
con.lf <- con.lf[con.lf$spec == 438 & con.lf$NAFOdiv == "3P", ]
ag <- ag[ag$spec == 438 & ag$NAFOdiv == "3P", ]
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

fwrite(out$raw.data$set.details, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/cod_3Ps_setdet.csv")
fwrite(out$raw.data$age.growth, file = "analysis/cod_sim_exports/2018-09-12_for_Geoff/cod_3Ps_samp.csv")


