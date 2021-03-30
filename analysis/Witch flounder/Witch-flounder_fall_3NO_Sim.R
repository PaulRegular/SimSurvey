# Simulate Fall 3NO witch flounder population
## See "imitate_witch_fall_3NO.R" file for details on parameter choices
# #-------------------------------------------------
# # R code for automatically loading and installing required packages.
#
# # Set repository options:
# local({r <- getOption("repos")
# r["CRAN"] <- "http://cran.stat.sfu.ca/"
# options(repos=r)})
#
#
# pkg_list2 = c("iterators","tidyr","crayon","data.table","foreach","doParallel","magrittr","dplyr",
#               "withr","ggplot2","plotly","progress")
# for (pkg2 in pkg_list2)
# {
# # # Try loading the library.
# if ( ! library(pkg2, logical.return=TRUE, character.only=TRUE,lib.loc = "Rlibs") )
# {
#   # # If the library cannot be loaded, install it; then load.
#   install.packages(pkg2,lib = "Rlibs" )
#   library(pkg2,logical.return=TRUE, character.only=TRUE,lib.loc= "Rlibs")
# }
# }
# # Set a vector of strings: package names to use (and install, if necessary)
# pkg_list = c("SimSurvey")
#
# for (pkg in pkg_list)
# {
#   library(pkg, logical.return=TRUE, character.only=TRUE,lib.loc = "Rlibs")
# }
# #install.packages("dlm", lib = "/usr/share/R/library")
###------------------------------
library(SimSurvey)
set.seed(890)
pop <- sim_abundance(ages = 1:16,
                     years = 1:24,
                     R = sim_R(log_mean = log(c(1500000,800000,2000000,4000000,
                                                2000000,1000000,3000000,3000000,5000000,
                                                5000000,5000000,5000000,18000000,22000000,
                                                8500000,5000000,8000000,5000000,3500000,
                                                2200000,3500000,3500000,6500000,10000000)),
                               log_sd = 0.7,
                               random_walk = T),
                     Z = sim_Z(log_mean = log(0.08),
                               log_sd = 0.75,
                               phi_age = 0.8,
                               phi_year = 0.9),
                     N0 = sim_N0(N0 = "exp", plot = FALSE),
                     growth = sim_vonB(Linf = 60, L0 = 5,
                                       K = 0.2, log_sd = 0.1,
                                       length_group = 2, digits = 0)) %>%
  sim_distribution(grid = make_grid(x_range = c(-184, 184),
                                    y_range = c(-184, 184),
                                    res = c(3.5, 3.5),
                                    shelf_depth = 60,
                                    shelf_width = 170,
                                    depth_range = c(0, 1600),
                                    n_div = 2,
                                    strat_breaks = seq(0, 1600, by = 65),
                                    strat_splits = 4,
                                    method = "bezier"),
                   ays_covar = sim_ays_covar(sd = 1.5,
                                             range = 800,
                                             lambda = .55,
                                             model = "matern",
                                             phi_age = 0.8,
                                             phi_year = 0.9,
                                             #group_ages = 8:16
                                             ),
                   depth_par = sim_parabola(mu = log(180),
                                            sigma = 0.8,
                                             log_space = TRUE))


## Test a series of surveys
## Simulate surveys and compare stratified estimates to the true index

surveys <- expand_surveys(set_den = c(0.5, 0.8, 1, 1.2, 1.5) / 1000,
                          lengths_cap = c(100, 200, 300, 400, 500),
                          ages_cap = c(8, 16, 24))
surveys[surveys$set_den == 0.001 &
          surveys$lengths_cap == 300 &
          surveys$ages_cap == 16, ]
sim <- test_surveys(pop,
                    surveys = surveys,
                    keep_details = 38,   ## survey 38 ~current protocol
                    n_sims = 5,
                    n_loops = 200,
                    cores = NULL,
                    q = sim_logistic(k = 2, x0 = 2),
                    export_dir = "C:/Users/fhate/Documents/SimSurvey-doc/Jan21/Sim-output")


# ## Visualize Results
# load("C:/Users/fhate/Documents/SimSurvey-doc/Jan21/Sim-output/test_output.RData")
#
plot_total_strat_fan(sim)
plot_length_strat_fan(sim, surveys = 1:75, years = 1:20, lengths = 0.5:144.5)
plot_age_strat_fan(sim, surveys = 1:75, years = 10, ages = 1:26)

plot_error_surface(sim)
plot_survey_rank(sim)


