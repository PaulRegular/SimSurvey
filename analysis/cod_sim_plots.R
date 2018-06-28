
library(SimSurvey)
library(plotly)
library(data.table)

load("analysis/cod_sim_exports/2018-06-27_test/test_output.RData")

vis_sim(sim)




plot_error_surface(sim)

plot_true_vs_est(sim, which_survey = 542, max_sims = 100, facet_by = "age")

plot_true_vs_est(sim, which_survey = 542, max_sims = 100, facet_by = "year")

plot_error_cross_sections(sim)

plot_samp_dist(sim, which_year = 14, which_sim = 4)

