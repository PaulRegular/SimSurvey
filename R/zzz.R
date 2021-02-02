
# pkgdown::build_site(examples = FALSE)
# devtools::check(env_vars = c(NOT_CRAN = "true", R_CHECK_DONTTEST_EXAMPLES = "false"))
# devtools::check_win_devel(email = "paul.regular@gmail.com")
# cran_check <- rhub::check_for_cran(email = "paul.regular@gmail.com", env_vars = c(`_R_CHECK_DONTTEST_EXAMPLES_` = "false", `_R_CHECK_FORCE_SUGGESTS_` = "false", `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true"))
# # cran_check <- rhub::get_check(c("SimSurvey_0.1.3.tar.gz-5c5c6c6d39094fb69979eff456c8af95", "SimSurvey_0.1.3.tar.gz-d8ce1e5cb7c54a5599a472a18022dfc3", "SimSurvey_0.1.3.tar.gz-a2249a8410f947d590bfe474b61e6f98"))
# cran_check$cran_summary()
# devtools::release()

release_questions <- function() {
  c(
    "Did you update the pkgdown site?",
    "Did you update cran-coments?"
  )
}
