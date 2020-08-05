
## pkgdown::build_site(examples = FALSE)
## devtools::check(env_vars = c(NOT_CRAN = "true", R_CHECK_DONTTEST_EXAMPLES = "false"))
## devtools::check_win_devel()
## cran_check <- rhub::check_for_cran(env_vars = c(`_R_CHECK_DONTTEST_EXAMPLES_` = "false", c(`_R_CHECK_FORCE_SUGGESTS_` = "true", `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true")))
## cran_check$cran_summary()

release_questions <- function() {
  c(
    "Did you update the pkgdown site?",
    "Did you update cran-coments?"
  )
}
