## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub fedora-clang-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)

## R CMD check results
> On windows-x86_64-devel (r-devel), fedora-clang-devel (r-devel), ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Paul Regular <Paul.Regular@dfo-mpo.gc.ca>'
  
  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/stable/

> On windows-x86_64-devel (r-devel)
  checking package dependencies ... NOTE
  Package suggested but not available for checking: 'INLA'

> On fedora-clang-devel (r-devel), ubuntu-gcc-release (r-release)
  checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                 user system elapsed
  sim_survey    4.029  0.084   8.099
  run_strat     3.401  0.113   7.143
  sim_abundance 3.189  0.151   6.778
  strat_error   2.532  0.064   5.030
  ** found \donttest examples: check also with --run-donttest

0 errors √ | 0 warnings √ | 3 notes x
