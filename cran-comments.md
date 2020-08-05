## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Paul Regular <Paul.Regular@dfo-mpo.gc.ca>'
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    al (16:21)
    et (16:18)
  Suggests or Enhances not in mainstream repositories:
    INLA
  
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/stable/

> On ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking examples ... NOTE
  Examples with CPU or elapsed time > 5s
                 user system elapsed
  sim_survey    4.975  0.102   9.224
  strat_error   4.599  0.037   9.197
  run_strat     4.556  0.078  11.762
  sim_abundance 3.391  0.228   7.234
  make_grid     3.151  0.014   5.905
  ** found \donttest examples: check also with --run-donttest

0 errors √ | 0 warnings √ | 2 notes x
