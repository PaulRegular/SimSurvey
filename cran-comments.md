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
    al (17:21)
    et (17:18)
    INLA
  Availability using Additional_repositories specification:
  
    INLA   yes   https://inla.r-inla-download.org/R/stable/
  Suggests or Enhances not in mainstream repositories:

> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
               user system elapsed
  run_strat   11.63   0.36   11.60
  strat_error 10.25   0.46   10.25
  make_grid    5.69   0.08    5.75

0 errors √ | 0 warnings √ | 2 notes x
