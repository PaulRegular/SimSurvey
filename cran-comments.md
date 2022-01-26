## Release summary

- In addition to some minor improvements, this version fixes a failing tests under r-devel ("Error: invalid 'size' argument"). Values supplied to the 'size' argument in sample.int are now length == 1, therefore passing a new error trap in sample.int.

## Test environments

- local R installation, R 4.1.2
- ubuntu 16.04 (on travis-ci), R 4.1.2
- win-builder (devel)
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
  * checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

> On fedora-clang-devel (r-devel)
  * checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
               user system elapsed
  sim_survey  9.644  0.471   1.965
  run_strat   7.009  0.228   1.536
  strat_error 5.742  0.183   0.760
  ** found \donttest examples: check also with --run-donttest
  
> On ubuntu-gcc-release (r-release)
  * checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
             user system elapsed
  run_strat 8.217  0.057   2.038
  ** found \donttest examples: check also with --run-donttest

0 errors √ | 0 warnings √ | 4 notes x



