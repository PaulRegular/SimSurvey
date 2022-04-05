## Release summary

- In addition to some minor improvements, this patch fixes an error associated with a donttest check, specifically with the use of a logical expression with length > 1. A plotly call was supplied a logical vector when a logical object of length 1 should have been provided. 

## Test environments

- local R installation, R 4.1.2
- Github Actions ubuntu-latest (release, devel, oldrel-1), 
- GitHub Actions windows-latest (release)
- win-builder (devel)
- R-hub windows-x86_64-devel (r-devel)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results

## R CMD check results
> On windows-x86_64-devel (r-devel), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Paul Regular <Paul.Regular@dfo-mpo.gc.ca>'
  
  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/stable/

> On windows-x86_64-devel (r-devel)
  checking package dependencies ... NOTE
  Package suggested but not available for checking: 'INLA'

> On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

> On fedora-clang-devel (r-devel)
  checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
            user system elapsed
  run_strat 6.06  0.275   2.349
  ** found \donttest examples: check also with --run-donttest

0 errors √ | 0 warnings √ | 4 notes x


