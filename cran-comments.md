## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  
  
  
  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/stable/
    et (17:18)
  Possibly mis-spelled words in DESCRIPTION:
  Maintainer: 'Paul Regular <Paul.Regular@dfo-mpo.gc.ca>'
  New submission
    al (17:21)

> On ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking examples ... NOTE
  Examples with CPU or elapsed time > 5s
              user system elapsed
  sim_survey 7.481  0.085   7.488

0 errors √ | 0 warnings √ | 2 notes x
