## Release summary

- This patch switches to dependence on the sf and stars package rather than the sp and raster package given the evolution of R spatial. We have also switched to new INLAspacetime functions for simulating from a barrier mesh.

## Test environments

- local R installation, R 4.3.0
- Github Actions ubuntu-latest (release, devel, oldrel-1), 
- GitHub Actions windows-latest (release)
- win-builder (devel)
- R-hub windows-x86_64-devel (r-devel)
- R-hub fedora-clang-devel (r-devel)

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


