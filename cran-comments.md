## Release summary

- This patch switches to dependence on the sf and stars package rather than the sp and raster package given the evolution of R spatial. We have also switched to new INLAspacetime functions for simulating from a barrier mesh.

## Test environments

- local R installation, R 4.3.0
- Github Actions ubuntu-latest (release, devel, oldrel-1), 
- win-builder (devel)
- R-hub windows-x86_64-oldrel (r-oldrel)
- R-hub windows-x86_64-devel (r-devel)

## R CMD check results

❯ On windows-x86_64-devel (r-devel)
  * checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Paul Regular <Paul.Regular@dfo-mpo.gc.ca>'
  
  New submission
  
  Package was archived on CRAN
  
  Possibly misspelled words in DESCRIPTION:
    al (16:21)
    et (16:18)
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2023-08-10 as issues were not corrected
      in time.
  
  Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/stable/
  * checking package dependencies ... NOTE
  Package suggested but not available for checking: 'INLA'

0 errors ✔ | 0 warnings ✔ | 2 notes ✖


