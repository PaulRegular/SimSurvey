## Release summary

- This patch switches to dependence on the sf and stars package rather than the sp and raster package given the evolution of R spatial. We have also switched to new INLAspacetime functions for simulating from a barrier mesh.

## Test environments

- local R installation, R 4.3.0
- Github Actions ubuntu-latest (release, devel, oldrel-1), 
- win-builder (devel)
- R-hub windows-x86_64-oldrel (r-oldrel)
- R-hub windows-x86_64-devel (r-devel)

## R CMD check results

❯ On windows-x86_64-oldrel (r-oldrel)
  checking dependencies in R code ... NOTE
    All declared Imports should be used.

❯ On windows-x86_64-devel (r-devel)
  checking CRAN incoming feasibility ... [14s] NOTE

0 errors ✔ | 0 warnings ✔ | 2 notes ✖


