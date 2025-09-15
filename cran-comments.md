## R CMD check results

All checks returned **no ERRORs, WARNINGs, or NOTEs**.

### ✅ Local (Windows)

- Platform: Windows 11 x64 (build 22631), x86_64-w64-mingw32
- R version: 4.4.2 (2024-10-31 ucrt)
- Status: OK

### ✅ Win-builder (R-devel)

- Platform: Windows Server 2022, x86_64-w64-mingw32
- R version: R Under development (unstable) (2025-07-31 r88477)
- Status: OK

### ✅ GitHub Actions

- **Linux (R-devel)**  
  - Ubuntu 24.04.2 LTS, R Under development (unstable) (2025-07-31 r88477)  
  - Status: OK

- **Linux (R 4.5.1)**  
  - Ubuntu 24.04.2 LTS, R 4.5.1 (2025-06-13)  
  - Status: OK

- **macOS (R 4.5.1)**  
  - macOS Sonoma 14.7.6, aarch64-apple-darwin20  
  - Status: OK

- **Windows (R 4.5.1)**  
  - Windows Server 2022, R 4.5.1 (2025-06-13 ucrt)  
  - Status: OK

Some runs returned `INFO` notes due to Suggests or Enhances not in mainstream repositories: INLA

---

## Test environments

- Local: R 4.4.2 on Windows 11
- GitHub Actions: Linux (devel, release), macOS (release), Windows (release)
- Win-builder: Windows (devel)

## Changes in this version (0.1.7.1)

- Add `Matrix` to `Imports` so lazy-loaded mesh data (`survey_mesh`, `survey_lite_mesh`) pass checks. 
  No user-visible changes; API and results unchanged.


## Reverse dependencies

There are **no reverse dependencies** on CRAN.

