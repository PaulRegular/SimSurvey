## R CMD check results

All checks returned **no ERRORs, WARNINGs, or NOTEs**.

### ✅ Local (Windows)

- Platform: Windows 11 x64 (build 22631), x86_64-w64-mingw32
- R version 4.5.1 (2025-06-13 ucrt) 
- Status: OK

### ✅ Win-builder (R-devel)

- Platform: Windows Server 2022, x86_64-w64-mingw32
- R version: R Under development (unstable) (2025-09-14 r88831 ucrt)
- Status: OK

### ✅ GitHub Actions

- **Linux (R-devel)**  
  - Ubuntu 24.04.3 LTS, R-devel (2025-09-14 r88831)
  - Status: OK

- **Linux (R 4.5.1)**  
  - Ubuntu 24.04.3 LTS
  - Status: OK

- **macOS (R 4.5.1)**  
  - macOS Sequoia 15.6 (aarch64) 
  - Status: OK

- **Windows (R 4.5.1, ucrt)**  
  - Windows Server 2022  
  - Status: OK
  
- **Linux (R 4.4.3, oldrel-1)**
  - Ubuntu 24.04.3 LTS
  - Status: OK

---

## Test environments

- Local: R 4.5.1 on Windows 11
- GitHub Actions: Linux (devel, release), macOS (release), Windows (release)
- Win-builder: Windows (devel)

## Changes in this version (0.1.8)

- Add `Matrix` to `Imports` so lazy-loaded mesh data (`survey_mesh`, `survey_lite_mesh`) pass checks. 
  No user-visible changes; API and results unchanged.

## Reverse dependencies

There are **no reverse dependencies** on CRAN.

