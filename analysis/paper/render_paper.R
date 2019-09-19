
## Render both word and html documents
rmarkdown::render("analysis/paper/sim_survey_paper.Rmd", output_format = "all")

## Replace the +'s in the tables with empty space
ms <- readLines("analysis/paper/sim_survey_paper.html")
ms <- gsub("<code>[+]", "<code>&nbsp;", ms)
writeLines(ms, "analysis/paper/sim_survey_paper.html")
