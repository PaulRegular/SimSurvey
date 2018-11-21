
## wrapper for orca that puts file somewhere other than the working dir
export_plot <- function(p, file = NULL, width = 700, height = 500, scale = 4, ...) {
  ext <- tools::file_ext(file)
  orca(p, paste0("tmp.", ext), width = width, height = height, scale = scale, ...)
  invisible(file.rename(paste0("tmp.", ext), file))
}

## wrapper of add_annotations with defaults set-up for plot labeling
add_labels <- function(p, text = NULL, x = 0, y = 1, xref = "paper", yref = "paper",
                       font = list(size = 13), showarrow = FALSE, ..., data = NULL) {
  add_annotations(p, text = text, x = x, y = y, xref = xref, yref = yref,
                  font = font, showarrow = showarrow, ..., data = data)
}
