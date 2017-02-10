plot_sim <- function(gird, zcol = NULL,
                     ramp = colorRampPalette(c("#FFFFFF", "#2166AC"))) {
  op <- par()
  layout(t(1:2), widths = c(0.9, 0.1))
  b <- seq(min(grid[[zcol]], na.rm = TRUE),
           max(grid[[zcol]], na.rm = TRUE),
           length.out = 225)
  col_ramp <- ramp(length(b))
  cols <- col_ramp[findInterval(grid[[zcol]], b)]
  par(mar = c(1, 1, 3, 0))
  plot(grid, col = cols, border = NA)
  box()
  par(mar = c(1, 0, 3, 3))
  image(1, b, t(seq_along(b)), col = col_ramp, axes = FALSE,
        main = zcol)
  axis(4, las = 2)
  box()
  on.exit(par(mar = op$mar, mfcol = c(1, 1)))
}
