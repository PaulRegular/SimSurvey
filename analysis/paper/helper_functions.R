
## http://www.moeding.net/archives/32-Metric-prefixes-for-ggplot2-scales.html
format_si <- function(...) {
  # Format a vector of numeric values according
  # to the International System of Units.
  # http://en.wikipedia.org/wiki/SI_prefix
  #
  # Based on code by Ben Tupper
  # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
  # Args:
  #   ...: Args passed to format()
  #
  # Returns:
  #   A function to format a vector of strings using
  #   SI prefix notation
  #

  function(x) {
    limits <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12,
                1e-9,  1e-6,  1e-3,  1e0,   1e3,
                1e6,   1e9,   1e12,  1e15,  1e18,
                1e21,  1e24)
    prefix <- c("y",   "z",   "a",   "f",   "p",
                "n",   "µ",   "m",   " ",   "k",
                "M",   "G",   "T",   "P",   "E",
                "Z",   "Y")

    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)

    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i == 0, which(limits == 1e0), i)

    paste0(format(round(x/limits[i], 1),
                  trim = TRUE, scientific = FALSE, ...),
           prefix[i])
  }
}
si <- format_si()

## little helper function for the fan plot
fan_plot <- function(d, x = "Age", which_strat = "Age") {

  if (which_strat == "Total") {
    ints_by = c("year", "set_den")
  }
  if (which_strat == "Age") {
    ints_by = c("year", "age", "set_den", "lengths_cap")
  }
  if (which_strat == "Length") {
    ints_by = c("year", "length", "set_den", "lengths_cap")
  }

  ## quantiles for fan plot
  quants <- seq(90, 10, by = -10)
  ints <- lapply(quants, function(q) {
    d[, list(prob = paste0(q, "%"),
             lower = quantile(I_hat, prob = (1 - q / 100) / 2),
             upper = quantile(I_hat, prob = 1 - (1 - q / 100) / 2)),
      by = ints_by]
  })
  ints <- rbindlist(ints)
  ints$prob <- factor(ints$prob, levels = paste0(quants, "%"))

  true <- d[d$sim %in% 1:3, ]

  if (x == "Age") {
    p <- ggplot(data = ints, aes(x = age))
    xlabs <- c("2", "4", "6", "8", "10+")
    xbreaks <- seq(2, 10, 2)
  }
  if (x == "Length") {
    p <- ggplot(data = ints, aes(x = length))
    xbreaks <- xlabs <- seq(0, 100, 20)
    xlabs[xlabs == 100] <- "100+"
  }
  if (x == "Year") {
    p <- ggplot(data = ints, aes(x = year))
    xlabs <- waiver()
    xbreaks <- waiver()
  }

  p <-  p +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = prob, alpha = prob)) +
    geom_line(aes(y = I), data = true, color = "black", size = 0.6) +
    scale_fill_viridis(name = "Pr", discrete = TRUE, begin = 0.1) +
    scale_alpha_discrete(name = "Pr", range = c(0.5, 1)) +
    xlab(x) + ylab("Abundance index") +
    scale_y_continuous(labels = format_si()) +
    scale_x_continuous(expand = c(0.01, 0),
                       labels = xlabs,
                       breaks = xbreaks) +
    guides(color = guide_legend(keywidth = 0.4, keyheight = 0.25,
                                default.unit = "cm", ncol = 1))

  if (which_strat == "Total") {
    p <- p + facet_wrap(~ set_den, nrow = 1,
                        labeller = label_bquote(cols = D[sets] == .(set_den)))
  } else {
    p <- p + facet_grid(lengths_cap ~ set_den,
                        labeller = label_bquote(cols = D[sets] == .(set_den),
                                                rows = M[lengths] == .(lengths_cap)))
  }

  p + theme_bw() +
    theme(legend.justification = c(0, 1),
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "darkgrey", size = 0.4),
          axis.ticks = element_line(color = "darkgrey", size = 0.4))

}


## wrapper for orca that puts file somewhere other than the working dir
export_plot <- function(p, file = NULL, width = 700, height = 500, scale = 4, ...) {
  ext <- tools::file_ext(file)
  orca(p, paste0("tmp.", ext), width = width, height = height, scale = scale, ...)
  invisible(file.rename(paste0("tmp.", ext), file))
}
