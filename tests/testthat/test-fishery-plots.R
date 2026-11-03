test_that("population trends support biological matrices and singleton dimensions", {
  withr::local_seed(731)
  sim <- sim_abundance(ages = 3:4, years = 2001:2002)
  for (metric in c("N", "B", "S")) {
    p <- plotly::plotly_build(plot_trend(sim, sum_ages = 3, mat = metric))
    expect_equal(as.numeric(p$x$data[[1]]$y), as.numeric(sim[[metric]][1, ]))
    expect_match(p$x$data[[1]]$text[1], if (metric == "N") "fish" else "kg")
  }
  one <- sim_abundance(ages = 1, years = 2001)
  expect_equal(as.numeric(plotly::plotly_build(plot_trend(one))$x$data[[1]]$y), as.numeric(one$N))
  expect_s3_class(plot_trend(one, 1, "red"), "plotly")
  expect_error(plot_trend(sim, mat = "missing"), "mat")
  expect_error(plot_trend(sim, sum_ages = 99), "population ages")
})

test_that("fishery maps preserve cell totals, zeros and aggregate compositions", {
  withr::local_seed(732)
  result <- fishery_run(n_sims = 2, lengths_cap = 4, ages_cap = 1)
  saved <- data.table::copy(result)
  for (metric in c("operations", "n", "kg", "effort")) {
    p <- plotly::plotly_build(plot_fishery(result, "fishery", which_sim = NULL, map_by = metric))
    map <- p$x$data[[1]]
    expected <- if (metric == "operations") nrow(result$fishery$setdet) else sum(result$fishery$setdet[[metric]])
    expect_equal(sum(map$z), expected)
    expect_equal(sum(map$z > 0), 1)
    cell <- result$grid_xy[result$grid_xy$cell == 1, ]
    expect_equal(as.numeric(map$z[match(cell$y, map$y), match(cell$x, map$x)]), expected)
    expect_true(any(map$z == 0))
  }
  p <- plotly::plotly_build(plot_fishery(result, "fishery", which_sim = NULL))
  bars <- Filter(function(trace) identical(trace$type, "bar"), p$x$data)
  expect_equal(sum(bars[[1]]$y), 100)
  expect_equal(sum(bars[[2]]$y), 100)
  expect_equal(sum(bars[[3]]$y), 100)
  expect_equal(sum(bars[[4]]$y), 100)
  expect_equal(result, saved)
  expect_error(plot_fishery(result, "missing"), "one simulated fishery")
  expect_error(plot_fishery(result, "fishery", which_year = 1999), "which_year")
  expect_error(plot_fishery(result, "fishery", which_sim = 99), "which_sim")
})

test_that("annual fishery plots retain replicate values, medians, targets and units", {
  withr::local_seed(733)
  result <- fishery_run(fishery_population(years = 2001:2002), n_sims = 2,
                         gear_unit = "net", time_unit = "hour")
  p <- plotly::plotly_build(plot_fishery_trend(result, "fishery", metric = "n"))
  expect_equal(as.numeric(p$x$data[[1]]$y), result$fishery$annual$n)
  expect_equal(as.numeric(p$x$data[[2]]$y), c(50, 50))
  expect_equal(as.numeric(p$x$data[[3]]$y), c(50, 50))
  p <- plotly::plotly_build(plot_fishery_trend(result, "fishery", metric = "reported_biomass_cpue"))
  expect_equal(as.numeric(p$x$data[[1]]$y), result$fishery$annual$reported_biomass_cpue)
  expect_match(paste(unlist(p$x$layout$yaxis$title), collapse = " "), "kg / net \\* hour")
  expect_error(plot_fishery_trend(result, "fishery", "age"), "arg")
})

test_that("empty fisheries and huge aggregate catches remain plottable", {
  withr::local_seed(734)
  empty <- sim_fishery(fishery_population(), 0, capture_all)
  p <- plotly::plotly_build(plot_fishery(empty, "fishery"))
  expect_equal(sum(p$x$data[[1]]$z), 0)
  bars <- Filter(function(trace) identical(trace$type, "bar"), p$x$data)
  expect_true(all(vapply(bars, function(trace) sum(trace$y) == 0, logical(1))))
  expect_s3_class(plotly::plotly_build(plot_fishery_trend(empty, "fishery")), "plotly")
  huge <- fishery_run(fishery_population(1e10), lengths_cap = 2, ages_cap = 1)
  expect_lte(nrow(huge$fishery$samp), 2)
  p <- plotly::plotly_build(plot_fishery(huge, "fishery", map_by = "n"))
  expect_equal(sum(p$x$data[[1]]$z), sum(huge$fishery$annual$n))
  expect_lt(length(p$x$data[[1]]$z), 100)
})
