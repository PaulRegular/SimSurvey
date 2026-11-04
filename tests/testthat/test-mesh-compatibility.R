test_that("legacy conversion preserves both bundled meshes and barrier indices", {
  skip_if_not_installed("fmesher")
  for (bundle in list(survey_mesh, survey_lite_mesh)) {
    original <- serialize(bundle, NULL)
    mesh <- .mesh_compat(bundle$mesh)
    expect_s3_class(mesh, "fm_mesh_2d")
    expect_identical(mesh$loc, bundle$mesh$loc)
    expect_identical(mesh$graph, bundle$mesh$graph)
    expect_identical(mesh$idx, bundle$mesh$idx)
    expect_identical(mesh$n, bundle$mesh$n)
    expect_identical(mesh$manifold, bundle$mesh$manifold)
    expect_identical(mesh$graph$tv[bundle$barrier_tri, , drop = FALSE],
                      bundle$mesh$graph$tv[bundle$barrier_tri, , drop = FALSE])
    expect_identical(.mesh_compat(mesh), mesh)
    expect_identical(serialize(bundle, NULL), original)
  }
})

test_that("mesh construction uses the current constructor without deprecations", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("fmesher")
  withr::local_seed(438)
  mesh <- make_mesh(small_grid(), max.edge = 10, bound.outer = 10, cutoff = 1)
  expect_s3_class(mesh, "fm_mesh_2d")
  expect_true(mesh$n > 0)
  expect_true(all(mesh$graph$tv >= 1 & mesh$graph$tv <= mesh$n))
})

test_that("both bundled legacy meshes work in standard SPDE simulation", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("fmesher")
  withr::local_seed(438)
  for (bundle in list(survey_mesh, survey_lite_mesh)) {
    original <- serialize(bundle, NULL)
    mesh <- .mesh_compat(bundle$mesh)
    locations <- bundle$mesh$loc[1:3, 1:2, drop = FALSE]
    projection <- INLA::inla.spde.make.A(mesh, locations)
    expect_equal(as.matrix(projection[, 1:3]), diag(3), ignore_attr = TRUE,
                 tolerance = 1e-10)
    expect_equal(as.numeric(Matrix::rowSums(projection)), rep(1, 3), tolerance = 1e-10)
    expect_equal(.Q(bundle$mesh, model = "spde"), .Q(mesh, model = "spde"))
    covariance <- sim_ays_covar_spde(mesh = bundle$mesh, model = "spde", sd = 0.1,
                                     group_ages = NULL)
    field <- covariance(locations, ages = 1:2, years = 2001:2002, cells = c(2, 4, 6))
    expect_equal(dim(field), c(2L, 2L, 3L))
    expect_equal(dimnames(field), list(age = c("1", "2"), year = c("2001", "2002"),
                                       cell = c("2", "4", "6")))
    expect_true(all(is.finite(field)))
    expect_identical(serialize(bundle, NULL), original)
  }
})

test_that("both bundled legacy meshes work with their existing barrier assignments", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("INLAspacetime")
  skip_if_not_installed("fmesher")
  withr::local_seed(438)
  for (bundle in list(survey_mesh, survey_lite_mesh)) {
    original <- serialize(bundle, NULL)
    mesh <- .mesh_compat(bundle$mesh)
    triangles <- bundle$barrier_tri
    expect_true(length(triangles) > 0)
    expect_true(all(triangles >= 1 & triangles <= nrow(mesh$graph$tv)))
    expect_equal(.Q(bundle$mesh, barrier.triangles = triangles, model = "barrier"),
                 .Q(mesh, barrier.triangles = triangles, model = "barrier"))
    covariance <- sim_ays_covar_spde(mesh = bundle$mesh, model = "barrier",
                                     barrier.triangles = triangles, sd = 0.1,
                                     group_ages = NULL)
    field <- covariance(bundle$mesh$loc[1:3, 1:2, drop = FALSE], ages = 1:2,
                         years = 2001:2002, cells = c(2, 4, 6))
    expect_equal(dim(field), c(2L, 2L, 3L))
    expect_true(all(is.finite(field)))
    expect_identical(serialize(bundle, NULL), original)
  }
})
