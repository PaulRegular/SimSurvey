test_that("generated grids retain the declared projected coordinate system", {
  grid <- small_grid()
  expected <- sf::st_crs("+proj=utm +zone=21 +datum=WGS84 +units=km +no_defs")
  expect_false(is.na(sf::st_crs(grid)))
  expect_equal(sf::st_crs(grid), expected)
  expect_false(sf::st_is_longlat(grid))
  expect_equal(abs(stars::st_res(grid)), c(x = 5, y = 5))
})
