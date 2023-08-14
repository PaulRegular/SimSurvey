
load("data/survey_grid.rda")
load("data/survey_mesh.rda")

l <- lapply(c("cell", "strat", "division", "depth"), function(nm) {
  stars::st_as_stars(survey_grid[[nm]])
})
survey_grid <- c(l[[1]], l[[2]], l[[3]], l[[4]])

load("data/land.rda")
land <- sf::st_as_sf(land)

load("data/bathy.rda")
bathy <- stars::st_as_stars(bathy)

survey_lite_mesh$barrier_poly <- sf::st_as_sf(survey_lite_mesh$barrier_poly)
survey_mesh$barrier_poly <- sf::st_as_sf(survey_mesh$barrier_poly)

save(survey_grid, file = "data/survey_grid.rda", compress = "xz")
save(land, file = "data/land.rda", compress = "xz")
save(bathy, file = "data/bathy.rda", compress = "xz")
save(survey_mesh, survey_lite_mesh, file = "data/survey_mesh.rda", compress = "xz")
