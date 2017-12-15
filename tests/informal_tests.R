
library(ggplot2)

sp_pop <- sim_distribution(pop = sim_abundance(),
                           grid = sim_grid(res = c(3.5, 3.5)),
                           space_covar = sim_sp_covar(range = 50, sd = 0.1),
                           ay_covar = sim_ay_covar(sd = 10,
                                                   phi_age = 0.5,
                                                   phi_year = 0.5),
                           depth_par = sim_parabola(alpha = 0, sigma = 50))
# sp_pop <- sim_distribution()
sp_pop$N
d <- sp_pop$sp_N
head(d)

p <- ggplot(d, aes(x = x, y = y)) +
  facet_grid(as.numeric(age) ~ as.numeric(year)) +
  scale_fill_viridis() + theme_void() +
  theme(panel.spacing = unit(-0.1, "lines"))
p + geom_tile(aes(fill = prob))
# p + geom_tile(aes(fill = N))

ggplot(d[sample(seq(nrow(d)), 10000), ]) + geom_point(aes(x = depth, y = N))
