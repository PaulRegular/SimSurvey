
library(ggplot2)


sp_pop <- sim_distribution(pop = sim_abundance(),
                           grid = sim_grid(res = c(3.5, 3.5)),
                           space_covar = sim_sp_covar(range = 100, sd = 0.1),
                           ay_covar = sim_ay_covar(sd = 2,
                                                   phi_age = 0.5,
                                                   phi_year = 0.05))
sp_pop$N
d <- sp_pop$sp_N
head(d)

plot(d$depth, d$N)


sim_pop <- sim_distribution()
ggplot(d[d$age %in% 1:5 & d$year %in% 1:5, ], aes(x = x, y = y)) +
  geom_tile(aes(fill = N)) + facet_grid(as.numeric(age) ~ as.numeric(year)) +
  scale_fill_viridis_c() + theme_void()
