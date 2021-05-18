
# Compare Index at Age

dat <- out$strat1$length$abundance$annual.totals

#dat <- dat %>% select(survey.year, length, totals)

dat_age <- dat %>%
  mutate(age = case_when(
    between(length, 1, 7.9) ~ 1,
    between(length, 8, 13.9) ~ 2,
    between(length, 14, 18.9) ~ 3,
    between(length, 19, 23.9) ~ 4,
    between(length, 24, 27.9) ~ 5,
    between(length, 28, 30.9) ~ 6,
    between(length, 31, 33.9) ~ 7,
    between(length, 34, 36.9) ~ 8,
    between(length, 37, 39.9) ~ 9,
    between(length, 40, 41.9) ~ 10,
    between(length, 42, 43.9) ~ 11,
    between(length, 44, 44.9) ~ 12,
    between(length, 45, 46.9) ~ 13,
    between(length, 47, 47.9) ~ 14,
    between(length, 48, 48.9) ~ 15,
    between(length, 49, 49.9) ~ 16,
    between(length, 50, 50.9) ~ 17,
    between(length, 51, 51.9) ~ 18,
    between(length, 52, 52.9) ~ 19,
    between(length, 53, 100.9) ~ 20,
  ))

dat_age <- rowMeans(dat_age[dat_age$age %in% survey$ages, grepl("y", names(dat_age))])

mean(dat_age)

plot_ly() %>%
  add_lines(x = seq_along(data_I), y = data_I, name = "real") %>%
  add_lines(x = seq_along(sim_I), y = sim_I, name = "simulated") %>%
  add_lines(x = seq_along(dat_age), y = dat_age, name = "length") %>%
  layout(title = "Average index at age", xaxis = list(title = "Age"))

dat2 <- out$strat1$length$abundance$details

dat2 <- dat2 %>% select(survey.year, length, obs, totals)
dat2$abundance <- dat2$obs*dat2$totals

dat_age2 <- dat2 %>%
  mutate(age = case_when(
    between(length, 1, 7.9) ~ 1,
    between(length, 8, 13.9) ~ 2,
    between(length, 14, 18.9) ~ 3,
    between(length, 19, 23.9) ~ 4,
    between(length, 24, 27.9) ~ 5,
    between(length, 28, 30.9) ~ 6,
    between(length, 31, 33.9) ~ 7,
    between(length, 34, 36.9) ~ 8,
    between(length, 37, 39.9) ~ 9,
    between(length, 40, 41.9) ~ 10,
    between(length, 42, 43.9) ~ 11,
    between(length, 44, 44.9) ~ 12,
    between(length, 45, 46.9) ~ 13,
    between(length, 47, 47.9) ~ 14,
    between(length, 48, 48.9) ~ 15,
    between(length, 49, 49.9) ~ 16,
    between(length, 50, 50.9) ~ 17,
    between(length, 51, 51.9) ~ 18,
    between(length, 52, 52.9) ~ 19,
    between(length, 53, 100.9) ~ 20,
  ))

dat_age2 <- dat_age2 %>%
  group_by(age) %>%
  summarize(sum(abundance))

mean(dat_age2$`sum(abundance)`)


plot_ly() %>%
  add_lines(x = seq_along(data_I), y = data_I, name = "real") %>%
  add_lines(x = seq_along(sim_I), y = sim_I, name = "simulated") %>%
  add_lines(x = seq_along(dat_age), y = dat_age, name = "length") %>%
  add_lines(x = seq_along(dat_age2$age), y = dat_age2$`sum(abundance)`, name = "age") %>%
  layout(title = "Average index at age", xaxis = list(title = "Age"))


## Compare age growth data

data_I <- out$raw.data$age.growth
sim_I <- survey$samp[aged == TRUE, ]
nrow(data_I)
nrow(sim_I)
nrow(lf)
hist(data_I$length, xlab = "length", main = "age growth data - real data", breaks = 100)
hist(sim_I$length, xlab = "length", main = "age growth data - simulated data", breaks = 100)
hist(lf$length, xlab = "length", main = "age growth data - simulated data", breaks = 100)

mean(data_I$length)
mean(sim_I$length)
mean(lf$length)

plot_ly() %>%
  add_histogram(x = data_I$length, name = "real") %>%
  add_histogram(x = sim_I$length, name = "simulated") %>%
  add_histogram(x = lf$length, name = "new") %>%
  layout(title = "Age Growth")

lf <- lf %>%  filter(freq>0)
hist(lf$length, xlab = "length", main = "age growth data - simulated data", breaks = 100)

plot_ly() %>%
  add_markers(x = data_I$age - 0.25, y = data_I$length, name = "real") %>%
  add_markers(x = sim_I$age + 0.25, y = sim_I$length, name = "simulated") %>%
  layout(title = "Age Growth",
         xaxis = list(title = "Age"),
         yaxis = list(title = "Length"))

# Age Growth
al_dat <- out$raw.data$age.growth
setDT(al_dat)
al_dat <- al_dat %>% select(year, length, age, sex)

al_dat <- al_dat %>%
  group_by(age) %>%
  summarise(length = mean(length)) %>%
  as.data.table
