
data_I <- setdet[setdet$number>0,]
hist(data_I$number, breaks = 100, xlab = "set catch", main = "set catch - real data")
sim_I <- survey$setdet[survey$setdet$n>0,]
hist(sim_I$n, breaks = 100, xlab = "set catch", main = "set catch - simulated data")

median(data_I$number)
median(sim_I$n)

plot_ly() %>%
  add_histogram(x = data_I$number, name = "real") %>%
  add_histogram(x = sim_I$n, name = "simulated") %>%
  layout(title = "Set catch")


### Compare annual index

data_I <- out$strat2$abundance$summary$total[survey$years]
names(data_I) <- survey$years
sim_I <- colSums(survey$I)
barplot(data_I, names.arg = names(data_I), xlab = "year", main = "annual index - real data")
barplot(sim_I, names.arg = names(sim_I), xlab = "year", main = "annual index - simulated data")

mean(data_I)
mean(sim_I)

plot_ly() %>%
  add_lines(data = out$strat2$abundance$summary,
            x = ~seq_along(survey.year), y = ~total, name = "real") %>%
  add_lines(x = seq(survey$years), y = colSums(survey$I), name = "simulated") %>%
  layout(title = "Annual index", xaxis = list(title = "Year"))


#### compare Index at length
##real_data

data_I <- out$strat1$length$abundance$details
sim_I <- survey$samp[aged == TRUE, ]
nrow(data_I)
nrow(sim_I)
hist(data_I$length, xlab = "length", main = "age growth data - real data", breaks = 100)
hist(sim_I$length, xlab = "length", main = "age growth data - simulated data", breaks = 100)

mean(data_I$length)
mean(sim_I$length)

plot_ly(alpha = 0.6, nbinsx = 100) %>%
  add_histogram(x = data_I$length, name = "real") %>%
  add_histogram(x = sim_I$length, name = "simulated") %>%
  layout(title = "Average Index at Length", xaxis = list(title = "Length"))


###### Compare relationship between catch and depth

data_I <- setdet
sim_I <- survey$setdet
plot(as.numeric(data_I$set.depth.min), data_I$number, xlab = "depth",
     ylab = "number", main = "real data", xlim = c(0, 1500))
plot(sim_I$depth, sim_I$n, xlab = "depth",
     ylab = "number", main = "simulated data", xlim = c(0,1500))

median(data_I$set.depth.mean)
median(sim_I$depth)

plot_ly() %>%
  add_markers(x = data_I$set.depth.mean, y = data_I$number, name = "real") %>%
  add_markers(x = sim_I$depth, sim_I$n, name = "simulated") %>%
  layout(title = "Compare Catch Depth", xaxis = list(title = "Depth"),
         yaxis = list(title = "Number"))

##### Compare intra-haul correlation
idvar <- c("vessel", "trip", "set", "year")
sub_lf <- merge(out$raw.data$set.details[, idvar], con.lf, by = idvar)
sub_lf$set <- as.numeric(as.factor(with(sub_lf, paste(vessel, trip, set, year))))
ind <- grep("^bin|^set$", names(con.lf)) # get length frequencies and expand to samples
lf <- sub_lf[, ind]
lf <- as.data.table(lf)
len_samp <- data.table::melt(lf, id.vars = "set")
len_samp <- len_samp[len_samp$value > 0, ]
len_samp$value <- round(len_samp$value)
len_samp$length <- as.numeric(gsub("bin", "", len_samp$variable))
len_samp <- as.data.frame(len_samp)
len_samp <- len_samp[rep(row.names(len_samp), len_samp$value), c("set", "length")]
len_samp <- data.table(len_samp)
sub_sets <- sample(len_samp$set, size = 10)
stripchart(length ~ set, data = len_samp[set %in% sub_sets, ],
           vertical = TRUE, pch = 1, cex = 0.5, method = "jitter", jitter = 0.2,
           main = "real data")
icc(len_samp$length, len_samp$set)

len_samp <- survey$samp[survey$samp$measured, ]
sub_sets <- sample(len_samp$set, size = 10)
stripchart(length ~ set, data = len_samp[set %in% sub_sets, ],
           vertical = TRUE, pch = 1, cex = 0.5, method = "jitter", jitter = 0.2,
           main = "simulated data")
icc(len_samp$length, len_samp$set)  #icc=how strongly measured length in the same set resemble each other

## Now size up the distribution

symbols(setdet$long.start, setdet$lat.start,
        circles = sqrt(setdet$number / pi),
        inches = 0.1, main = "size of distribution - real data",
        xlab = "x", ylab = "y")
symbols(survey$setdet$x, survey$setdet$y,
        circles = sqrt(survey$setdet$n / pi),
        inches = 0.1, main = "size of distribution - simulated data",
        xlab = "x", ylab = "y")

## Relationship of catch and depth by length
real_l <- data.frame(lf[lf$freq>0])
real_l  %>%
  ggplot(aes(x=set.depth.mean, y=freq, col=length))+
  geom_point() + scale_color_gradientn(colours = rainbow(5)) + theme_bw()

real_l  %>%
  ggplot(aes(x=set.depth.mean, y=freq)) +
  geom_point() + facet_wrap(~length)

# Examine the distributions of fish by length and year. Hold any year or any length constant to examine
# how they distibute. This will help determine phi_age and phi_year in sim_ays_covar

plot_ly(data = lf[lf$survey.year == 1999,]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~length,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

plot_ly(data = lf[lf$length == 40,]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)
