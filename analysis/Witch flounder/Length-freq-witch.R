
#####################Length frequency for Witch flounder as we lack of age data############# 
#first run imitate_witch_fall_3NO

## Melt length frequency data
lf <- data.table::melt(setdet,
                       id.vars = c("survey.year", "vessel", "trip", "set", "easting", "northing", "set.depth.mean"),
                       measure.vars = names(setdet)[grepl("^lf", names(setdet))],
                       variable.name = "length", value.name = "freq")
lf <- lf[lf$length != "lfNA", ]
lf$length <- as.integer(gsub("lf", "", lf$length))

## Real by length
real_l <- data.frame(lf[lf$freq>0])
real_l  %>%
  ggplot(aes(x=set.depth.mean, y=freq, col=length))+
  geom_point() + scale_color_gradientn(colours = rainbow(5)) + theme_bw()

real_l  %>%
  ggplot(aes(x=set.depth.mean, y=freq)) +
  geom_point() + facet_wrap(~length)

## Real data (hold length or year and animate the other)
plot_ly(data = lf[lf$length == 30, ]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 5)
plot_ly(data = lf[lf$survey.year == 2003,]) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~length,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

# Examine at the age dimension, with frequency scaled by age to allow for
# distribution shifts at older ages to be visible
lf[lf$survey.year == 2011, ] %>%
  group_by(length) %>%
  mutate(scaled_freq = scale(freq)) %>%
  plot_ly() %>%
  add_markers(x = ~easting, y = ~northing, size = ~scaled_freq, frame = ~length,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

## Real data all length and years
plot_ly(data = lf) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 5)
plot_ly(data = lf) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~length,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

## Again, scale within length
lf %>%
  group_by(length) %>%
  mutate(scaled_freq = scale(freq)) %>%
  plot_ly() %>%
  add_markers(x = ~easting, y = ~northing, size = ~scaled_freq, frame = ~length,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)

###############Bubble plot: length at survey.year/set/set.depth.mean

ggplot(real_l, aes(x=survey.year, y=length, size = freq)) +
  geom_point(alpha=0.5)+
scale_size(range = c(0,7), name="Frequency")


