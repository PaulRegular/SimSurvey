#################################### Witch flounder #####################################################
####################################### 3NO Fall #################################################

##  Make grid function for each Division---------------
### Abundance and distribution--------------------
#### Simulated Data-----------------------------------------
##### Comparisons between real and simulated data  BETWEEN REAL: set catch, annual index,..----

#######################################################################################################

## Make grid function for each Division -----------------------------------------------------------------
setwd("~/Github/SimSurvey")
library(sf)
library(sp)
library(raster)
library(SimSurvey)
library(dplyr)
library(ggplot2)
library(plotly)

## UTM projection for Newfoundland
utm_proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

## Import strata shape file

strat_polys <- st_read("~/Github/SimSurvey/data-raw/DFO_NL_survey_strat/2HJ3KLNOP_Strata_Polygons_WGS84.shp",
                       layer = "2HJ3KLNOP_Strata_Polygons_WGS84")

## Subset to Division
strat_polys_div <- strat_polys %>%
                   filter (DIV == "3N" | DIV == "3O")

strat_polys_div_utm <- strat_polys_div %>% st_transform(utm_proj)

## Import bathy data
library(ncdf4)
strat_bathy <- raster("~/SimSurvey-doc/2HJ3KLNOP_GEBCO_bathy.nc")
strat_bathy <- raster::mask(strat_bathy, strat_polys_div)

## Strata-by-strata depths
strat_no <- unique(strat_polys_div$STRAT)
means <- numeric(length(strat_no))
mins <- numeric(length(strat_no))
maxs <- numeric(length(strat_no))
names(means) <- names(mins) <- names(maxs) <- strat_no

for (s in strat_no) {

   temp <- raster::mask(strat_bathy, strat_polys_div[strat_polys_div$STRAT == s, ])

   cat(s, "\n")
   cat("mean:", mean(temp[], na.rm = TRUE), "\n")
   cat("min:", min(temp[], na.rm = TRUE), "\n")
   cat("max:", max(temp[], na.rm = TRUE), "\n")

   means[as.character(s)] <- mean(temp[], na.rm = TRUE)
   mins[as.character(s)] <- min(temp[], na.rm = TRUE)
   maxs[as.character(s)] <- max(temp[], na.rm = TRUE)

 }

 depth_by_strata <- data.frame(strat = strat_no,
                               mean = means,
                               min = mins, max = maxs)

## Survey area
survey_area <- sum(sf::st_area(strat_polys_div_utm))
survey_area
# 135191.5 [km^2]

## Area of the strata
strat_sums <- strat_polys_div_utm %>%
  mutate(area = sf::st_area(strat_polys_div_utm)) %>%
  group_by(STRAT) %>%
  summarize(strat_area = sum(area))
range(strat_sums$strat_area)
# Units: [km^2]
# 207.7644 10359.0056

## Number of strata
strata <- length(unique(strat_polys_div_utm$STRAT))
strata
# 71

## Mean area
mean_area <- mean(strat_sums$strat_area)
mean_area
# 1904.106 [km^2]

## Range depth in the survey area
depth_range <- range(strat_bathy[], na.rm = TRUE)
depth_range
# -2102    -2

## Mean depth in the survey area
depth_mean <- mean(strat_bathy[], na.rm = TRUE)
depth_mean
# Units: [km^2]
# -206.2931

depth_median <- median(strat_bathy[], na.rm = TRUE)
depth_median
# Units: [km^2]
# -78

## Determines x_y_range
sqrt(survey_area)/2
# 183.842 [km]

## TO DO: Modify make_grid arguments to create a similar survey area
## shelf depth/width modifies mean grid depth, start breaks/splits modifies
## strata (length), method determines slope of shelf
library(bezier)
grid <- make_grid(x_range = c(-184, 184),
                  y_range = c(-184, 184),
                  res = c(3.5, 3.5),
                  shelf_depth = 60,
                  shelf_width = 170,
                  depth_range = c(0, 1600),
                  n_div = 2,
                  strat_breaks = seq(0, 1600, by = 65),
                  strat_splits = 4,
                  method ="bezier" )

tibble::lst(survey_area, strata, mean_area, depth_range, depth_mean, depth_median)

prod(res(grid)) * ncell(grid)
# 135056.2
length(unique(grid$strat))
# 72
mean(table(values(grid$strat)) * prod(res(grid)))
# 1875.781
range(values(grid$depth), na.rm = TRUE)
# 6 1415
mean(values(grid$depth), na.rm = TRUE)
# 208.0762
median(values(grid$depth), na.rm = TRUE)
# 83

xyz <- data.frame(rasterToPoints(grid$depth))
plot_ly(data = xyz, x = ~x, y = ~-depth) %>% add_lines()

plot(grid)
plot(grid$depth)
plot(grid$strat)
plot(rasterToPolygons(grid$strat, dissolve = TRUE), col = "grey")

grid_dat <- data.frame(raster::rasterToPoints(grid))  #11025 grid

strat_depths <- grid_dat %>%
  group_by(strat) %>%
  summarize(mean = mean(depth), min = min(depth), max = max(depth))

## Compare real vs simulated depths

real_depths <- data.frame(depth_by_strata)  #71 strata
real_depths <- abs(real_depths)
sim_depths <- data.frame(strat_depths)  #72 strata
all_depths <- rbind(real_depths, sim_depths)

all_depths <- all_depths %>% mutate(grid = factor(ifelse(strat > 300, "real", "sim")))

all_depths %>%
  filter(!is.na(grid)) %>%
  ggplot(aes(x=mean, color=grid, fill=grid)) +
  geom_histogram(position="dodge", bins=10, alpha = 1) + theme_bw()

plot_grid(grid)


### Abundance and distribution --------------------------------------------------

library(Rstrap)
library(plotly)
library(data.table)
library(tidyr)

### REAL DATA
# Subset data to Witch
con.setdet <- con.setdet[(con.setdet$rec == 5 |(con.setdet$rec == 6 & con.setdet$spec == 890)), ]
con.lf <- con.lf[con.lf$spec == 890, ]
ag <- ag[ag$spec == 890, ]
ag_fall<-subset(ag, ag$season=="fall")
sort(unique(ag_fall$age))    # range of witch age in the fall
rv_data <- list(setdet = con.setdet, lf = con.lf, ag = ag)

####### no investigation Rstrap output for both males and females

## Save Rstrap output by extracting data
out <- strat.fun(setdet = rv_data$setdet, lf = rv_data$lf, ag = rv_data$ag,
                 data.series = "Campelen", program = "strat2 & strat1", which.survey = "multispecies",
                 species = 890,
                 season = "fall",
                 NAFOdiv = c("3N", "3O"), strat = NULL,
                 sex = c("male","female","unsexed"),
                 length.group = 2, length.weight = NULL,
                 group.by = "length",        # No age-growth data
                 export = NULL, plot.results = FALSE)

## Convert lat and lon to UTM
setdet <- data.table(out$raw.data$set.details)
st_xy <- st_as_sf(data.frame(long = -setdet$long.start, lat = setdet$lat.start),
                  coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs(strat_polys_div_utm)) %>%
  st_coordinates() %>% data.frame(.)
names(st_xy) <- c("easting", "northing")
setdet <- cbind(setdet, st_xy)

## Set density
den <- setdet[rec == 5, list(n = .N), by = c("survey.year", "NAFOdiv", "strat", "strat.area")]
den <- den[, list(strat_area_nm = sum(strat.area),
                  strat_area_km = sum(strat.area * 3.4299),
                  den_km = sum(n) / sum(strat.area * 3.4299),
                  den_nm = sum(n) / sum(strat.area / 200)), by = c("survey.year")]
den

## ~ 0.001 sets per km^2
## ~ 1 sets per 200 sq. NM

### Melt length frequency data
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



#### Simulated Data ---------------------------------------------------------------------------

## Simulate data for comparison

# Abundance parameters and catchability curve roughly based on NCAM estimates
# Distribution parameters manually tweaked until results roughly corresponded
# to observations from 3NO witch
set.seed(890)
system.time(pop <- sim_abundance(ages = 1:16,
                                 years = 1:24,
                                 R = sim_R(log_mean = log(c(1500000,800000,2000000,4000000,
                                                            2000000,1000000,3000000,3000000,5000000,
                                                            5000000,5000000,5000000,18000000,22000000,
                                                            8500000,5000000,8000000,5000000,3500000,
                                                            2200000,3500000,3500000,6500000,10000000)),
                                           log_sd = 0.7,
                                           random_walk = T),
                                 Z = sim_Z(log_mean = log(0.08),
                                           log_sd = 0.75,
                                           phi_age = 0.8,
                                           phi_year = 0.9),
                                 N0 = sim_N0(N0 = "exp", plot = F),
                                 growth = sim_vonB(Linf = 60, L0 =5,  # max age and k from Bowering 1976
                                                   K = 0.2, log_sd = 0.1,
                                                   length_group = 2, digits = 0)) %>%
              sim_distribution(grid,
                               ays_covar = sim_ays_covar(sd = 1.5,
                                                         range = 800,
                                                         lambda = .55,
                                                         model = "matern",
                                                         phi_age = 0.6,
                                                         phi_year = 0.9,
                                                         #group_ages = 8:16
                                                         ),
                               depth_par = sim_parabola(mu = log(180),  #increase mu for catch and depth in simulated data
                                                        sigma = 0.8,
                                                        log_space = TRUE )))
## Quick look at distribution
sp_N <- data.frame(merge(pop$sp_N, pop$grid_xy, by = "cell"))    #sp_N, abundance (N) split by age,year and cell


for (j in rev(pop$ages)) {
  z <- xtabs(N ~ x + y, subset = age == j & year == 1, data = sp_N)
  image(z = z, axes = FALSE, col = viridis::viridis(100), main = paste("age", j))
}
for (i in rev(pop$years)) {
  z <- xtabs(N ~ x + y, subset = age == 1 & year == i, data = sp_N)
  image(z = z, axes = FALSE, col = viridis::viridis(100), main = paste("year", i))
}

###Sim SUrvey
system.time(survey <- sim_survey(pop,
                     n_sims = 1,
                     q = sim_logistic(k = 2, x0 = 2),
                     trawl_dim = c(1.5, 0.02),
                     resample_cells = FALSE,
                     binom_error = TRUE,
                     min_sets = 2,
                     set_den = 1/1000,
                     lengths_cap = 300,   # survey protocol for witch
                     ages_cap = 16,
                     age_sampling = "stratified",
                     age_length_group = 2,
                     age_space_group = "division",
                     light = FALSE)
)


##### Comparisons between real and simulated data----------------------------------
## Compare set catch

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
data_I<-out$strat1$length$abundance$annual.totals
data_I_l <- rowMeans(data_I[data_I$length %in% survey$lengths, grepl("y", names(data_I))])
barplot(data_I_l, names.arg = data_I$length , xlab = "length", main = "Index @ length - real data")
#sim data
sim_I_l <- rowMeans(survey$I_at_length)
barplot(sim_I_l, names.arg = names(sim_I_l), xlab = "length", main = " Index @ length- simulated data")

plot_ly() %>%
  add_lines(x = seq_along(data_I_l), y = data_I_l, name = "real") %>%
  add_lines(x = seq_along(sim_I_l), y = sim_I_l, name = "simulated") %>%
  layout(title = "Average index at length", xaxis = list(title = "Length"))
#
data_I <- out$strat1$length$abundance$details
sim_I <- survey$samp[aged == TRUE, ]
nrow(data_I)
nrow(sim_I)
hist(data_I$length, xlab = "length", main = "age growth data - real data", breaks = 100)
hist(sim_I$length, xlab = "length", main = "age growth data - simulated data", breaks = 100)

mean(data_I$length)
mean(sim_I$length)

plot_ly(alpha = 0.6, nbinsx = 30) %>%
  add_histogram(x = data_I$length, name = "real") %>%
  add_histogram(x = sim_I$length, name = "simulated") %>%
  layout(title = "Average Index at Length", xaxis = list(title = "Length"))

###### Compare relationship between catch and depth

data_I <- setdet
sim_I <- survey$setdet
plot(as.numeric(data_I$set.depth.min), data_I$number, xlab = "depth",
     ylab = "number", main = "real data", xlim = c(0, 2000))
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

