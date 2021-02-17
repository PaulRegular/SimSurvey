####Age frequency within Engel data series, witch flounder in Div. 3NO fall######################

setwd("~/Github/SimSurvey")
library(sf)
library(sp)
library(raster)
#library(SimSurvey)
library(dplyr)
library(ggplot2)
library(Rstrap)
library(plotly)
library(data.table)
library(tidyr)

## UTM projection for Newfoundland
utm_proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"

## Import strata shape file

strat_polys <- st_read("~/Github/SimSurvey/data-raw/DFO_NL_survey_strat/2HJ3KLNOP_Strata_Polygons_WGS84.shp",
                       layer = "2HJ3KLNOP_Strata_Polygons_WGS84")

## Subset to Division
strat_polys_div <- strat_polys %>%
  filter (DIV == "3N" | DIV == "3O")

strat_polys_div_utm <- strat_polys_div %>% st_transform(utm_proj)
###
con.setdet <- con.setdet[(con.setdet$rec == 5 |(con.setdet$rec == 6 & con.setdet$spec == 890)), ]
con.lf <- con.lf[con.lf$spec == 890, ]
ag <- ag[ag$spec == 890, ]
ag_fall<-subset(ag, ag$season=="fall")
sort(unique(ag_fall$age))    # range of witch age in the fall
ag_spring<-subset(ag, ag$season=="spring")
sort(unique(ag_spring$age))    # range of witch age in the spring
rv_data <- list(setdet = con.setdet, lf = con.lf, ag = ag)

out <- strat.fun(setdet = rv_data$setdet,lf = rv_data$lf, ag = rv_data$ag,
                 data.series = "Engel", program = "strat2 & strat1", which.survey = "multispecies",
                 species = 890, survey.year = c(1967:1993),
                 season = "fall",
                 NAFOdiv = c("3N", "3O"), strat = NULL,
                 sex = c("male","female","unsexed"),
                 length.group = 2, length.weight = NULL,
                 group.by = "length & age",
                 export = NULL, plot.results = FALSE)

## Convert lat and lon to UTM
setdet <- data.table(out$raw.data$set.details)
st_xy <- st_as_sf(data.frame(long = -setdet$long.start, lat = setdet$lat.start),
                  coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs(strat_polys_div_utm)) %>%
  st_coordinates() %>% data.frame(.)
names(st_xy) <- c("easting", "northing")
setdet <- cbind(setdet, st_xy)

## Melt age frequency data
af <- data.table::melt(setdet,
                       id.vars = c("survey.year", "vessel", "trip", "set", "easting", "northing", "set.depth.mean"),
                       measure.vars = names(setdet)[grepl("^af", names(setdet))],
                       variable.name = "age", value.name = "freq")
af <- af[af$age != "afNA", ]
af$age <- as.integer(gsub("af", "", af$age))

## Real by age
real_a <- data.frame(af[af$freq>0])
real_a  %>%
  ggplot(aes(x=set.depth.mean, y=freq, col=age))+
  geom_point() + scale_color_gradientn(colours = rainbow(5)) + theme_bw()

real_a  %>%
  ggplot(aes(x=set.depth.mean, y=freq)) +
  geom_point() + facet_wrap(~age)

## size up the distribution

symbols(setdet$long.start, setdet$lat.start,
        circles = sqrt(setdet$number / pi),
        inches = 0.1, main = "size of distribution - real data",
        xlab = "x", ylab = "y")
## Real data all ages and years
plot_ly(data = af) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~survey.year,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 5)
plot_ly(data = af) %>%
  add_markers(x = ~easting, y = ~northing, size = ~freq, frame = ~age,
              sizes = c(5, 1000), showlegend = FALSE) %>%
  animation_opts(frame = 500)
