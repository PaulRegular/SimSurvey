library(SimSurvey)
library(e1071)
library(igraph)
library(raster)
data("survey_grid")


xy <- coordinates(survey_grid)
d <- as.matrix(dist(xy))
n <- 1000
d <- d[1:n, 1:n]

plot(x = xy[1:n, 1], y = xy[1:n, 2], cex = 5)
text(x = xy[1:n, 1], y = xy[1:n, 2], label = as.character(1:n))

g <- graph.adjacency(d, weighted=TRUE)
shortd1 <- shortest.paths(g, algorithm = "dijkstra")
shortd2 <- allShortestPaths(d)

d[1:10, 1:10]
shortd1[1:10, 1:10]
shortd2$length[1:10, 1:10]

d[c(1,83), c(1, 83)]
shortd1[c(1,83), c(1, 83)]
shortd2$length[c(1,83), c(1, 83)]

from <- 1
to <- 61
d[from, to]
shortd1[from, to]
shortd2$length[from, to]


extractPath(shortd2, from, to)



library(gdistance)

#example equivalent to that in the documentation on r.cost in GRASS
r <- raster(nrows=6, ncols=7, xmn=0, xmx=7, ymn=0, ymx=6, crs="+proj=utm +units=m")

r[] <- c(2, 2, 1, 1, 5, 5, 5,
         2, 2, 8, 8, 5, 2, 1,
         7, 1, 1, 8, 2, 2, 2,
         8, 7, 8, 8, 8, 8, 5,
         8, 8, 1, 1, 5, 3, 9,
         8, 1, 1, 2, 5, 3, 9)

T <- transition(r, function(x) 1/mean(x), 8)
# 1/mean: reciprocal to get permeability
T <- geoCorrection(T)

c1 <- c(5.5,1.5)
c2 <- c(1.5,5.5)

A <- accCost(T, c1)
plot(A)
text(A)





## build a graph with 5 nodes
x <- matrix(NA, 5, 5)
diag(x) <- 0
x[1,2] <- 30; x[1,3] <- 10
x[2,4] <- 70; x[2,5] <- 40
x[3,4] <- 50; x[3,5] <- 20
x[4,5] <- 60
x[5,4] <- 10
print(x)


## compute all path lengths
z <- allShortestPaths(x)
print(z)








g <- make_ring(10)
g <- make_lattice(c(10, 10))
g <- make_tree(20, mode = "undirected")
plot(g)
distances(g)
shortest_paths(g, 5)
all_shortest_paths(g, 1, 6:8)
mean_distance(g)



## Weighted shortest paths
el <- matrix(nc=3, byrow=TRUE,
             c(1,2,5, 1,3,5, 1,4,5, 2,3,5, 2,5,5, 2,6,5,
               3,7,5, 4,3,5, 4,7,5, 5,6,5, 5,8,5, 6,3,5, 6,7,5, 6,9,5,
               6,10,5, 8,6,5, 8,9,5, 9,10,5) )
g2 <- add_edges(make_empty_graph(10, directed = FALSE), t(el[,1:2]), weight=el[,3])
plot(g2)
distances(g2)



n <- 10
temp <- xy[1:n, ]
d <- dist(temp)
m <- data.frame(t(combn(rownames(as.matrix(d)),2)), as.numeric(d))
names(m) <- c("c1", "c2", "distance")


g2 <- add_edges(make_empty_graph(10), t(m[,1:2]), weight=m[,3])
plot(g2)



## Make a graph weighted by distances of adjacent cells

n <- 1000
grid <- survey_grid[1:n, ]
xy <- coordinates(grid)

adj_mat <- rgeos::gTouches(grid, byid = TRUE)
dist_mat <- d <- as.matrix(dist(xy))
dist_mat[!adj_mat] <- NA

#plot(grid)
plot(x = xy[1:n, 1], y = xy[1:n, 2], cex = 5)
text(x = xy[1:n, 1], y = xy[1:n, 2], label = as.character(1:n))

up_dist_mat <- dist_mat
up_dist_mat[lower.tri(up_dist_mat)] <- NA
m <- expand.grid(c1 = seq.int(nrow(dist_mat)), c2 = seq.int(nrow(dist_mat)))
m$dist <- as.numeric(up_dist_mat)
m <- na.omit(m)
g <- add_edges(make_empty_graph(n, directed = FALSE), t(m[, c("c1", "c2")]), weight = m$dist)
#plot(g)

system.time(shortd1 <- distances(g, algorithm = "dijkstra"))
system.time(shortd2 <- allShortestPaths(dist_mat))

all.equal(shortd1, shortd2$length)

dist_mat[1:10, 1:10]
d[1:10, 1:10]
shortd1[1:10, 1:10]
shortd2$length[1:10, 1:10]

from <- 1
to <- 61
d[from, to]
shortd2$length[from, to]
extractPath(shortd2, from, to)


## Conclusion ----------


library(SimSurvey)
library(igraph)
library(raster)
data("survey_grid")


n <- nrow(survey_grid)
grid <- survey_grid[1:n, ]
xy <- coordinates(grid)

adj_mat <- rgeos::gTouches(grid, byid = TRUE)
dist_mat <- d <- as.matrix(dist(xy))
dist_mat[!adj_mat] <- NA

#plot(grid)
#plot(x = xy[1:n, 1], y = xy[1:n, 2], cex = 5)
#text(x = xy[1:n, 1], y = xy[1:n, 2], label = as.character(1:n))

up_dist_mat <- dist_mat
up_dist_mat[lower.tri(up_dist_mat)] <- NA
m <- expand.grid(c1 = seq.int(nrow(dist_mat)), c2 = seq.int(nrow(dist_mat)))
m$dist <- as.numeric(up_dist_mat)
m <- m[!is.na(m$dist), ]
g <- add_edges(make_empty_graph(n, directed = FALSE), t(m[, c("c1", "c2")]), weight = m$dist)
#plot(g)

system.time(shortd1 <- distances(g, algorithm = "dijkstra"))

d <- shortd1
#d <- as.matrix(dist(xy))
W <- exp(- d / 21)
sigma <- chol(W)
grid$e <- sigma %*% rnorm(n)
plot_sim(grid, zcol = "e")

## not working like I pictured...



## IPDW testing ------------------

library("ipdw")

help(package = "ipdw")
vignette("ipdw", package = "ipdw")


library("ipdw")
library("geoR")
data(kattegat)
katproj<-c("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs")
pols1<-Polygons(list(Polygon(kattegat$dk[1])),"pol1")
pols2<-Polygons(list(Polygon(kattegat$dk[2])),"pol2")
pols3<-Polygons(list(Polygon(kattegat$dk[3])),"pol3")
pols4<-Polygons(list(Polygon(kattegat$dk[4])),"pol4")
pols5<-Polygons(list(Polygon(kattegat$dk[5])),"pol5")
pols6<-Polygons(list(Polygon(kattegat$dk[6])),"pol6")
pols7<-Polygons(list(Polygon(kattegat$dk[7])),"pol7")
pols8<-Polygons(list(Polygon(kattegat$dk[8])),"pol8")
pols9<-Polygons(list(Polygon(kattegat$dk[9])),"pol9")
pols10<-Polygons(list(Polygon(kattegat$dk[10])),"pol10")
pols11<-Polygons(list(Polygon(kattegat$dk[11])),"pol11")
pols12<-Polygons(list(Polygon(kattegat$dk[12])),"pol12")
pols<-SpatialPolygons(list(pols1,pols2,pols3,pols4,pols5,pols6,
                           pols7,pols8,pols9,pols10,pols11,pols12),1:12)


projection(pols)<-katproj
costras<-costrasterGen(kattegat$coords,pols,extent="pnts",katproj)
#insert contiguous barrier
costras[160:170,1:80] <- 10000



#find average nearest neighbor
library(spatstat)
W=owin(range(kattegat$coords[,1]),range(kattegat$coords[,2]))
kat.pp<-ppp(kattegat$coords[,1],kattegat$coords[,2],window=W)
mean.neighdist<-mean(nndist(kat.pp))
#grid building
gridsize<-mean.neighdist*2
grainscale.fac<-gridsize/res(costras)[1]
gridras<-aggregate(costras,fact=grainscale.fac)
gridpol<-rasterToPolygons(gridras)
gridpol$value<-row.names(gridpol)
#spatial join
kat.df<-data.frame(kattegat)
coordinates(kat.df)<-~x.utm+y.utm
projection(kat.df)<-katproj
fulldataset.over<-over(kat.df,gridpol)
fulldataset.over<-cbind(data.frame(fulldataset.over),data.frame(kat.df))
#grid selection
set.seed(2)
gridlev<-unique(fulldataset.over$value)
for(i in 1:length(gridlev)){
  activesub<-subset(fulldataset.over,fulldataset.over$value==gridlev[i])
  selectnum<-gdata::resample(1:nrow(activesub),1)
  if(i==1){
    training<-activesub[selectnum,]
  }
  else{
    training<-rbind(training,activesub[selectnum,])
  }
}



validate<-fulldataset.over[!(row.names(fulldataset.over) %in% row.names(training)),]
xy<-cbind(training$x.utm,training$y.utm)
training<-SpatialPointsDataFrame(xy,training)
xy<-cbind(validate$x.utm,validate$y.utm)
validate<-SpatialPointsDataFrame(xy,validate)
projection(training)<-katproj
projection(validate)<-katproj



paramlist <- c("data")

final.ipdw <- ipdw(training, costras, range = mean.neighdist * 10, paramlist)
















## gdistance testing --------------------------

library("ipdw")
library("geoR")
data(kattegat)
katproj<-c("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs")
pols1<-Polygons(list(Polygon(kattegat$dk[1])),"pol1")
pols2<-Polygons(list(Polygon(kattegat$dk[2])),"pol2")
pols3<-Polygons(list(Polygon(kattegat$dk[3])),"pol3")
pols4<-Polygons(list(Polygon(kattegat$dk[4])),"pol4")
pols5<-Polygons(list(Polygon(kattegat$dk[5])),"pol5")
pols6<-Polygons(list(Polygon(kattegat$dk[6])),"pol6")
pols7<-Polygons(list(Polygon(kattegat$dk[7])),"pol7")
pols8<-Polygons(list(Polygon(kattegat$dk[8])),"pol8")
pols9<-Polygons(list(Polygon(kattegat$dk[9])),"pol9")
pols10<-Polygons(list(Polygon(kattegat$dk[10])),"pol10")
pols11<-Polygons(list(Polygon(kattegat$dk[11])),"pol11")
pols12<-Polygons(list(Polygon(kattegat$dk[12])),"pol12")
pols<-SpatialPolygons(list(pols1,pols2,pols3,pols4,pols5,pols6,
                           pols7,pols8,pols9,pols10,pols11,pols12),1:12)



r <- raster(extent(pols), nrow = 50, ncol = 50)
r <- rasterize(pols, r)
r[!is.na(r)] <- Inf
r[is.na(r)] <- 1
r[3000:3300] <- Inf # make things interesting; add a road
plot(r)
xy <- rasterToPoints(r)
#xy <- data.frame(xy[sample(which(r[] == 1), 100), ])
xy <- data.frame(xy[which(r[] == 1), ])
xy <- SpatialPointsDataFrame(cbind(xy[1], xy[2]), xy[3])
proj4string(xy) <- proj4string(r)
plot(xy, add = TRUE)

trans <- gdistance::transition(r, function(x) 1/mean(x), directions = 16)
trans <- geoCorrection(trans)
trans[1:10, 1:10]

costs <- gdistance::accCost(trans, xy)
plot(costs) ## cool!
plot(pols, add = TRUE)
#text(coordinates(xy)[, 1], coordinates(xy)[, 2], seq_along(xy), cex = 0.5)

costd <- costDistance(trans, xy)
crowd <- dist(coordinates(xy))

as.matrix(costd)[1:10, 1:10]
as.matrix(crowd)[1:10, 1:10]

from <- 5
to <- 37
as.matrix(costd)[from, to]
as.matrix(crowd)[from, to]


d <- as.matrix(costd)
#d <- as.matrix(crowd)
n <- nrow(d)
W <- exp(- d / 20)
sigma <- chol(W)
xy$e <- t(sigma) %*% rnorm(n)
er <- rasterize(xy, r)$e
#er <- disaggregate(er, fact = 8, method = "bilinear")
plot(er)


focal <- sample(seq_along(xy), 1)
xy$focal <- W[focal,]
fr <- rasterize(xy, r)$focal
#fr <- disaggregate(fr, fact = 4, method = "bilinear")
plot(fr)

plot(d[focal, ], W[focal, ], xlim = c(0, 500))


## TO DO: find out why "the leading minor of order ___ is not positive definite"!?


# ## Matern covariance
#
# dmat <- as.matrix(costd)
# sigma2e <- 0.1; sigma2x <- 20; kappa <- 0.2; nu <- 1
#
# mcor <- as.matrix(2^(1-nu)*(kappa*dmat)^nu*
#                     besselK(dmat*kappa,nu)/gamma(nu))
# diag(mcor) <- 1; mcov <- sigma2e*diag(n) + sigma2x*mcor
#
# L <- chol(mcov)
# xy$e <- drop(rnorm(n) %*% L)
#
# er <- rasterize(xy, r)$e
# er <- disaggregate(er, fact = 8, method = "bilinear")
# plot(er)
#
#
# focal <- sample(seq_along(xy), 1)
# xy$focal <- mcor[focal,]
# fr <- rasterize(xy, r)$focal
# #fr <- disaggregate(fr, fact = 4, method = "bilinear")
# plot(fr)
# plot(dmat[focal, ], mcor[focal, ], xlim = c(0, 500))
#
# ## Looks like the Matern covariance function is not going to work using least cost distance



## Size up the Egan values of the exponential autocovarianc function

d <- as.matrix(costd)
h <- seq(10, 200, by = 10)
e <- rep(NA, length(h))
for(i in seq_along(h)) {
  W <- exp(- d / h[i])
  e[i] <- min(eigen(W)$values)
}
plot(h, e, type = "o", pch = 16)
abline(h = 0, lty = 3)

## least-cost distance matrix is not positive-difinitive...
## don't know how to make a positive-difinitive least-cost covariance matrix...yet


