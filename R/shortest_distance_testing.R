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
