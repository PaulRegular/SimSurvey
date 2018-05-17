
# Simple description of covariance:
# In 2D: X_ay = X_a,y-1 + X_a-1,y - X_a-1,y-1 + error
# at 1st age it is random walk in year,
# and first year it is random walk in age

library(SimSurvey)


sim_ay_covar <- function(sd = 0.5, phi_age = 0.7, phi_year = 0.9) {
  function(ages = NULL, years = NULL, cells = NULL, w = NULL) {

    na <- length(ages)
    ny <- length(years)
    nc <- length(cells)
    E <- array(rep(NA, nc * na * ny), dim = c(nc, na, ny),
               dimnames = list(cell = cells, age = ages, year = years))
    pc_age <- sqrt(1 - phi_age * phi_age)
    pc_year <- sqrt(1 - phi_year * phi_year)
    for (j in seq_along(years)) {
      for (i in seq_along(ages)) {
        if ((i == 1) & (j == 1)) {
          m <- 0
          s <- sd / pc_age
          E[, i, j] <- w %*% rnorm(nc, m, s)
        }
        if ((i > 1) & (j == 1)) {
          m <- phi_age * E[, i - 1, j]
          s <- sd / pc_year
          E[, i, j] <- w %*% rnorm(nc, m, s)
        }
        if ((i == 1) & (j > 1)) {
          m <- phi_year * E[, i, j - 1]
          s <- sd / pc_age
          E[, i, j] <- w %*% rnorm(nc, m, s)
        }
        if ((i > 1) & (j > 1)) {
          m <- phi_year * E[, i, j - 1] + phi_age * (E[, i - 1, j] - phi_year * E[, i - 1, j - 1])
          s <- sd
          E[, i, j] <- w %*% rnorm(nc, m, s)
        }
      }
    }
    E
  }
}


Zfun <- sim_arZ(mean = )
Zfun2 <- sim_Z()
years <- 1:20
ages <- 1:10
X <- Zfun(years = years, ages = ages)
Y <- Zfun2(years = years, ages = ages)

image(years, ages, t(X), col = viridis::viridis(100))
image(years, ages, t(Y), col = viridis::viridis(100))
#contour(years, ages, t(X), add = TRUE, col = "white")


mu <- matrix(rep(10, length(years) * length(ages)), nrow = length(ages), ncol = length(years))
mu[1,] <- 50
Zfun3 <- sim_arZ(mean = mu, phi_year = 0.9)
Z <- Zfun3(years = years, ages = ages)
image(years, ages, t(Z), col = viridis::viridis(100))















## agexyear correlated process error nll

sim_bpe <- bpe
sim_bpe[] <- NA
sim_bpe[, 1] <- 0
sim_bpe[1, ] <- 0
phi_age <- exp(logit_ar_pe_age) / (1 + exp(logit_ar_pe_age))
phi_year <- exp(logit_ar_pe_year) / (1 + exp(logit_ar_pe_year))
pc_age <- sqrt(1 - (phi_age * phi_age))
pc_year <- sqrt(1 - (phi_year * phi_year))

for(j in 2:Y){
  for(i in 2:A){
    if((i > 1) & (j == 1)){
      m <- phi_age * sim_bpe[i - 1, j]
      s <- std_pe / pc_year
    }
    if((i == 1) & (j > 1)){
      m <- phi_year * sim_bpe[i, j - 1]
      s <- std_pe / pc_age
    }
    if((i > 1) & (j > 1)){
      m <- phi_year * sim_bpe[i, j - 1] + phi_age *
        (sim_bpe[i - 1, j] - phi_year * sim_bpe[i - 1, j - 1])
      s <- std_pe
    }
    sim_bpe[i, j] <- rnorm(1, m, s)
  }
}
# image(1:Y, 1:A, t(bpe), col = colorRampPalette(c("white", "red"))(100),
#       xlab = "Year", ylab = "Age")
# image(1:Y, 1:A, t(sim_bpe), col = colorRampPalette(c("white", "red"))(100),
#       xlab = "Year", ylab = "Age")
# i <- 1
# plot(bpe[i, ], type = "b", ylim = range(c(bpe[i,], sim_bpe[i,])), ylab = i)
# points(sim_bpe[i, ], type = "b", col = "red")
# plot(year, bpe[1, ], ylim = range(bpe), pch = NA, ylab = "pe")
# for(i in 1:A) points(year, bpe[i, ], type = "l", col = i)
# plot(year, sim_bpe[1, ], ylim = range(sim_bpe), pch = NA, ylab = "pe")
# for(i in 1:A) points(year, sim_bpe[i, ], type = "l", col = i)

