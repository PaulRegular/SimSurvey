library(FSA)
library(FSAdata)
library(nlstools)

##Test
crm <- data(Croaker2)


typ <- vbStarts(tl~age,data=crm,type="original", plot = TRUE)
unlist(typ)



dat <- out$raw.data$age.growth %>%
  select(survey.year, strat, sex, length, age)

hist(dat$length, xlab = "length", main = "length frequency", breaks = 100)
hist(dat$age, xlab = "age", main = "age frequency", breaks = 30)

plot_ly() %>%
  add_markers(x = dat$age, y = dat$length) %>%
  layout(title = "Age Growth",
         xaxis = list(title = "Age"),
         yaxis = list(title = "Length"))

svOriginal <- vbStarts(length~age,data=dat,type="original", methLinf="oldAge", plot = TRUE)
unlist(svOriginal)
