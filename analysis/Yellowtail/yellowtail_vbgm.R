library(FSA)
library(FSAdata)
library(nlstools)

## Attempt to follow:
## http://derekogle.com/fishR/examples/oldFishRVignettes/VonBertalanffy.pdf
## http://derekogle.com/NCNRS349/modules/Growth/BKG

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

vbO <- vbFuns("Original")

svOriginal <- vbStarts(length~age,data=dat,type="Original", methLinf="oldAge", plot = TRUE)
unlist(svOriginal)
#       Linf           K          L0
# 50.00000000  0.01811917  7.74968669


nlsO <- nls(length~vbO(age,Linf,L0,K),data=dat,start=svOriginal)



#nlsO <- nls(length~vbO(age,Linf,L0,K),data=dat,start=svOriginal)
#cbind(Ests=coef(nlsO),confint(nlsO))
