#**************************************************************************
# Title: bootstrappingHeightData.R
# Author: William Murrah
# Description:
# Created: Sunday, 13 January 2019
# R version: R version 3.5.2 (2018-12-20)
# Directory: /home/wmmurrah/Projects/BayesianAnalysis
#**************************************************************************
# packages used -----------------------------------------------------------
library(tidyverse)  
library(car)
library(rethinking)
library(cpar)


data("Howell1")
d2 <- Howell1[Howell1$age >= 18, ]

d2 <- d2[sample(nrow(d2), 20), ]
d2$weight.c <- d2$weight - mean(d2$weight)

mod4.3 <- lm(height ~ weight.c, d2)

mod4.3boot <- Boot(mod4.3)


bs_coefs <- mod4.3boot$t

plot(height ~ weight.c, data = d2, ylim = c(136, 180),
     xlim = c(-13, 19))
abline(a = coef(mod4.3)[1], b = coef(mod4.3)[2], col = "red")
for(i in seq_along(bs_coefs)) {
  abline(a = bs_coefs[i ,], b = bs_coefs[i ,2],
         col = col.alpha("black", 0.03), add = TRUE)
}


d2$bmi <- d2$weight/(d2$height/100)^2
#----------------------------------------

summary(mod4.3)
vcov(mod4.3)

sqrt(diag(vcov(mod4.3)))


mod4.3 %>% 
  vcov %>% 
  diag %>% 
  sqrt()



m4.3 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight.c,
    a ~ dnorm(156, 100),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = d2
)

precis(m4.3, corr = TRUE)