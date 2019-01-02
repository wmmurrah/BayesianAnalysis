#**************************************************************************
# Title: Gaussian_Height.R
# Author: William Murrah
# Description:
# Created: Thursday, 27 December 2018
# R version: R version 3.5.2 (2018-12-20)
# Directory: /home/wmmurrah/Projects/BayesianAnalysis
#**************************************************************************
# packages used -----------------------------------------------------------
library(tidyverse)  
library(rethinking)

data("Howell1")
d <- Howell1

str(d)

d$height
d2 <- d[d$age >= 18, ]

simplehist(d$height[d$age >=18])


with(d, simplehist(height[age >= 18]))

mean(d$height)



curve(dnorm(x, 178, 20), from = 100, to = 250)
curve(dunif(x, 0, 50), -10, 60)
curve(dt(x, 1, 1), -10, 60)

mod <- lm(height ~ 1, d2)
summary(mod)


# Rstan -------------------------------------------------------------------
library(rstan)

dat <- list(N = nrow(d2), 
            height = d2$height)


ht0 <- "
data {
  int<lower=0> N;
  vector[N] height;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  height ~ normal(mu, sigma);
  mu ~ normal(178, 20);
  sigma ~ uniform(0, 50);
}"

fit <- stan(model_code = ht0, data = dat)

summary(fit)
fit


# Mplus -------------------------------------------------------------------
library(MplusAutomation)

mplus_h0 <- mplusObject(usevariables = "height", 
  ANALYSIS = "estimator = bayes;",
  MODEL = "[height](mu);",
  rdata = d2,
  MODELPRIORS = "mu~N(178,20);"
  )

mplusModeler(mplus_h0, dataout = "mplus/height/height0.dat", 
             modelout = "mplus/height/height0.inp", 
             writeData = 'always',
             hashfilename = FALSE)
