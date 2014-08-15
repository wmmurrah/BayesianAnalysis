#**************************************************************************
# Title:   simple linear regression ---------------------------------------
# Author:  William Murrah
# License: GPLv3
#**************************************************************************

# Data --------------------------------------------------------------------

wais <- read.table('data/wais.txt', header=TRUE)
wais <- wais[ ,-1]


# lm ----------------------------------------------------------------------

arith1.lm <- lm(arithmetic ~ picture.completion, wais)
summary(arith1.lm)

# Stan --------------------------------------------------------------------
library(rstan)

N <- 37
y <- wais$arithmetic
x <- wais$picture.completion

arith1.data <- c('N', 'y', 'x')

arith1.mod <- '
/*
 * ------------------------
 * Simple linear regression
 * ------------------------
 */

data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real b0;
  real b1;
  real<lower=0, upper=100> sigma;
}
model {
    b0 ~ normal(0, 100);
    b1 ~ normal(0, 100);
    y ~ normal(b0 + b1 * x, sigma);
}
'

arith1.fit <- stan(model_code = arith1.mod, data=arith1.data, 
                   iter = 10000, chains= 4)
arith1.fit
traceplot(arith1.fit)
plot(arith1.fit)

# runjags -----------------------------------------------------------------
library(runjags)
data.list <- list(y=y, x=x, N=N)
arith1.rjmod <- '
model {
  for (i in 1:N){ 
    mu[i] <- b0 + b1 * x[i] 
    y[i] ~ dnorm(mu[i], tau) 
  }
b0 ~ dunif(-1000, 1000)
b1 ~ dunif(-1000, 1000)
tau ~ dexp(1)
}
'
rj.fit <- run.jags(arith1.rjmod, data=data.list, n.chains = 4, 
                   monitor = c('b0', 'b1', 'tau'))
summary(rj.fit)



# MCMCpack ----------------------------------------------------------------
library(MCMCpack)

mc.fit <- MCMCregress(arithmetic ~ picture.completion, data=wais)
summary(mc.fit)
