#**************************************************************************
# Title: Sampling_grid_approximate_posterior.R
# Author: William Murrah
# Description:
# Created: Friday, 21 December 2018
# R version: R version 3.5.2 (2018-12-20)
# Directory: /home/wmmurrah/Projects/BayesianAnalysis
#**************************************************************************
# packages used -----------------------------------------------------------
library(rethinking)  
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


p_grid <- seq(0, 1, length.out = 1e3)
prior <- rep(1, 1e3)
prior <- ifelse(p_grid < .5, 0, 1)
likelihood <- dbinom(8, 15, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)

matplot(posterior, type = "l")

samples <- sample(p_grid, prob = posterior, 1e4, replace = TRUE)

plot(density(samples), xlim = c(0,1))
plot(samples, col = alpha("blue", .4))

dens(samples)
matplot(samples, type = "l", col = "blue")
which.max(samples)


# add up posterior probability where p < 0.5
sum(posterior[p_grid < 0.5])
sum(samples < 0.5)/1e4
quantile(samples, 0.8)

quantile(samples, c(0.1, 0.9))

PI(samples)
HPDI(samples)

w <- rbinom(1e4, size = 15, prob = samples)

simplehist(w)


plot(density(w/1e4))
barplot(table((w)))
simplehist(w)


samples2 <- dbinom(8, 15, prob = p_grid)        
table(samples2)

dens(samples2)
