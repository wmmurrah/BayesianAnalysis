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
likelihood <- dbinom(6, 9, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)

matplot(posterior, type = "l")

samples <- sample(p_grid, prob = posterior, 1e4, replace = TRUE)

plot(density(samples))
plot(samples, col = alpha("blue", .4))

dens(samples)
matplot(samples, type = "l", col = "blue")
which.max(samples)


# add up posterior probability where p < 0.5
sum(posterior[p_grid < 0.5])
sum(samples < 0.5)/1e4
quantile(samples, 0.5)
