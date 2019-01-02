#**************************************************************************
# Title: GlobeExample.R
# Author: William Murrah
# Description:
# Created: Thursday, 27 December 2018
# R version: R version 3.5.2 (2018-12-20)
# Directory: /home/wmmurrah/Projects/BayesianAnalysis
#**************************************************************************
# packages used -----------------------------------------------------------
library(tidyverse)  
library(rethinking)

w <- 6; n <- 9;
p_grid <- seq(0, 1, length.out = 100)
posterior <- dbinom(w, n, prob = p_grid)*dunif(p_grid, 0, 1)
posterior <- posterior/sum(posterior)

