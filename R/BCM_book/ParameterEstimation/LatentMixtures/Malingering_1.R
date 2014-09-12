rm(list=ls())
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/LatentMixtures")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

k <- c(45,45,44,45,44,45,45,45,45,45,30,20,6,44,44,27,25,17,14,27,35,30)
p <- length(k) # number of people
n <- 45        # number of questions

data    <- list("p", "k", "n") # to be passed on to WinBUGS
myinits <- list(
  list(psi = c(0.7,0.5), z = round(runif(p)))) # Initial group assignment

# parameters to be monitored:	
parameters <- c("psi","z")

samples <- bugs(data, inits=myinits, parameters,
               model.file="Malingering_1.txt",
               n.chains=1, n.iter=1000, n.burnin=1, n.thin=1,
               DIC=T, bugs.directory=bugsdir,
               codaPkg=F, debug=T)

