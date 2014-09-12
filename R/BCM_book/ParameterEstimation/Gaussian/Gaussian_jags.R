# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Gaussian")

library(R2jags)

x <- c(1.1, 1.9, 2.3, 1.8)
n <- length(x)

data <- list("x", "n") # to be passed on to JAGS
myinits <- list(
  list(mu = 0, sigma = 1))

# parameters to be monitored:	
parameters <- c("mu", "sigma")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file="Gaussian.txt", n.chains=1, n.iter=1000, 
         n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

mu    <- samples$BUGSoutput$sims.list$mu
sigma <- samples$BUGSoutput$sims.list$sigma 

