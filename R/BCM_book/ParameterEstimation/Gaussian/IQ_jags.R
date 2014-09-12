# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Gaussian")

library(R2jags)

x <- matrix(c(90,95,100,105,110,115,150,155,160),nrow=3,ncol=3,byrow=T) 
x

n <- nrow(x) # number of people
m <- ncol(x) # number of repeated measurements

data <- list("x", "n", "m") # to be passed on to JAGS
myinits <- list(
  list(mu = rep(100,n), sigma = 1))

# parameters to be monitored:	
parameters <- c("mu", "sigma")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file ="IQ.txt", n.chains=1, n.iter=1000, 
         n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.


