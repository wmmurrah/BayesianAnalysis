# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Gaussian")

library(R2jags)

x <- c(-27.020,3.570,8.191,9.898,9.603,9.945,10.056)
n <- length(x)

data <- list("x", "n") # to be passed on to JAGS
myinits <- list(
  list(mu = 0, lambda = rep(1,n)))

# parameters to be monitored:	
parameters <- c("mu", "sigma")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file ="SevenScientists.txt", n.chains=1, n.iter=1000, 
         n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.


