# JAGS may complain but the output looks good
# clears workspace:  
rm(list=ls()) 

setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/LatentMixtures")
library(R2jags)

k <- c(45,45,44,45,44,45,45,45,45,45,30,20,6,44,44,27,25,17,14,27,35,30)
p <- length(k) # number of people
n <- 45        # number of questions

data <- list("p", "k", "n") # to be passed on to JAGS
myinits <- list(
  list(z = round(runif(p)), mudiff=0.2),
  list(z = round(runif(p)), mudiff=0.3),
  list(z = round(runif(p)), mudiff=0.4)) 

# parameters to be monitored:	
parameters <- c("theta","z","mubon","lambdabon",
               "mumal","lambdamal","mudiff","phi")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
               model.file ="Malingering_2.txt", n.chains=3, n.iter=10000, 
               n.burnin=1000, n.thin=1, DIC=T)

plot(samples)	 			