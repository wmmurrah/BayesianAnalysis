# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/LatentMixtures")
library(R2jags)

dset <- 1

if (dset==1)
{
k <- c(1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      0,1,1,0,0,1,0,0,
      0,1,1,0,0,1,1,0,
      1,0,0,1,1,0,0,1,
      0,0,0,1,1,0,0,1,
      0,1,0,0,0,1,1,0,
      0,1,1,1,0,1,1,0)
k <- matrix(k, nrow=8, byrow=T)
}

if (dset==2)
{
k <- c(1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      0,1,1,0,0,1,0,0,
      0,1,1,0,0,1,1,0,
      1,0,0,1,1,0,0,1,
      0,0,0,1,1,0,0,1,
      0,1,0,0,0,1,1,0,
      0,1,1,1,0,1,1,0,
      1,0,0,1,NA,NA,NA,NA,
      0,NA,NA,NA,NA,NA,NA,NA,
      NA,NA,NA,NA,NA,NA,NA,NA)
k <- matrix(k, nrow=11, byrow=T)
}

if (dset==3)
{
k <- c(1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      1,0,0,1,1,0,0,1,
      0,1,1,0,0,1,0,0,
      0,1,1,0,0,1,1,0,
      1,0,0,1,1,0,0,1,
      0,0,0,1,1,0,0,1,
      0,1,0,0,0,1,1,0,
      0,1,1,1,0,1,1,0,
      1,0,0,1,NA,NA,NA,NA,
      0,NA,NA,NA,NA,NA,NA,NA,
      NA,NA,NA,NA,NA,NA,NA,NA)
k <- matrix(k, nrow=21, byrow=T)
}

nx <- nrow(k)
nz <- ncol(k)

data <- list("nx","nz","k") # to be passed on to JAGS
inits <-	list(
  list(z = round(runif(nz)), x = round(runif(nx)), alpha=0.5, beta=0.5))

if (dset==1)
{
# parameters to be monitored:	
  parameters <- c("z", "x", "alpha", "beta")

  # The following command calls JAGS with specific options.
  # For a detailed description see the R2jags documentation.
  samples <- jags(data, inits, parameters,
    	 			 model.file ="TwoCountryQuiz.txt", n.chains=1, n.iter=2000, 
             n.burnin=1000, n.thin=1, DIC=T)
  # Now the values for the monitored parameters are in the "samples" object, 
  # ready for inspection.
}

if (dset==2)
{
# parameters to be monitored:	
  parameters <- c("z", "x", "alpha", "beta", "NA.LP1", "NA.LP2", "NA.LP3")

  # The following command calls JAGS with specific options.
  # For a detailed description the R2jags documentation.
  samples <- jags(data, inits, parameters,
    	 			 model.file ="TwoCountryQuiz_NA1.txt", n.chains=1, 
             n.iter=2000, n.burnin=1000, n.thin=1, DIC=T)
  # Now the values for the monitored parameters are in the "samples" object, 
  # ready for inspection.
}

if (dset==3)
{
# parameters to be monitored:	
  parameters <- c("z", "x", "alpha", "beta", "NA.LP1", "NA.LP2", "NA.LP3")

  # The following command calls JAGS with specific options.
  # For a detailed description the R2jags documentation.
  samples <- jags(data, inits, parameters,
    	 			 model.file ="TwoCountryQuiz_NA2.txt", n.chains=1, 
             n.iter=2000, n.burnin=1000, n.thin=1, DIC=T)
  # Now the values for the monitored parameters are in the "samples" object, 
  # ready for inspection.
}
 
plot(samples)

