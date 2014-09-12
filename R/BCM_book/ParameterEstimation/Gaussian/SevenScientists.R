# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Gaussian")

library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

x <- c(-27.020,3.570,8.191,9.898,9.603,9.945,10.056)
n <- length(x)

data <- list("x", "n") # to be passed on to WinBUGS
myinits <- list(
  list(mu = 0, lambda = rep(1,n)))

# parameters to be monitored:	
parameters <- c("mu", "sigma")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="SevenScientists.txt",
	 			n.chains=1, n.iter=1000, n.burnin=1, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir,
	 			codaPkg=F, debug=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.


