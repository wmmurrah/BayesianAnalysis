# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/DataAnalysis")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

# choose a data set:
# Influenza 
y <- c(14, 4, 5, 210)
#Hearing Loss 
#y <- c(20, 7, 103, 417)
# Rare Disease
#y <- c(0, 0, 13, 157)

n <- sum(y) # number of people/units measured

data <- list("y", "n") # to be passed on to WinBUGS
myinits <-	list(
  list(alpha = 0.5, beta = 0.5, gamma = 0.5))

# parameters to be monitored:	
parameters <- c("kappa","xi","psi","alpha","beta","gamma","pi")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="Kappa.txt",
	 			n.chains=1, n.iter=2000, n.burnin=1, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir,
	 			codaPkg=F, debug=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

plot(samples)

samples$mean$kappa

# Compare to Cohen's point estimate
p0 <- (y[1]+y[4])/n
pe <- (((y[1]+y[2]) * (y[1]+y[3])) + ((y[2]+y[4]) * (y[3]+y[4]))) / n^2
kappa.Cohen <- (p0-pe) / (1-pe) 
kappa.Cohen
