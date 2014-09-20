#**************************************************************************
# Rate_1_runjags.R --------------------------------------------------------
# William Murrah
# Description: This script is a modification of Rate_1_jags.R to use the 
#              runjags package instead of the R2jags package.
# License: GPLv3
#**************************************************************************

# clears workspace:  
rm(list=ls()) 

# sets working directories:
# setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/GettingStarted")

library(runjags)

# Data
k <- 5
n <- 10

model <- '
# Inferring a Rate
model{
   # Prior Distribution for Rate Theta
   theta ~ dbeta(1,1)
   # Observed Counts
   k ~ dbin(theta,n)
}
'
data <- list(k=k, n=n) # to be passed on to JAGS

myinits <-	list(
  list(theta = 0.1), #chain 1 starting value
  list(theta = 0.9)) #chain 2 starting value

# parameters to be monitored:	
parameters <- c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the runjags documentation.
samples <- run.jags(model = model, monitor='theta', data=data, 
                    inits=myinits, 
                    n.chains=2, sample=20000, 
                    burnin=1, thin=1, plots=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# The commands below are useful for a quick overview:
samples         # a rough summary
plot(samples)   # a visual representation
plot(samples, type='trace') # traceplot 
plot(samples, type='all')
#more info on what is returned:
summary(samples)

# Use some plots in coda
# first use as.mcmmc to convert runjags object into mcmc.list:
samples.mcmc <- as.mcmc(samples)
# then use the plotting methods from coda:
xyplot(samples.mcmc)
densityplot(samples.mcmc)

# Use mcmc object to plot histogram
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(samples.mcmc, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), xlab="Rate", ylab="Posterior Density") 
