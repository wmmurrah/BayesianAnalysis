# In order to get JAGS to work:
# (1) Install the latest version of JAGS from sourceforge (e.g., 
# http://sourceforge.net/projects/mcmc-jags/files/)
# (2) Install the latest version of rjags. Do **not** try to install this via the
# usual route (i.e., R -> Packages -> Install package(s)) as this may install
# an old version of rjags that expects an old version of JAGS.
# Instead, you want to Google for rjags CRAN, go to a site 
# such as http://cran.r-project.org/web/packages/rjags/index.html, and --when
# using Windows-- download the .zip file. Then, in R, go to Packages ->
# Install package(s) from local zip file...
# To check, type library(rjags) at the R prompt.
# (3) Install R2jags. This does work via the usual route 
# R -> Packages -> Install package(s).
# http://yusung.blogspot.nl/2008/05/r2jags-package-for-running-jags-from-r.html
# Now you are all set.

# When you work through the code for the first time, 
# execute each command one at a time to better understand
# what it does.

# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/GettingStarted")

library(R2jags)

k <- 5
n <- 10

data <- list("k", "n") # to be passed on to JAGS

myinits <-	list(
  list(theta = 0.1), #chain 1 starting value
  list(theta = 0.9)) #chain 2 starting value

# parameters to be monitored:	
parameters <- c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			  model.file ="Rate_1.txt", n.chains=2, n.iter=20000, 
          n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# The commands below are useful for a quick overview:
print(samples)  # a rough summary
plot(samples)   # a visual representation
traceplot(samples) # traceplot (press <enter> repeatedly to see the chains)

#more info on what is returned:
summary(samples)
summary(samples$BUGSoutput)

samples$BUGSoutput$sims.array[1:15,,2]# array: sample, chain, parameter 

# Collect posterior samples across all chains:
theta <- samples$BUGSoutput$sims.list$theta
 
# Now let's plot a histogram for theta. 
# NB. Some the plots will not look good in RStudio.
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), xlab="Rate", ylab="Posterior Density") 

# Additional option: use some plots in coda
# first use as.mcmmc to convert rjags object into mcmc.list:
samples.mcmc <- as.mcmc(samples)
# then use the plotting methods from coda:
xyplot(samples.mcmc)
densityplot(samples.mcmc)
