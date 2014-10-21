# R2jags Example;
# When you work through the code for the first time, 
# execute each command one at a time to better understand
# what it does.

# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Binomial")

library(runjags)

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

myinits <- list(
  list(theta = 0.1), #chain 1 starting value
  list(theta = 0.9)) #chain 2 starting value

# parameters to be monitored:	
parameters <- c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- run.jags(model=model, monitor='theta', data=data, inits=myinits, 
                     n.chains=2, sample=20000, 
                  burnin=1, thin=1)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# The commands below are useful for a quick overview:
print(samples)  # a rough summary
plot(samples)   # a visual representation
traceplot(samples) # traceplot (press <enter> repeatedly to see the chains)

#more info on what is returned:
summary(samples)
summary(samples$BUGSoutput)

chain <- 1
samples$BUGSoutput$sims.array[1:15,chain,]# array: element, chain, column (theta/deviance) 

# Collect posterior samples across all chains:
theta <- samples$BUGSoutput$sims.list$theta

# Now let's plot a histogram for theta. 
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), ylim=c(0,10), xlab="Rate", ylab="Posterior Density") 
# NB. ylim=c(0,10) defines the range of the y-axis. Adjust the upper value
# in case your posterior distribution falls partly outside this range.
