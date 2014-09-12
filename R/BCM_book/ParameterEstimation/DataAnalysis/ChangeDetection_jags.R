# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/DataAnalysis")

library(R2jags)

c <- scan("changepointdata.txt")
n <- length(c)
t <- 1:n

data <- list("c", "n", "t") # to be passed on to JAGS
myinits <- list(
  list(mu = c(1,1), lambda = 1, tau = n/2))

# parameters to be monitored:	
parameters <- c("mu","sigma","tau")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file ="ChangeDetection.txt", n.chains=1, n.iter=1000, 
         n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

plot(samples)
mean.tau <- samples$BUGSoutput$mean$tau
mean.mu1 <- samples$BUGSoutput$mean$mu[1]
mean.mu2 <- samples$BUGSoutput$mean$mu[2]

#some plotting options to make things look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
# the plotting: 
plot(c, type="l", main="", ylab="Values", xlab="Samples")
lines(c(1, mean.tau), c(mean.mu1,mean.mu1), lwd=2, col="red")
lines(c(mean.tau+1,length(c)), c(mean.mu2,mean.mu2), lwd=2, col="red")
