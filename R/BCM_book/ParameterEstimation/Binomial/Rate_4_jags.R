# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Binomial")

library(R2jags)

k <- 1
n <- 15
# Uncomment for Trompetter Data
# k <- 24
# n <- 121

data <- list("k", "n") # to be passed on to JAGS
myinits <-	list(
  list(theta = 0.5, thetaprior = 0.5))

# parameters to be monitored:	
parameters <- c("theta", "thetaprior", "postpredk", "priorpredk")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file ="Rate_4.txt", n.chains=1, n.iter=5000, 
         n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

######################Plots########################################################################
theta      <- samples$BUGSoutput$sims.list$theta
thetaprior <- samples$BUGSoutput$sims.list$thetaprior 
postpredk  <- samples$BUGSoutput$sims.list$postpredk 
priorpredk <- samples$BUGSoutput$sims.list$priorpredk

layout(matrix(c(1,2),2,1))
layout.show(2)
#Prior and posterior of theta
plot(density(theta, from=0, to=1), zero.line=F, axes=F, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0,6))
axis(1, at=c(0,0.2,0.4,0.6,0.8,1), lab=c("0","0.2","0.4","0.6","0.8","1"),cex.axis=0.8)
mtext("Rate", side=1, line=2.25, cex=1.2)
axis(2, at=c(0,2,4,6),cex.axis=0.8)
mtext("Density", side=2, line=2.25, cex=1.2)
lines(density(thetaprior, from=0, to=1), lty=3, col="gray")
legend(0.6,5.75, c("Prior", "Posterior"), lty=c(3,1), col=c ("grey", "black"))

#Prior and posterior predictive
mybreaks <- seq(from=-.5,to=n+1,by=1)
my.at    <- seq(from=0,to=n,by=1)
hist(postpredk,breaks=mybreaks,freq=F, right=F, ylab="", xlab="", ylim=c(0,0.3),main="", axes=F )
axis(1, at=my.at,lab=my.at,cex.axis=0.8)
mtext("Success Count", side=1, line=2.25, cex=1.2)
axis(2,at=c(0,0.1,0.2,0.3),lab=c("0","0.1","0.2","0.3"),cex.axis=0.8)
mtext("Mass", side=2, line=2.25, cex=1.2)
hist(priorpredk, breaks=mybreaks,freq=F,right=F,add=T, lty=3,border="grey")
legend(8,0.3, c("Prior", "Posterior"), lty=c(3,1),col=c("grey", "black"))

