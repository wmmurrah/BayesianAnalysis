# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Binomial")

library(R2jags)

k1 <- 0
k2 <- 10
n1 <- 10
n2 <- 10

data <- list("k1", "k2", "n1", "n2") # to be passed on to JAGS
myinits <-	list(
  list(theta = 0.5))

# parameters to be monitored:	
parameters <- c("theta", "postpredk1", "postpredk2")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file ="Rate_5.txt", n.chains=1, n.iter=1000, 
         n.burnin=1, n.thin=1, DIC=T)

theta      <- samples$BUGSoutput$sims.list$theta 
postpredk1 <- samples$BUGSoutput$sims.list$postpredk1
postpredk2 <- samples$BUGSoutput$sims.list$postpredk2
 	 			 	 			
# Two-panel plot. 
layout(matrix(c(1,2),1,2))
layout.show(2)
# First, a histogram for theta.
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), ylim=c(0,10), xlab="Theta", ylab="Density") 
# let's plot a density estimate over this:
lines(density(theta), col="red", lwd=2)

# Second plot, the data space (predictives)
plot(k1,k2,type="p", pch=4, cex=2, lwd=2, xlab="Success Count 1", ylab="Success Count 2",
     xlim=c(-1, n1+1), ylim=c(-1,n2+1))        
nsamples <- length(theta)
sc <- 10
for (i in 0:n1)
{
  for (j in 0:n2) 
  {
    match.preds <- sum(postpredk1==i & postpredk2==j)/nsamples
    if (match.preds > 0)
    {
      points(i,j, pch=0, cex=sc*sqrt(match.preds)) 
    }
  }
}

