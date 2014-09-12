# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Binomial")

library(R2jags)

k1 <- 5
k2 <- 7
n1 <- 10
n2 <- 10

data <- list("k1", "k2", "n1", "n2") # to be passed on to JAGS
myinits <-	list(
  list(theta = 0.5))

# parameters to be monitored:	
parameters <- c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file ="Rate_3.txt", n.chains=1, n.iter=1000, 
         n.burnin=1, n.thin=1, DIC=T)

theta <- samples$BUGSoutput$sims.list$theta
 	 			
# Now let's plot a histogram for theta. 
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(theta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(0,1), ylim=c(0,10), xlab="Rate", ylab="Posterior Density") 

