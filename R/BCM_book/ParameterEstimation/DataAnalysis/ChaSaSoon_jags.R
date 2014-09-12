#### After 949 failed attempts, Cha Sa-Soon from Korea finally passed   
#### her theoretical drivers exam. What can we infer about theta, the 
#### probability of answering any one question correctly?
     
# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/DataAnalysis")

library(R2jags)

nattempts  <- 950   
nfails     <- 949   
n          <- 50    # Number of questions
z = c(rep(NA,nfails),30)
y = c(rep(1,nfails),0)

data <- list("nattempts","n","z","y") # to be passed on to JAGS

myinits <- list(
  list(theta = 0.5))

parameters <- c("theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file ="ChaSaSoon_jags.txt", n.chains=1, n.iter=4100, 
         n.burnin=100, n.thin=1, DIC=T)
# Now the values for theta are in the "samples" object, ready for inspection.

# Collect all samples in "theta":
theta <- samples$BUGSoutput$sims.list$theta 

# Plot the posterior for theta:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(density(theta), ylim=c(0,200), xlim=c(0.2,0.5), lty=1, lwd=2, 
     axes=F, main=" ", xlab=" ", ylab="Posterior Density")
axis(1, at = c(0.20, 0.30, 0.40, 0.50), lab=c("0.20", "0.30", "0.40", "0.50"))
axis(2)
mtext(expression(paste(theta," for Cha Sa-soon")), side=1, line = 2.8, cex=2)

# plot 95% confidence interval
x0 <- quantile(theta,p=c(.025,.975))[[1]]
x1 <- quantile(theta,p=c(.025,.975))[[2]]
arrows(x0, 150, x1, 150, length = 0.05, angle = 90, code = 3, lwd=2)
text((x0+x1)/2, 170, labels="95%", cex = 1.5) 
text(x0-.03, 150, labels=as.character(round(x0,2)), cex = 1.5) 
text(x0+.04, 150, labels=as.character(round(x1,2)), cex = 1.5) 

