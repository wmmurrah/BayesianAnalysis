# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Binomial")

library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

k1 <- 5
k2 <- 7
n1 <- 10
n2 <- 10

data <- list("k1", "k2", "n1", "n2") # to be passed on to WinBUGS
myinits <-	list(
  list(theta1 = 0.1, theta2 = 0.9))

# parameters to be monitored:	
parameters <- c("delta", "theta1", "theta2")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="Rate_2.txt",
	 			n.chains=1, n.iter=10000, n.burnin=1, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir,
	 			codaPkg=F, debug=F)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# Collect posterior samples:
delta <- samples$sims.list$delta 

# Now let's plot a histogram for delta. 
# First, some options to make the plot look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 80
y       <- hist(delta, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(-1,1), ylim=c(0,10), xlab="Difference in Rates", ylab="Posterior Density") 
 
# mean of delta:
mean(delta)
# median of delta:
median(delta)
# mode of delta, estimated from the "density" smoother:
density(delta)$x[which(density(delta)$y==max(density(delta)$y))]
# 95% credible interval for delta:
quantile(delta, c(.025,.975))