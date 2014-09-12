# R2WinBUGS Example;
# When you work through the code for the first time, 
# execute each command one at a time to better understand
# what it does.

# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/Binomial")

library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

k <- 5
n <- 10

data <- list("k", "n") # to be passed on to WinBUGS

myinits <- list(
  list(theta = 0.1), #chain 1 starting value
  list(theta = 0.9)) #chain 2 starting value

# parameters to be monitored:	
parameters <- c("theta")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="Rate_1.txt",
	 			n.chains=2, n.iter=20000, n.burnin=1, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir,
	 			codaPkg=F, debug=F)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# The commands below are useful for a quick overview:
print(samples)  # a rough summary
plot(samples)   # a visual representation
names(samples)  # summarizes the variables
samples$summary # more detailed summary
chain <- 1
samples$sims.array[1:15,chain,]# array: element, chain, column (theta/deviance)

# Collect posterior samples across all chains:
theta <- samples$sims.list$theta 

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
