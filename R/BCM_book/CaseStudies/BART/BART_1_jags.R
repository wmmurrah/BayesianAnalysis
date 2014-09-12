# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/BART")
library(R2jags)

p       <- .15	# (Belief of) bursting probability
ntrials <- 90 # Number of trials for the BART

Data   <- matrix (data = as.numeric (as.matrix (read.table ("GeorgeSober.txt"))[-1,]), ntrials, 8)
d      <- matrix (, ntrials, 30) # Data in binary format

cash    <-(Data[,7]!=0)*1	# Cash or burst?
npumps  <- Data[,6]				# Nr. of pumps
options <- cash + npumps  # Nr. of decision possibilities

for (j in 1:ntrials)
{
	if (npumps[j]>0) {d[j, 1:npumps[j]] <- rep (0, npumps[j])}
	if (cash[j]==1) {d[j, (npumps[j]+1)] <- 1}
}

data <- list("ntrials", "p", "options", "d") # to be passed on to JAGS

myinits <-	list(
  list(gplus = 1.2, beta = 0.5))

# parameters to be monitored:
parameters <- c("gplus", "beta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters, model.file = "BART_1.txt", 
	         n.chains=1, n.iter=5000, n.burnin=2000, n.thin=1)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

gplus <- samples$BUGSoutput$sims.list$gplus
beta  <- samples$BUGSoutput$sims.list$beta

#################### PLOT RESULTS

par (cex.main = 2.5, cex.lab = 2, cex.axis = 1.5, mar = c(5, 5, 4, 0), las = 1)
layout (matrix (1:3, 1, 3))

hist(npumps, main = " ", xlab="Number of Pumps", ylab="Frequency", breaks=c(0:7), 
     col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,6.5,by=1), labels = c("1", "2", "3", "4", "5", "6", "7"))
     
plot(density(gplus), xlab = expression (gamma^'+'),
	main = " ", bty = 'n', lwd=2, ylab="Posterior Density")
plot(density(beta), xlab = expression (beta),
	main = " ", bty = 'n', lwd=2, lab = c(5,3,5), ylab="Posterior Density")

