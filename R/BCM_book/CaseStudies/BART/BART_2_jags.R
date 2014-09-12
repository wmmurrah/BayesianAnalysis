# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/BART")
library(R2jags)

p       <- .15	# (Belief of) bursting probability
ntrials <- 90 # Number of trials for the BART

################### READ IN THE DATA

Data <- list(
	matrix (data = as.numeric (as.matrix (read.table ("GeorgeSober.txt"))[-1,]), ntrials, 8),
	matrix (data = as.numeric (as.matrix (read.table ("GeorgeTipsy.txt"))[-1,]), ntrials, 8),
	matrix (data = as.numeric (as.matrix (read.table ("GeorgeDrunk.txt"))[-1,]), ntrials, 8)
)

nconds <- length(Data)						             
cash   <- npumps <- matrix (, nconds, ntrials) # Cashes and nr. of pumps
d <- array (, c(nconds, ntrials, 30)) 			   # Data in binary format

for (i in 1:nconds)
{
	cash[i,]   <- (Data[[i]][,7]!=0)*1	# Cash or burst?
	npumps[i,] <- Data[[i]][,6]				  # Nr. of pumps

	for (j in 1:ntrials)
	{
		if (npumps[i,j]>0) {d[i, j, 1:npumps[i,j]] <- rep(0, npumps[i,j])}
		if (cash[i,j]==1) {d[i, j, (npumps[i,j]+1)] <- 1}
	}
}
options <- cash + npumps	# Nr. of decision possibilities

data <- list("nconds", "ntrials", "p", "options", "d") # to be passed on to JAGS
myinits <- list(
  list(mug = 1.2, sigmag = 0.1, mub = 0.8, sigmab = 0.8),
  list(mug = 1.5, sigmag = 0.2, mub = 1.0, sigmab = 1.2))

# parameters to be monitored:
parameters <- c("gplus", "beta", "mug", "sigmag", "mub", "sigmab")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters, model.file = "BART_2_jags.txt", 
	         n.chains=2, n.iter=5000, n.burnin=2000, n.thin=1)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

gplus.sober <- samples$BUGSoutput$sims.list$gplus[,1]
gplus.tipsy <- samples$BUGSoutput$sims.list$gplus[,2]
gplus.drunk <- samples$BUGSoutput$sims.list$gplus[,3]

beta.sober  <- samples$BUGSoutput$sims.list$beta[,1]
beta.tipsy  <- samples$BUGSoutput$sims.list$beta[,2]
beta.drunk  <- samples$BUGSoutput$sims.list$beta[,3]

#################### PLOT SOME RESULTS
par (cex.main = 2.5, cex.lab = 2, cex.axis = 1.5, mar = c(5, 5, 4, 0), las = 1)
layout (matrix (1:9, 3, 3, byrow = T))

hist(npumps[1,], xlab=" ", main = "Sober: # pumps", breaks=c(0:max(npumps[1,])), xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))

plot(density(gplus.sober), xlab = expression (gamma^'+'),
	main = expression (paste ("Sober: Posterior ", gamma^'+')), xlim = c(0.6,1.8), bty = 'n')
plot (density(beta.sober), xlab = expression (beta),
	main = expression (paste ("Sober: Posterior ", beta)), xlim = c(0.2,1.4), bty = 'n')

hist(npumps[2,], xlab=" ", main = "Tipsy: # pumps", breaks=c(0:max(npumps[2,])), xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))
plot(density (gplus.tipsy), xlab = expression (gamma^'+'),
	main = expression (paste ("Tipsy: Posterior ", gamma^'+')), xlim = c(0.6,1.8), bty = 'n')
plot(density (beta.tipsy), xlab = expression (beta),
	main = expression (paste ("Tipsy: Posterior ", beta)), xlim = c(0.2,1.4), bty = 'n')

hist(npumps[3,], xlab="Number of Pumps", main = "Drunk: # pumps", breaks=c(0:max(npumps[3,])), xlim = range(c(0:9)), col="lightblue", axes=F)
axis(2)
axis(1, at=seq(.5,8.5,by=1), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))
plot(density(gplus.drunk), xlab = expression (gamma^'+'),
	main = expression(paste ("Drunk: Posterior ", gamma^'+')), xlim = c(0.6,1.8), bty = 'n')
plot(density(beta.drunk), xlab = expression (beta),
	main = expression(paste ("Drunk: Posterior ", beta)), xlim = c(0.2,1.4), bty = 'n')

