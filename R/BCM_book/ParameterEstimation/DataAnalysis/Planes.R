# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/DataAnalysis")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

x <- 10 # number of captures
k <- 4  # number of recaptures from n
n <- 5  # size of second sample
tmax <- 50 # maximum population size

data <- list("x","k","n","tmax") # to be passed on to WinBUGS

myinits <- list(
  list(t = 25))

parameters <- c("t")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="Planes.txt",
	 			n.chains=1, n.iter=4100, n.burnin=1, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir, #DIC=T or else error message
	 			codaPkg=F, debug=T)
# Now the values for theta are in the "samples" object, ready for inspection.

# Collect all samples in "t":
t <- samples$sims.list$t 

# Plot the posterior for theta:
windows(width=9,height=6) #Works only under Windows
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
hist(t, xlim=c(x+n-k,tmax), lty=1, lwd=2, col="grey", prob=T, breaks=(((x+n-k-1):tmax)+.5),
     axes=T, main=" ", xlab="Number of Planes", ylab="Posterior Mass")

