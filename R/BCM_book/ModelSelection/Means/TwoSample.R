# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ModelSelection/Means")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

x <- c(70,80,79,83,77,75,84,78,75,75,78,82,74,81,72,70,75,72,76,77)
y <- c(56,80,63,62,67,71,68,76,79,67,76,74,67,70,62,65,72,72,69,71)

n1 <- length(x)
n2 <- length(y)

# Rescale
y <- y - mean(x)
y <- y/sd(x)
x <- (x-mean(x))/sd(x); 

data <- list("x", "y", "n1", "n2") # to be passed on to WinBUGS

myinits <- list(
  list(delta = rnorm(1,0,3), mu = rnorm(1,0,1), sigmatmp = runif(1,0,5)),
  list(delta = rnorm(1,0,3), mu = rnorm(1,0,1), sigmatmp = runif(1,0,5)),
  list(delta = rnorm(1,0,3), mu = rnorm(1,0,1), sigmatmp = runif(1,0,5)))

# Parameters to be monitored
parameters <- c("delta")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="TwoSample.txt",
	 			n.chains=3, n.iter=10000, n.burnin=5000, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir,
	 			codaPkg=F, debug=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

samples$summary # overview

# Collect posterior samples across all chains:
delta.posterior <- samples$sims.list$delta  

#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(delta.posterior)

# 95% confidence interval:
x0 <- qlogspline(0.025,fit.posterior)
x1 <- qlogspline(0.975,fit.posterior)

posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior     <- dcauchy(0)                   # height of order--restricted prior at delta = 0
BF01      <- posterior/prior
BF01

#============ Plot Prior and Posterior  ===========================
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xlow  <- -3
xhigh <- 3
yhigh <- 2
Nbreaks <- 80
y       <- hist(delta.posterior, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=2,
     xlim=c(xlow,xhigh), ylim=c(0,yhigh), xlab=" ", ylab="Density", axes=F) 
axis(1, at = c(-4,-3,-2,-1,0,1,2,3,4), lab=c("-4","-3","-2","-1","0", "1", "2", "3", "4"))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
par(new=T)
plot ( function( x ) dcauchy( x, 0, 1 ), xlow, xhigh, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lwd=2, lty=1, ylab=" ", xlab = " ", axes=F) 
axis(1, at = c(-4,-3,-2,-1,0,1,2,3,4), lab=c("-4","-3","-2","-1","0", "1", "2", "3", "4"))
axis(2)
points(0, dcauchy(0), pch=19, cex=2)

