# The original WinBUGS code gave "undefined real result" -- the JAGS equivalent does work,
# and the WinBUGS also works when called from Matlab. Alexander Ly discovered that the 
# censoring on node delta does not work using the I(,0) operator; hence the 
# trick, implemented in OneSampleOrderRestricted_AL.txt, to multiply delta by
# -1 and use the I(0,) operator instead.

# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ModelSelection/Means")

library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

# Read data Dr. Smith
Winter <- c(-0.05,0.41,0.17,-0.13,0.00,-0.05,0.00,0.17,0.29,0.04,0.21,0.08,0.37,0.17,0.08,-0.04,-0.04,0.04,-0.13,-0.12,0.04,0.21,0.17,
       0.17,0.17,0.33,0.04,0.04,0.04,0.00,0.21,0.13,0.25,-0.05,0.29,0.42,-0.05,0.12,0.04,0.25,0.12)
 
Summer <- c(0.00,0.38,-0.12,0.12,0.25,0.12,0.13,0.37,0.00,0.50,0.00,0.00,-0.13,-0.37,-0.25,-0.12,0.50,0.25,0.13,0.25,0.25,0.38,0.25,0.12,
      0.00,0.00,0.00,0.00,0.25,0.13,-0.25,-0.38,-0.13,-0.25,0.00,0.00,-0.12,0.25,0.00,0.50,0.00)

x <- Winter-Summer # allowed because it is a within-subjects design
x <- x/sd(x)       # standardize

ndata <- length(Winter) # number of subjects

data  <- list("x", "ndata") # to be passed on to WinBUGS

myinits <- list(
  list(deltatmp = runif(1), sigmatmp = runif(1)),
  list(deltatmp = runif(1), sigmatmp = runif(1)),
  list(deltatmp = runif(1), sigmatmp = runif(1)))

# Parameters to be monitored
parameters <- c("delta")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file="OneSampleOrderRestricted_AL.txt",
	 			n.chains=3, n.iter=10000, n.burnin=1000, n.thin=1,
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

#============ BFs based on logspline fit ===========================
fit.posterior <- logspline(delta.posterior,ubound=0) # NB. note the bound

# 95% confidence interval:
x0 <- qlogspline(0.025,fit.posterior)
x1 <- qlogspline(0.975,fit.posterior)

posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior     <- 2*dcauchy(0)                 # height of order--restricted prior at delta = 0
BF01      <- posterior/prior
BF01

#============ Plot Prior and Posterior  ===========================
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xlow  <- -3
xhigh <- 0
yhigh <- 12
Nbreaks <- 80
y <- hist(delta.posterior, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=2,
     xlim=c(xlow,xhigh), ylim=c(0,yhigh), xlab=" ", ylab="Density", axes=F) 
axis(1, at = c(-3,-2,-1,0), lab=c("-3","-2","-1","0"))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
par(new=T)
plot ( function( x ) 2*dcauchy( x, 0, 1 ), xlow, xhigh, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lwd=2, lty=1, ylab=" ", xlab = " ", axes=F) 
axis(1, at = c(-3,-2,-1,0), lab=c("-3","-2","-1","0"))
axis(2)
points(0, 2*dcauchy(0), pch=19, cex=2)

