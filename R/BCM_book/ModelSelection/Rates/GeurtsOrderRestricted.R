#### R2Winbugs Code for a Bayesian hierarchical test of the Geurts data.
#### Do children with ADHD perform worse on the Wisconsin Card Sorting Task
#### than children who develop typically (i.e., normal controls)?
   
rm(list=ls()) # clears workspace

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ModelSelection/Rates")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

### Geurts data:
# Normal Controls:         
num.errors <- c(15,10, 61,11, 60, 44, 63, 70, 57,11, 67, 21, 89,12, 63, 11, 96,10, 37,19, 44,18, 78, 27, 60,14)
nc         <- c(89,74,128,87,128,121,128,128,128,78,128,106,128,83,128,100,128,73,128,86,128,86,128,100,128,79)
kc         <- nc - num.errors
nsc        <- length(kc)
# ADHD:
num.errors <- c(88, 50, 58,17, 40, 18,21, 50, 21, 69,19, 29,11, 76, 46, 36, 37, 72,27, 92,13, 39, 53, 31, 49, 
                57,17,10,12,21, 39, 43, 49,17, 39,13, 68, 24, 21,27, 48, 54, 41, 75, 38, 76,21, 41, 61,24, 28,21)
na         <- c(128,128,128,86,128,117,89,128,110,128,93,107,87,128,128,113,128,128,98,128,93,116,128,116,128,
               128,93,86,86,96,128,128,128,86,128,78,128,111,100,95,128,128,128,128,128,128,98,127,128,93,110,96)
ka         <- na - num.errors
nsa        <- length(ka)
           
data <- list("nc","kc","nsc","na","ka","nsa") # to be passed on to WinBUGS

myinits <- list(
  list(mu = 0, sigma = 1, delta = .3),
  list(mu = -.8, sigma = 2, delta = .5),
  list(mu = .8, sigma = 1.5, delta = .7))

parameters <- c("delta")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="GeurtsOrderRestricted.txt",
	 			n.chains=3, n.iter=100000, n.burnin=1000, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir, #DIC=T or else error message
	 			codaPkg=F, debug=T)
# Now the values for delta are in the "samples" object, ready for inspection.

# Collect posterior samples across all chains:
delta.posterior <- samples$sims.list$delta      

#============ Order-restricted analysis      ===========================
#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(delta.posterior, lbound=0)

posterior <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
prior     <- 2*dnorm(0)                   # height of order--restricted prior at delta = 0
BF01      <- posterior/prior
# BF01 =  5.07, BF10 = 0.20 

#####################################################################
### Plot Prior and Posterior
#####################################################################
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
ymax <- 6
Nbreaks <- 80
y <- hist(delta.posterior, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=2,
     xlim=c(0,3), ylim=c(0,ymax), xlab=" ", ylab="Density", axes=F) 
lines(c(0,0), c(0,4), col="white", lwd=2)
axis(1, at = c(0,1,2,3), lab=c("0", "1", "2", "3"))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,ymax), xlim=c(0,3), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
par(new=T)
plot ( function( x ) 2*dnorm( x, 0, 1 ), 0, 3, xlim=c(0,3), ylim=c(0,ymax), 
       lwd=2, lty=1, ylab=" ", xlab = " ", axes=F) 
points(0, 2*dnorm(0), pch=19, cex=2)

text(0.3, 1.5, labels = "Posterior", cex = 1.5, pos=4)
text(1, .75, labels = "Prior", cex=1.5, pos=4)


