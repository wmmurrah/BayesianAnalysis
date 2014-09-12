# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ModelSelection/Rates")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

# pledger data:
s1 <- 424
s2 <- 5416
n1 <- 777
n2 <- 9072

data  <- list("s1","s2","n1","n2") # to be passed on to WinBUGS

myinits <- list(
  list(thetap = c(-.8,-.4)),
  list(thetap = c(-.5,-.25)),
  list(thetap = c(-.2,-.1)))

parameters <- c("delta", "deltaprior")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(	data, inits=myinits, parameters,
	 			model.file ="Pledgers_2.txt",
	 			n.chains=3, n.iter=20000, n.burnin=1000, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir, #DIC=T or else error message
	 			codaPkg=F, debug=T)
# Now the values for delta and deltaprior are in the "samples" object, ready for inspection.

######################################################
# Order-restriction. H2: delta < 0
######################################################
# Collect posterior samples across all chains:
delta.posterior  <- samples$sims.list$delta      
delta.prior      <- samples$sims.list$deltaprior

#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R

fit.posterior <- logspline(delta.posterior, lbound=-1, ubound=0)
posterior     <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
BF02          <- posterior/2 # 0.26, BF20 = 3.78

#======= Plot Order-Restricted Prior and Posterior ======================
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
Nbreaks <- 20
y <- hist(delta.prior, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(-1,0), ylim=c(0,25), xlab=" ", ylab="Density", axes=F) 
axis(1, at = c(-1, -0.5, 0), lab=c("-1", "-0.5", "0"))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
par(new=T)
x <- hist(delta.posterior, Nbreaks, plot=F)
plot(c(x$breaks, max(x$breaks)), c(0,x$density,0), type="S", lwd=2, lty=2,
     xlim=c(-1,0), ylim=c(0,25), xlab=" ", ylab="Density", main ="Full Scale", axes=F) 
lines(c(0,0), c(0,2), col="white", lwd=2)
axis(1, at = c(-1, -0.5, 0), lab=c("-1", "-0.5", "0"))
axis(2)
#now bring in log spline density estimation:
par(new=T)
# plot the prior:
lines(c(-1,0),c(0,2), lty=1, lwd=1)
par(new=T)
plot(fit.posterior, ylim=c(0,25), xlim=c(-1,0), lty=1, lwd=1, axes=F)

text(-0.42, 20, labels = "Posterior", cex = 1.5, pos=4)
text(-0.5, 2.5, labels = "Prior", cex=1.5, pos=4)

######## Second plot, zoom in:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xmin <- -0.05
xmax <- 0
ymax <- 5
plot(0,0, ylim=c(0,ymax), xlim=c(xmin,xmax), lwd=2, lty=3, ylab="Density", xlab=" ", main="Zoomed in", axes=F, col="white") 
#white makes this invisible
axis(1, at = c(xmin, -0.025, xmax), lab=c(paste(xmin), -0.025, paste(xmax)))
axis(2)
mtext(expression(delta), side=1, line = 2.8, cex=2)
par(new=T)
plot(fit.posterior, ylim=c(0,ymax), xlim=c(xmin,xmax), lty=2, lwd=2, axes=F)
lines(c(-1,0),c(0,2), lty=1, lwd=2)
points(0, 2, pch=19, cex=2)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)

text(-0.04, 1.7, labels = "Prior", cex=1.5, pos=4)
text(-0.037, 4, labels = "Posterior", cex = 1.5, pos=4)

