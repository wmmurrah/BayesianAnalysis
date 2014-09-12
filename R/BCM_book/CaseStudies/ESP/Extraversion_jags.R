# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/ESP")
library(R2jags)

k <- c(36, 32, 36, 36, 28, 40, 40, 24, 36, 36, 28, 40, 28, 
       36, 20, 24, 24, 16, 20, 32, 40, 32, 36, 24, 28, 44,
       40, 36, 40, 32, 32, 40, 28, 20, 24, 32, 24, 24, 20, 
       28, 24, 28, 28, 32, 20, 44, 16, 36, 32, 28, 24, 32,
       40, 28, 32, 32, 28, 24, 28, 40, 28, 20, 20, 20, 24,
       24, 36, 28, 20, 20, 40, 32, 20, 36, 28, 28, 24, 20,
       28, 32, 48, 24, 32, 32, 40, 40, 40, 36, 36, 32, 20,
       28, 40, 32, 20, 20, 16, 16, 28, 40)
       
x <- c(50, 80, 79, 56, 50, 80, 53, 84, 74, 67, 50, 45, 62, 
       65, 71, 71, 68, 63, 67, 58, 72, 73, 63, 54, 63, 70, 
       81, 71, 66, 74, 70, 84, 66, 73, 78, 64, 54, 74, 62, 
       71, 70, 79, 66, 64, 62, 63, 60, 56, 72, 72, 79, 67, 
       46, 67, 77, 55, 63, 44, 84, 65, 41, 62, 64, 51, 46,
       53, 26, 67, 73, 39, 62, 59, 75, 65, 60, 69, 63, 69, 
       55, 63, 86, 70, 67, 54, 80, 71, 71, 55, 57, 41, 56, 
       78, 58, 76, 54, 50, 61, 60, 32, 67)
       
nsubjs  <- length(k)
ntrials <- 60
sigmax  <- 3
lambdax <- 1/(sigmax^2)

data <- c("k", "x", "lambdax", "ntrials", "nsubjs") # to be passed on to JAGS
myinits <-	list(
  list(r = -0.1, mu = c(0.1,-0.1), lambda = c(.5,1)),
  list(r = 0.1, mu = c(0.2,-0.2), lambda = c(1,.5))
  )
# parameters to be monitored:	
parameters <- c("r", "mu", "sigma", "theta")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			model.file ="Extraversion.txt",
	 			n.chains=2, n.iter=10000, n.burnin=5000, n.thin=1)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

r <- samples$BUGSoutput$sims.list$r

#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(r)

# 95% confidence interval:
x0 <- qlogspline(0.025,fit.posterior)
x1 <- qlogspline(0.975,fit.posterior)

posterior     <- dlogspline(0, fit.posterior) # this gives the pdf at point r = 0
prior         <- .5                           # height of prior at r = 0
BF10          <- prior/posterior
BF10
# 0.37

#============ Plot Prior and Posterior  ===========================
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xlow  <- -1
xhigh <- 1
yhigh <- 5
Nbreaks <- 80

y <- hist(r, Nbreaks, plot=F)
plot(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1,
     xlim=c(xlow,xhigh), ylim=c(0,yhigh), xlab=" ", ylab=" ", main = " ", axes=F) 
axis(1)
axis(2)
mtext("Correlation r", side=1, line = 2.8, cex=2)
mtext("Density", side=2, line = 2.8, cex=2, las=0)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
# plot the prior:
lines(c(-1,1),c(0.5,0.5),lwd=2)
points(0, .5, pch=19, cex=2)
###########################################################################

