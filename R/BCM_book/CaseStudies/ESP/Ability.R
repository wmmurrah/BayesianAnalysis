# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/ESP")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

# Proportion correct on erotic pictures, block 1 and block 2:
prc1.ero <- c(0.6000000, 0.5333333, 0.6000000, 0.6000000, 0.4666667, 
             0.6666667, 0.6666667, 0.4000000, 0.6000000, 0.6000000,
             0.4666667, 0.6666667, 0.4666667, 0.6000000, 0.3333333,
             0.4000000, 0.4000000, 0.2666667, 0.3333333, 0.5333333,
             0.6666667, 0.5333333, 0.6000000, 0.4000000, 0.4666667, 
             0.7333333, 0.6666667, 0.6000000, 0.6666667, 0.5333333,
             0.5333333, 0.6666667, 0.4666667, 0.3333333, 0.4000000,
             0.5333333, 0.4000000, 0.4000000, 0.3333333, 0.4666667,
             0.4000000, 0.4666667, 0.4666667, 0.5333333, 0.3333333,
             0.7333333, 0.2666667, 0.6000000, 0.5333333, 0.4666667,
             0.4000000, 0.5333333, 0.6666667, 0.4666667, 0.5333333,
             0.5333333, 0.4666667, 0.4000000, 0.4666667, 0.6666667,
             0.4666667, 0.3333333, 0.3333333, 0.3333333, 0.4000000,
             0.4000000, 0.6000000, 0.4666667, 0.3333333, 0.3333333,
             0.6666667, 0.5333333, 0.3333333, 0.6000000, 0.4666667,
             0.4666667, 0.4000000, 0.3333333, 0.4666667, 0.5333333,
             0.8000000, 0.4000000, 0.5333333, 0.5333333, 0.6666667,
             0.6666667, 0.6666667, 0.6000000, 0.6000000, 0.5333333,
             0.3333333, 0.4666667, 0.6666667, 0.5333333, 0.3333333,
             0.3333333, 0.2666667, 0.2666667, 0.4666667, 0.6666667)

prc2.ero <- c(0.3333333, 0.6000000, 0.5333333, 0.2666667, 0.6666667,
             0.5333333, 0.6666667, 0.4666667, 0.4666667, 0.6666667,
             0.4000000, 0.6666667, 0.2666667, 0.4000000, 0.4666667,
             0.3333333, 0.5333333, 0.6000000, 0.3333333, 0.4000000,
             0.4666667, 0.4666667, 0.6000000, 0.5333333, 0.5333333,
             0.6000000, 0.5333333, 0.6666667, 0.6000000, 0.2666667,
             0.4666667, 0.4000000, 0.6000000, 0.5333333, 0.4000000,
             0.4666667, 0.5333333, 0.3333333, 0.4000000, 0.4666667,
             0.8000000, 0.6000000, 0.2000000, 0.6000000, 0.4000000,
             0.4000000, 0.2666667, 0.2666667, 0.6000000, 0.4000000,
             0.4000000, 0.4000000, 0.4000000, 0.4000000, 0.6666667,
             0.7333333, 0.5333333, 0.5333333, 0.3333333, 0.6000000,
             0.5333333, 0.5333333, 0.4666667, 0.5333333, 0.4666667,
             0.5333333, 0.4000000, 0.4000000, 0.4666667, 0.6000000,
             0.6000000, 0.6000000, 0.4666667, 0.6000000, 0.6666667,
             0.5333333, 0.4666667, 0.6000000, 0.2000000, 0.5333333,
             0.4666667, 0.4000000, 0.5333333, 0.5333333, 0.5333333,
             0.5333333, 0.6000000, 0.6666667, 0.4000000, 0.4000000,
             0.5333333, 0.8000000, 0.6000000, 0.4000000, 0.2000000,
             0.6000000, 0.6666667, 0.4666667, 0.4666667, 0.4666667)             

x <- matrix(cbind(prc1.ero,prc2.ero),nrow=100) 
n <- nrow(x) # number of participants

data <- list("x", "n") # to be passed on to WinBUGS
myinits <-	list(
  list(r = 0, mu = c(0,0), lambda = c(1,1)))
# parameters to be monitored:	
parameters <- c("r", "mu", "sigma")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="Correlation_1.txt",
	 			n.chains=1, n.iter=5000, n.burnin=5, n.thin=1,
	 			DIC=T, bugs.directory=bugsdir,
	 			codaPkg=F, debug=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

r <- samples$sims.list$r

#Method 1, select only positive r's:
r.pos <- r[which(r>0)]
#============ BFs based on logspline fit ===========================
library(polspline) # this package can be installed from within R
fit.posterior <- logspline(r.pos, lbound=0, ubound=1)

# 95% confidence interval:
x0 <- qlogspline(0.025,fit.posterior)
x1 <- qlogspline(0.975,fit.posterior)

posterior     <- dlogspline(0, fit.posterior) # this gives the pdf at point r = 0
prior         <- 1                            # height of prior at r = 0
BF10.M1       <- prior/posterior
BF10.M1 # 0.36

#Method 2, renormalization:
fit.posterior <- logspline(r)
posterior     <- dlogspline(0, fit.posterior) # this gives the pdf at point delta = 0
# renormalize:
area            <- sum(posterior > 0)/length(posterior)
posterior.OR.M2 <- posterior/area
BF10.M2         <- prior/posterior.OR.M2
BF10.M2 # 0.49

# Compare to exact solution Jeffreys (numerical integration):
BF10.J.exact.positive <- function(n, r)
{
  integrand <- function(rho) {((1-rho^2)^(n/2)) / ((1-rho*r)^(n-.5))}
  BF10      <- integrate(integrand, lower=0, upper=1)$value
  return(BF10)
}
BF10.J.exact.positive(n=length(prc1.ero)-1, r=cor(prc1.ero,prc2.ero)) # 0.46
#====================================================================

#============ Plot Prior and Posterior  ===========================
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
xlow  <- 0
xhigh <- 1
yhigh <- 5
Nbreaks <- 80
y <- hist(r.pos, Nbreaks, prob=T, border="white", ylim=c(0,yhigh), xlim=c(xlow,xhigh), lwd=2, lty=1, 
         ylab=" ", xlab=" ", main=" ", axes=F) 
#white makes the original histogram -- with unwanted vertical lines -- invisible
lines(c(y$breaks, max(y$breaks)), c(0,y$density,0), type="S", lwd=2, lty=1) 
axis(1)
axis(2)
mtext("Positive correlation r", side=1, line = 2.8, cex=2)
mtext("Density", side=2, line = 2.8, cex=2, las=0)
#now bring in log spline density estimation:
par(new=T)
plot(fit.posterior, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F)
points(0, dlogspline(0, fit.posterior),pch=19, cex=2)
points(0, posterior.OR.M2, pch=19, cex=2, col="grey")
# plot the prior:
lines(c(0,1),c(1,1),lwd=2)
points(0, 1, pch=19, cex=2)
###########################################################################
