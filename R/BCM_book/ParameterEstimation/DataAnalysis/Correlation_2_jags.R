# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/DataAnalysis")

library(R2jags)

x <- matrix(c(.8,102, 1,98, .5,100, 0.9,105, .7,103, 
             0.4,110, 1.2,99, 1.4,87, 0.6,113, 1.1,89, 1.3,93),
             nrow=11,ncol=2,byrow=T) 
n <- nrow(x) # number of people/units measured

#precision of measurement:
sigmaerror <- c(.03, 1)
# sigmaerror = c(.03, 10)
lambdaerror <- 1/sigmaerror^2

data <- list("x", "n", "lambdaerror") # to be passed on to JAGS
myinits <-	list(
  list(r = 0, mu = c(0,0), lambda = c(1,1)))
# parameters to be monitored:	
parameters <- c("r", "mu", "sigma")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			 model.file ="Correlation_2.txt", n.chains=1, n.iter=5000, 
         n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

r <- samples$BUGSoutput$sims.list$r

#Frequentist point-estimate of r:
freq.r <- cor(x[,1],x[,2])

#make the two panel plot:
windows(width=9,height=6) #this command works only under Windows!
layout(matrix(c(1,2),1,2))
layout.show(2)
#some plotting options to make things look better:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
# data panel:    
plot(x[,1],x[,2], type="p", pch=19, cex=1)
for (i in 1:n)
{
  lines(c(x[i,1]-sigmaerror[1],x[i,1]+sigmaerror[1]), c(x[i,2],x[i,2]))
  lines(c(x[i,1],x[i,1]), c(x[i,2]-sigmaerror[2],x[i,2]+sigmaerror[2]))
}
# correlation panel:
plot(density(r, from=-1,to=1), main="", ylab="Posterior Density", xlab="Correlation", lwd=2)
lines(c(freq.r, freq.r), c(0,100), lwd=2, lty=2)

