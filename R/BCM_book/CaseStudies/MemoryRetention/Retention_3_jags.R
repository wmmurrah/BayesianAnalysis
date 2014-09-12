# The original JAGS code is very similar to the WinBUGS code. It works fine
# when called from Matlab. When called from R, however, the JAGS code does not 
# produce sufficient skrinkage for the missing participant. For some reason, 
# replacing the T(0,1) constraint on the alpha and beta by T(0,) alleviates 
# the problem. This is the difference between Retention_3J.txt and 
# Retention_3_jags.txt.
# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/MemoryRetention")
library(R2jags)

t     <- c(1, 2, 4, 7, 12, 21, 35, 59, 99, 200)
nt    <- length(t)
slist <- 1:4
ns    <- length(slist)

k <- matrix(c(18, 18, 16, 13, 9, 6, 4, 4, 4, NA,
             17, 13,  9,  6, 4, 4, 4, 4, 4, NA,
             14, 10,  6,  4, 4, 4, 4, 4, 4, NA,
             NA, NA, NA, NA,NA,NA,NA,NA,NA, NA), nrow=ns, ncol=nt, byrow=T)
k

n <- 18

data <- list("k", "n", "t", "ns", "nt") # to be passed on to JAGS
myinits <-	list(
  list(alphamu = 0.5, alphalambda = 1, betamu = 0.5, betalambda = 1, alpha = rep(0.5,ns), beta = rep(0.1,ns)))  

# parameters to be monitored:	
parameters <- c("alpha", "beta", "predk")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			model.file ="Retention_3_jags.txt",
	 			n.chains=1, n.iter=10000, n.burnin=1, n.thin=1, DIC=T)
# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

##Figure 10.8
n.iter <- 10000
keepi <- 500
keep <- sample(n.iter, keepi)

alpha1 <- samples$BUGSoutput$sims.array[,,1]
alpha2 <- samples$BUGSoutput$sims.array[,,2]
alpha3 <- samples$BUGSoutput$sims.array[,,3]
alpha4 <- samples$BUGSoutput$sims.array[,,4]

beta1 <- samples$BUGSoutput$sims.array[,,5]
beta2 <- samples$BUGSoutput$sims.array[,,6]
beta3 <- samples$BUGSoutput$sims.array[,,7]
beta4 <- samples$BUGSoutput$sims.array[,,8]
d.beta1 <- density(beta1)
d.beta2 <- density(beta2)
d.beta3 <- density(beta3)
d.beta4 <- density(beta4)

layout(matrix(c(1,2,3,0),2,2,byrow=T), width=c(2/3, 1/3), heights=c(2/3,1/3))
#layout.show()

par(mar=c(2,2,1,0))
plot(alpha1[keep],beta1[keep], xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes=F)
points(alpha2[keep],beta2[keep], col="red")
points(alpha3[keep],beta3[keep], col="green")
points(alpha4[keep],beta4[keep],col="blue")
box(lty=1)

par(mar=c(2,1,1,4))
plot(d.beta1$y, d.beta1$x, ylim=range(c(0,1)), xlim=c(12,0),type='l', axes=F, xlab="", ylab="")
#plot(d.beta1$y, d.beta1$x, ylim=range(c(0,1)), xlim=rev(range(d.beta1$y)),type='l', axes=F, xlab="", ylab="")
lines(d.beta2$y, d.beta2$x, col="red")
lines(d.beta3$y, d.beta3$x, col="green")
lines(d.beta4$y, d.beta4$x, col="blue")
axis(4, at=c(0,1))
mtext(expression(beta), side=4,line=1, cex=1.3)
box(lty=1)

par(mar=c(6,2,0,0))
plot(density(alpha1),zero.line=F ,main="", ylab="", xlab="", cex.lab=1.3,xlim=c(0,1), axes=F)
lines(density(alpha2), col="red")
lines(density(alpha3), col="green")
lines(density(alpha4),col="blue")
axis(1,at=c(0,1))
mtext(expression(alpha), side=1.2,line=1, cex=1.3)
box(lty=1)

##Figure 10.9
#close previous graph window before running this code!
layout(matrix(c(1:4),2,2,byrow=T))
#layout.show()
sc <- 3.5
jj <- numeric()
xx <- numeric()

for (i in 1:ns) {
	plot(-1,100,xlim=c(0,10),ylim=c(0,18), main=(paste("Subject", i)),xlab=("Time Lags"), ylab=("Retention Count"),cex.lab=1.3, axes=F)
	axis(1, at=c(1,2,3,4,5,6,7,8,9,10), lab=c("1","2","4","7","12","21","35","59","99","200"),cex.axis=0.7)
	axis(2, at=c(0,18),lab=c("0","18"),cex.axis=0.7)
	box(lty=1)
	for (j in 1:nt) {
		j1 <- i+9 + (j-1)*4
		count <- hist(samples$BUGSoutput$sims.array[,,j1],c(0:n),plot=F)
		count <- count$counts
		count <- count/sum(count)
		for (x in 1:n){
			if (count[x]>0){
				points(j,x,pch=22, col="black",cex=sc*sqrt(count[x]))
				if (!is.na(k[i,j]) && k[i,j]==x){
					points(j,x,pch=22,bg="black",cex=sc*sqrt(count[x]))
					jj <- c(jj,j)
					xx <- c(xx,x)
				}
			}
		}
	}
	coords <- list(x=jj, y=xx)
	lines(coords,lwd=2)
	jj <- numeric()
	xx <- numeric()
}

