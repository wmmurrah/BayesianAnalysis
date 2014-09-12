rm(list=ls()) 
library(R2WinBUGS)
library(MCMCpack)

setwd("C:/Dropbox/My Documents/MPT_Project/S&B/Chapter/Hier")
bugsdir = "C:/Program Files/WinBUGS14"
working.directory=getwd()

### Riefer et al (2002) data:
trial_1 = list( subjs = 21, items=20, response = structure(.Data = c(2,4,4,10,2,1,3,14,2,2,5,11,6,0,4,10,1,0,4,15,1,0,2,17,1,2,4,13,4,1,6,9,5,1,4,10,1,0,9,10,5,0,3,12,0,1,6,13,1,5,7,7,1,1,4,14,2,2,3,13,2,1,5,12,2,0,6,12,1,0,5,14,2,1,8,9,3,0,2,15,1,2,3,14),.Dim = c(4, 21)))
trial_2 = list( subjs = 21, items=20, response = structure(.Data = c(7,5,3,5,5,2,3,10,6,2,7,5,9,4,2,5,2,2,7,9,1,3,3,13,5,0,5,10,7,3,4,6,7,3,6,4,4,1,10,5,9,1,2,8,3,1,6,10,3,5,9,3,2,0,6,12,8,0,3,9,3,2,7,8,7,1,5,7,2,1,6,11,5,3,5,7,5,0,6,9,6,2,2,10),.Dim = c(4, 21)))
trial_6 = list( subjs = 21, items=20, response = structure(.Data = c(14,3,1,2,12,3,1,4,18,0,1,1,15,3,0,2,7,1,10,2,3,6,11,0,8,4,3,5,17,1,1,1,13,4,3,0,11,6,1,2,16,1,2,1,10,1,3,6,7,13,0,0,8,4,3,5,16,1,1,2,5,4,7,4,15,0,5,0,6,3,6,5,17,2,0,1,17,1,0,2,8,3,6,3),.Dim = c(4, 21)))

response_1 = t(trial_1$response)
response_2 = t(trial_2$response)
response_6 = t(trial_6$response)

I = 21		# Number of participant
J1 = rep(20,I) 	# Number of word pairs per participant	
W = diag(3)		# Identity matrix for Invere-Wishart
P = 3			# Number of free parameters per participant: c_i, r_i, u_i 

###Trial 1
n1 = response_1
data  = list("n1","J1","I","W","P")  

inits=function()
{
    	list(T.prec=rwish(P+1,diag(P)),xi=runif(P), mu=rnorm(P),delta.raw=array(rnorm(I*P),c(I,P)))                           
}

parameters = c("mu","sigma","rho")

samples_1 = bugs(data, inits, parameters,
	 			model.file ="LatentTraitMPT.txt",
	 			n.chains=3, n.iter=20000, n.burnin=2000, n.thin=3,
	 			DIC=T, bugs.directory=bugsdir, #DIC=T or else error message
	 			codaPkg=F, debug=T,working.directory=working.directory)


###Trial 2
n1 = response_2
data  = list("n1","J1","I","W","P")  

samples_2 = bugs(data, inits, parameters,
	 			model.file ="LatentTraitMPT.txt",
	 			n.chains=3, n.iter=20000, n.burnin=2000, n.thin=3,
	 			DIC=T, bugs.directory=bugsdir, #DIC=T or else error message
	 			codaPkg=F, debug=T,working.directory=working.directory)


###Trial 6
n1 = response_6
data  = list("n1","J1","I","W","P")  

samples_6 = bugs(data, inits, parameters,
	 			model.file ="LatentTraitMPT.txt",
	 			n.chains=3, n.iter=20000, n.burnin=2000, n.thin=3,
	 			DIC=T, bugs.directory=bugsdir, #DIC=T or else error message
	 			codaPkg=F, debug=T,working.directory=working.directory)


#### Plots posteriors of the group--level c, r, and u parameters

#postscript("Posteriors_Riefer_Hierarchical.eps",paper = "special" ,width=5, height=2.5,horizontal=F)
layout(matrix(1:3,1,3,byrow=T))
par(mar = c(4.25, 3, 1, 1))
par(cex.axis=0.9)
plot(density(pnorm(samples_1$sims.list$mu[,1])),ylim = c(0,10),xlim=c(0,1), axes=F,xlab="c",ylab="",main="")
axis(1,at=seq(0,1,0.5))
axis(2,at = c(0,12),labels=F,lwd.ticks=0)
mtext("Density",2,1,cex=0.55,las=0)
par(cex=0.60)
legend(0.3, 10, c("Trial1","Trial 2","Trial 3"),lty = c(1,2,3),col=c("black"),text.col = "black")
par(cex=0.65)
lines(density(pnorm(samples_2$sims.list$mu[,1])),lty=2)
lines(density(pnorm(samples_6$sims.list$mu[,1])),lty=3)

plot(density(pnorm(samples_1$sims.list$mu[,2]),to=0.75),ylim = c(0,13),xlim=c(0,1), axes=F,xlab="r",ylab="",main="")
axis(1,at=seq(0,1,0.5))
axis(2,at = c(0,13),labels=F,lwd.ticks=0)
mtext("Density",2,1,cex=0.55,las=0)
lines(density(pnorm(samples_2$sims.list$mu[,2])),lty=2)
lines(density(pnorm(samples_6$sims.list$mu[,2])),lty=3)

plot(density(pnorm(samples_1$sims.list$mu[,3])),ylim = c(0,7.5),xlim=c(0,1), axes=F,xlab="u",ylab="",main="")
axis(1,at=seq(0,1,0.5))
axis(2,at = c(0,10),labels=F,lwd.ticks=0)
mtext("Density",2,1,cex=0.55,las=0)
lines(density(pnorm(samples_2$sims.list$mu[,3])),lty=2)
lines(density(pnorm(samples_6$sims.list$mu[,3])),lty=3)
#dev.off()

















