# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/PsychophysicalFunctions")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

x       <- matrix(NA,8,28)
x[]     <- as.matrix(read.table("data_x.txt",sep="\t"))
n       <- matrix(NA,8,28)
n[]     <- as.matrix(read.table("data_n.txt",sep="\t"))
r       <- matrix(NA,8,28)
r[]     <- as.matrix(read.table("data_r.txt",sep="\t"))
rprop   <- matrix(NA,8,28)
rprop[] <- as.matrix(read.table("data_rprop.txt",sep="\t"))

xmean <- c(318.888,311.0417,284.4444,301.5909,296.2000,305.7692,294.6429,280.3571)
nstim <- c(27, 24, 27, 22, 25, 26, 28, 28)
nsubjs <- 8

data <- list("x","xmean","n","r","nsubjs","nstim") # to be passed on to WinBUGS
myinits <-	list(
  list(alpha = runif(nsubjs,-2,2), beta = runif(nsubjs,0,.5)),
  list(alpha = runif(nsubjs,-2,2), beta = runif(nsubjs,0,.5)),
  list(alpha = runif(nsubjs,-2,2), beta = runif(nsubjs,0,.5)))

################################################
# Part for model without contamination parameter
################################################

# parameters to be monitored:	
parameters <- c("alpha","beta")

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters,
	 			model.file ="Psychophysical_1.txt",
	 			n.chains=3, n.iter=10000, n.burnin=5000, n.thin=1, 
	 			DIC=T, bugs.directory=bugsdir,
	 			codaPkg=F, debug=T)

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

# Extracting the parameters
alpha	    <- samples$sims.list$alpha
beta  	  <- samples$sims.list$beta
alphaMAP  <- c(rep(0,nsubjs))
betaMAP   <- c(rep(0,nsubjs))
alpha_sel <- matrix(NA,20,8) 
beta_sel  <- matrix(NA,20,8) 

# Constructing MAP-estimates and alpha/beta range
for (i in 1:nsubjs)
{
	alphaMAP[i]   <- density(alpha[,i])$x[which(density(alpha[,i])$y==max(density(alpha[,i])$y))]
	betaMAP[i]    <- density(beta[,i])$x[which(density(beta[,i])$y==max(density(beta[,i])$y))]
	alpha_sel[,i] <- sample(alpha[,i],20)
	beta_sel[,i]  <- sample(beta[,i],20)
}

############################## PSYCHOMETRIC FUNCTIONS ##############################

F1 <- function(X,s) # only the MAP estimate; use this to plot psychometric functions
{
  exp(alphaMAP[s] + betaMAP[s]*(X - xmean[s]))/(1+exp(alphaMAP[s] + betaMAP[s]*(X - xmean[s])))
}

F1inv <- function(Y,s)
{
  (log(-Y/(Y-1))-alphaMAP[s])/betaMAP[s]
}

F2 <- function(X,s) # function for all the posterior alpha/beta values; use this to calculate JND posterior
{
  exp(alpha[,s] + beta[,s]*(X - xmean[s]))/(1+exp(alpha[,s] + beta[,s]*(X - xmean[s])))
}
F2inv <- function(Y,s)
{
  (log(-Y/(Y-1))-alpha[,s])/beta[,s]
}

F3 <- function(X,s,g) # function for 20 grabbed posterior alpha/beta values; use this to plot overlapping sigmoids to visualize variance
{
  exp(alpha_sel[g,s] + beta_sel[g,s]*(X - xmean[s]))/(1+exp(alpha_sel[g,s] + beta_sel[g,s]*(X - xmean[s])))
}

##################################### JND/PSE calculation ########################################
JND 	 <- F2inv(0.84,c(1:nsubjs))-F2inv(0.5,c(1:nsubjs))
JNDmap <- F1inv(0.84,c(1:nsubjs))-F1inv(0.5,c(1:nsubjs))
								  				             
PSE 	 <- F2inv(0.5,c(1:nsubjs))+xmean
PSEmap <- F1inv(0.5,c(1:nsubjs))+xmean
		
################## PLOTS ####################

### Figure 12.2

dev.new(width=10,height=5)
layout(matrix(1:nsubjs,2,4,byrow=T))
par(mar=c(1,2,2,0),oma=c(5,5,1,1))
for (i in 1:nsubjs)
{
	scale <- seq(x[i,1],x[i,nstim[i]], by=.1)
	plot(x[i,],rprop[i,],main=paste("Subject",as.character(i)),xlab="",ylab="",pch=15,col="dark grey",ylim=c(0,1),yaxt="n",xaxt="n")
	lines(scale,F1(scale,i),type="l")
	segments(x0=x[i,1],x1=PSEmap[i]+JNDmap[i],y0=0.84,lty=2)
	segments(x0=x[i,1],x1=PSEmap[i],y0=0.5,lty=2)
	segments(y0=0,y1=0.84,x0=PSEmap[i]+JNDmap[i],lty=2)
	segments(y0=0,y1=0.5,x0=PSEmap[i],lty=2)
	if (i==1 | i==5) 
  {
		axis(2,las=1,yaxp=c(0,1,2))
		axis(2,at=0.84,las=1)
	}
	if (i>4) axis(1)
}
mtext("Proportion 'Long' Response",side=2,line=2,outer=T,cex=1.4)
mtext("Test Interval (ms)",side=1,outer=T,line=3,cex=1.4)

### WARNING: Do not close R window.

### NOTE:	 Answers to the exercises can be found in PsychometricFunction1_Answers.R
