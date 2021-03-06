# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/Simple")
library(R2jags)

# read the data
y          <- matrix(scan("k_M.txt", sep=","), ncol=40, nrow=6, byrow=T) 
n          <- c(1440,1280,1520,1520,1200,1280)
listlength <- c(10,15,20,20,30,40)
pc         <- matrix(scan("pc_M.txt", sep=","), ncol=40, nrow=6, byrow=T) 
labs       <- scan("labs_M.txt", sep="", what="character")
dsets <- 6 
gsets <- 9

# Set Dataset to Use
w <- array()
l <- array()
m <- matrix(rep(0,50*9), ncol=50, nrow=9, byrow=T) 
for (dset in 1:gsets)
{
  if (dset==1)
  {
    nwords <- 10
    lag    <- 2
    offset <- 15
  } 
  if (dset==2)
  {
    nwords <- 15
    lag    <- 2
    offset <- 20
  } 
  if (dset==3)
  {
    nwords <- 20
    lag    <- 2
    offset <- 25
  } 
  if (dset==4)
  {
    nwords <- 20
    lag    <- 1
    offset <- 10
  } 
  if (dset==5)
  {
    nwords <- 30
    lag    <- 1
    offset <- 15
  } 
  if (dset==6)
  {
    nwords <- 40
    lag    <- 1
    offset <- 20
  }
  #Generalization: 
  if (dset==7)
  {
    nwords <- 10
    lag    <- 1
    offset <- 5
  } 
  if (dset==8)
  {
    nwords <- 25
    lag    <- 1
    offset <- 12.5
  } 
  if (dset==9)
  {
    nwords <- 50
    lag    <- 1
    offset <- 25
  } 
  # Temporal Offset For Free Recall
  m[dset,1:nwords] <- offset+seq(from=(nwords-1)*lag, by=-lag, to=0)
  w[dset]          <- nwords 
  l[dset]          <- lag
  listlength[dset] <- nwords
}

n[(dsets+1):gsets] <- 1200
labs = c(labs,c("10-1","25-1","50-1"))

y <- t(y)
m <- t(m)
data <- list("gsets", "y", "n", "listlength", "m", "w", "dsets") # to be passed on to JAGS

myinits <-	list(
  list(c = 20.5, a = c(-.003, .63), s = 9.5))  
 
# parameters to be monitored:	
parameters <- c("c","s","t","predpc[1:50,1:9]","a")

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters,
	 			       model.file ="SIMPLE_2R_jags.txt", n.chains=1, 
               n.iter=2000, n.burnin=1000, n.thin=1)
# NB. The predpc values of -999 are obviously dummy values;
# they should not appear in any of the plots.
# NB2. For better convergence, you might want to take a lunch break and run
# multiple chains, for more iterations, and do some thinning.
# NB3. Sometimes JAGS complains about initial values. Even setting the initial values
# to the posterior means from a WinBUGS run did not help. The new JAGS file makes the
# model more robust by having an upper bound for theta of .999 instead of 1.

###Figure 15.5

layout(matrix(c(
	10,1,2,3,
	10,4,5,6,
      10,7,8,9,
	11,11,11,11
	),4,4,byrow=T), c(1,2,2,2), c(2,2,2,1))
layout.show(11)
hm <- 20
ll <- listlength

for (dset in 1:gsets) {
	plot(-1,-1,xlim=c(0,50),ylim=c(0,1),xlab="",ylab="",las=1)
	text(47,.96,labs[dset])
	for (i in 1:ll[dset]) {                 
		data <- samples$BUGSoutput$sims.array[,,paste("predpc[",i,",",dset,"]",sep="")]
		points(i+runif(hm,0,1)*.1,data[ceiling(runif(hm,0,1)*samples$BUGSoutput$n.iter)],col="grey")
	}
	if (dset <= 6) {
		points(1:ll[dset],pc[dset,1:ll[dset]],xlim=c(0,50),ylim=c(0,1))
		lines(1:ll[dset],pc[dset,1:ll[dset]])
	}
	box("plot")
}
par(mar=c(rep(0,4)))
plot.new()
text(.45,.5,"Probability Correct",cex=2.5,srt=90)
plot.new()
text(.5,.5,"Serial Position",cex=2.5,mar=c(rep(0,4)))


###Figure 15.6
layout(matrix(1:3,1,3))
layout.show(3)
#Threshold
epss <- .05
sx  <- seq(9,11,epss)
sxe <- seq(9+epss/2,11-epss/2,epss)
S <- samples$BUGSoutput$sims.array[,,"s"]
count <- hist(S,breaks=sx,plot=F)
count <- count$counts
count <- count/sum(count)*epss
keep <- which(count>1e-12)

plot(sxe[keep],count[keep],type="l",ylim=c(0,0.015),xlim=c(9,11),xlab="Threshold Noise (s)", ylab="Posterior Density", cex.lab=1.2,axes=F)
axis(1)
axis(2,labels=F,lwd.ticks=0)
box("plot")

#Distinctiveness
epsc <- .17
cx <- seq(18,24,epsc) #mids are cxe in R
cxe <- seq(18+epsc/2,24-epsc/2,epsc)
C <- samples$BUGSoutput$sims.array[,,"c"]
count <- hist(C,breaks=cx,plot=F)
count <- count$counts
count <- count/sum(count)*epsc
keep <- which(count>1e-12)

plot(cxe[keep],count[keep],type="l",ylim=c(0,0.06),xlim=c(18,23),xlab="Distinctiveness (c)", ylab="Posterior Density", cex.lab=1.2, axes=F)
axis(1)
axis(2,labels=F,lwd.ticks=0)
box("plot")

###Threshold Parameter as a Function Of List Length
howmany <- 50
nsamples <- 1000
keep <- ceiling(runif(howmany,min=0,max=1)*nsamples)
wdom <- seq(1,50,1)
plot(-1, -1, xlim=c(1,50), ylim=c(0,1),xlab="Item List Length (W)", ylab="Treshold (t)")
for (i in 1:howmany){
	predt <- samples$BUGSoutput$sims.array[,,"a[1]"][[keep[i]]]*wdom+samples$BUGSoutput$sims.array[,,"a[2]"][[keep[i]]]
	points(wdom,predt,col="grey")
}
