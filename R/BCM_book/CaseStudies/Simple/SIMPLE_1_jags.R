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

dsets <- 6 

# Loop over conditions
m <- y*0
for (dset in 1:dsets)
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
 # Temporal Offset For Free Recall
  m[dset,1:nwords] <- offset+seq(from=(nwords-1)*lag, by=-lag, to=0)
}

y <- t(y)
m <- t(m)
data <- list("y", "n", "listlength", "m", "dsets") # to be passed on to JAGS

#myinits <-	list(
#  list(c = c(14.16,17.38,17.42,15.41,17.98,21.55), t = c(.60,.59,.55,.54,.50,.47), 
#       s = c(7.61,8.07,9.43,11.36,11.73,14.23)))  

myinits <-	list(
  list(c = rep(15,dsets), t = rep(.5,dsets), s = rep(10,dsets)))  

# parameters to be monitored:
parameters <- c("c", "s", "t", "predpc[1:40,1:6]")
 # Monitoring parameter "predpc" in SIMPLE_1.txt gives an error upon return to R. 
# The following code (SIMPLE_1R.txt) fixes this by assigning dummy values to 
# unused entries in the predpc matrix.

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples <- jags(data, inits=myinits, parameters, model.file ="SIMPLE_1R_jags.txt",
	 			   n.chains=1, n.iter=2000, n.burnin=1000, n.thin=1)
# NB. The predpc values of -999 are obviously dummy values;
# they should not appear in any of the plots.
# NB2. For better convergence, you might want to take a lunch break and run
# multiple chains, for more iterations, and do some thinning.
# NB3. Sometimes JAGS complains about initial values. Even setting the initial values
# to the posterior means from a WinBUGS run did not help. The new JAGS file makes the
# model more robust by having an upper bound for theta of .999 instead of 1.

#Figure 15.2

layout(matrix(c(
	7,1,2,3,
	7,4,5,6,
	8,8,8,8
	),3,4,byrow=T), c(1,2,2,2), c(2,2,.5))
layout.show(8)
hm <- 20
ll <- listlength

for (dset in 1:dsets) {
	plot(-1,-1,xlim=c(0,40),ylim=c(0,1),xlab="",ylab="",las=1)
	for (i in 1:ll[dset]) { 
		data <- samples$BUGSoutput$sims.array[,,paste("predpc[",i,",",dset,"]",sep="")]
		points(i+runif(hm,0,1)*.1,data[ceiling(runif(hm,0,1)*samples$BUGSoutput$n.iter)],col="grey")
	}
	points(1:ll[dset],pc[dset,1:ll[dset]],xlim=c(0,40),ylim=c(0,1))
	lines(1:ll[dset],pc[dset,1:ll[dset]])

	box("plot")
}
par(mar=c(rep(0,4)))
plot.new()
text(.45,.5,"Probability Correct",cex=2.5,srt=90)
plot.new()
text(.5,.5,"Serial Position",cex=2.5,mar=c(rep(0,4)))

