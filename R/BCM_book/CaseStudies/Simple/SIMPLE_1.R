# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/Simple")
library(R2WinBUGS)
bugsdir <- "C:/Program Files/WinBUGS14"

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
data <- list("y", "n", "listlength", "m", "dsets") # to be passed on to WinBUGS

myinits <-	list(
  list(c = rep(15,dsets), t = rep(0.5,dsets), s = rep(10,dsets)))  

# parameters to be monitored:
parameters <- c("c", "s", "t", "predpc")
 # Monitoring parameter "predpc" in SIMPLE_1.txt gives an error upon return to R. 
# The following code (SIMPLE_1R.txt) fixes this by assigning dummy values to 
# unused entries in the predpc matrix.

# The following command calls WinBUGS with specific options.
# For a detailed description see Sturtz, Ligges, & Gelman (2005).
samples <- bugs(data, inits=myinits, parameters, model.file ="SIMPLE_1R.txt",
	 			   n.chains=1, n.iter=2000, n.burnin=1000, n.thin=1,
	 			   bugs.directory=bugsdir, codaPkg=F, debug=T, DIC=T)
# NB. The predpc values of -999 are obviously dummy values;
# they should not appear in any of the plots.
# NB2. For better convergence, you might want to take a lunch break and run
# multiple chains, for more iterations, and do some thinning.

round(samples$summary,3)

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
		data <- samples$sims.array[,,paste("predpc[",i,",",dset,"]",sep="")]
		points(i+runif(hm,0,1)*.1,data[ceiling(runif(hm,0,1)*samples$n.iter)],col="grey")
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

