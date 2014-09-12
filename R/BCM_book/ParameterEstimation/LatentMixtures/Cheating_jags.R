# This code makes JAGS complain but the results do match those of WinBUGS.
# Part of this code is courtesy of Wolf Vanpaemel.                                          
# clears workspace:  
rm(list=ls()) 

setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/ParameterEstimation/LatentMixtures")
library(R2jags)

cheat.dat  <- read.table("cheat.csv",header=F,sep=",")
cheatt.dat <- read.table("cheatt.csv",header=F,sep="")
truth <- cheatt.dat$V1 #truth = 1 if cheater
k <- apply(cheat.dat,1,sum) # total correct per participant
p <- length(k) # number of people
n <- 40        # total trials

data <- list("p", "k", "n", "truth") # to be passed on to JAGS
myinits <- list(
  list(z = round(runif(p)), mudiff=0.1, phi=0.5, mubon=0.5, lambdabon=30, lambdache=25),
  list(z = round(runif(p)), mudiff=0.15, phi=0.5, mubon=0.5, lambdabon=25, lambdache=30)
  ) 

#myinits <- list(
#  list(z = round(runif(p)), mudiff=0.1, phi=0.5, mubon=0.5, lambdabon=10, lambdache=10),
#  list(z = round(runif(p)), mudiff=0.2, phi=0.5, mubon=0.5, lambdabon=10, lambdache=10)
#  ) 

# parameters to be monitored:	
parameters <- c("theta","z","mubon","lambdabon",
               "muche","lambdache","mudiff","phi","alpha","beta","pc")

set.seed(3) # some chains result in "undefined real result -- see box & exercise"

# The following command calls JAGS with specific options.
# For a detailed description see the R2jags documentation.
samples = jags(data, inits=myinits, parameters,
               model.file ="Cheating.txt", n.chains=2, n.iter=6000, 
               n.burnin=3000, n.thin=1, DIC=T)
               
samples$summary
pc <- samples$BUGSoutput$sims.list$pc/p #to get proportion correct
mean(pc)

# plot 6.9
#make the two panel plot:
windows(width=8,height=6) #this command works only under Windows!
layout(matrix(c(1,2),2,1))
layout.show(2)
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
bins <- c(-1:n)+.5
bonafide <- hist(k[truth==0], breaks=bins, plot=F)$counts
cheat    <- hist(k[truth==1], breaks=bins, plot=F)$counts

counts <- rbind(bonafide, cheat)
barplot(counts, main=" ", xlab=" ", col=c("grey","white"),
  legend.text = c("Bona Fide","Cheater"), args.legend = list(x="topleft"),
  beside=TRUE, axes=F)
# bottom panel:
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
pc.line <- array()
for (i in 1:41)
{
  pc.line[i] <- mean((k>=(i-1))==truth)
}

dev.new() # so the plot below does not overwrite the plot above

plot(c(0:40), pc.line, type="l", lwd=2, xlim=c(0,40), ylim=c(0.4,1), 
     xlab="Number of Items Recalled Correctly", 
     ylab=" ", axes=F)
axis(1, at=c(0,seq(from=5,by=5,to=40)))
axis(2, at=c(.5,.75,1))
par(las=0)
mtext("Prop. Correct",side=2, line=2.5,cex=1.5)
# Now add the distribution:
pc.dens <- density(pc)
polygon(c(0,pc.dens$y,0,0), c(pc.dens$x[1]-.01,pc.dens$x,pc.dens$x[1]+.01,pc.dens$x[1]-.01), col="green")


# plot 6.10
par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
    font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
plot(k,samples$BUGSoutput$mean$z,ylim=c(0,1),xlim=c(0,n), xlab= "Number of Items Recalled Correctly", 
     ylab="Cheater Classification", lwd=2, pch=4) 
#in the code, z=0 is bonafide and z=1 is cheating
#so z gives the prob of being assigned to cheating group

