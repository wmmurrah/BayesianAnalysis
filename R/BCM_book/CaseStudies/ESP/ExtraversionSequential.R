# Sequential analysis, not discussed in the book
# clears workspace:  
rm(list=ls()) 

# sets working directories:
setwd("C:/Users/EJ/Dropbox/EJ/temp/BayesBook/test/CaseStudies/ESP")

# Extraversion and Performance on erotic pictures in session 1.
extrav   <- c(50, 80, 79, 56, 50, 80, 53, 84, 74, 67,
             50, 45, 62, 65, 71, 71, 68, 63, 67, 58,
             72, 73, 63, 54, 63, 70, 81, 71, 66, 74, 
             70, 84, 66, 73, 78, 64, 54, 74, 62, 71,
             70, 79, 66, 64, 62, 63, 60, 56, 72, 72,
             79, 67, 46, 67, 77, 55, 63, 44, 84, 65,
             41, 62, 64, 51, 46, 53, 26, 67, 73, 39,
             62, 59, 75, 65, 60, 69, 63, 69, 55, 63,
             86, 70, 67, 54, 80, 71, 71, 55, 57, 41,
             56, 78, 58, 76, 54, 50, 61, 60, 32, 67)
  
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

BF10.J.exact <- function(n, r)
{
  integrand <- function(rho) {((1-rho^2)^(n/2)) / ((1-rho*r)^(n-.5))}
  BF10      <- integrate(integrand, lower=-1, upper=1)$value/2
  return(BF10)
}

BF10 <- array(dim=100)
BF10[1] <- 1
for (i in 2:100)
{
  BF10[i] <- BF10.J.exact(n=i-1, r=cor(extrav[1:i],prc1.ero[1:i]))
} 

#============ log Bayes factors  ===========================
# with thanks to Ruud Wetzels and Benjamin Scheibehenne
windows(8,6)
par(cex.main = 1.3, mar = c(4.5, 6, 4, 7)+.1, mgp = c(3, 1, 0),   #bottom, left, top, right
  cex.lab = 1.3, font.lab = 2, cex.axis = 1.3, las=1)

plot(log(BF10), xlim=c(1,100), ylim=c(-1*log(200),log(200)), xlab="No. of Participants",
           ylab="", cex.lab=1.3,cex.axis=1.3, las =1, yaxt="n", bty = "n",
           type="p", pch=21, bg="grey")

labelsUpper <- log(c(100,30,10,3,1))
labelsLower <- -1*labelsUpper
criticalP <- c(labelsLower,0,labelsUpper)
for (idx in 1:length(criticalP)) 
{
  abline(h=criticalP[idx],col='darkgrey',lwd=1,lty=2)
}
abline(h=0)
axis(side=4, at=criticalP,tick=T,las=2,cex.axis=1, labels=F)
axis(side=4, at=labelsUpper+.602, tick=F, cex.axis=1, labels=c("Extreme","Very strong", "Strong","Substantial", "Anecdotal"))
axis(side=4, at=labelsLower-.602,tick=F, cex.axis=1, labels=c("Extreme","Very strong", "Strong","Substantial", "Anecdotal"))

axis(side=2, at=c(criticalP),tick=T,las=2,cex.axis=1,
labels=c("-log(100)","-log(30)","-log(10)","-log(3)","log(1)","", "log(100)","log(30)","log(10)","log(3)",""))
 
mtext("log BF10", side=2, line=4.5, las=0, cex=1.3)
mtext("Evidence", side=4, line=5.5, las=0, cex=1.3)

arrows(20, -log(10), 20, -log(100), length=.25, angle=30, code=2, lwd=2)
arrows(20, log(10), 20, log(100), length=.25, angle=30, code=2, lwd=2)
text(25, -log(70), "Evidence for H0", pos=4, cex=1.3)
text(25, log(70), "Evidence for H1", pos=4, cex=1.3)
###########################################################################
