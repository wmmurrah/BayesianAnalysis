#----------------------------------------------------------------
# Example 9.2 Path Analysis--Noninformative Prior
#
# Programming steps similar to Example 9.1
#----------------------------------------------------------------

install.packages("MCMCpack")
install.packages("lavaan")
install.packages("rjags")
require(MCMCpack)
require(lavaan)
require(rjags) 


##Specify model
modelstring = "								
model {										 		
#subject-level models
for (i in 1:nData) {
x[i,1:3]~ dmnorm(v[1:3], PI[1:3,1:3])# obtain covariance among X's.

y[i]~dnorm(mu1[i],tau1)
mu1[i]<- a1[1] + a1[2]*momedu[i] + a1[3]*dadedu[i]+a1[4]*perteach[i]+
            a1[5]*importnt[i]+a1[6]*enjoy[i]

enjoy[i]~dnorm(mu2[i], tau2)
mu2[i]<-b1[1]+b1[2]*perteach[i]

importnt[i]~dnorm(mu3[i], tau3)
mu3[i]<-c1[1]+c1[2]*momedu[i]+c1[3]*perteach[i]+c1[4]*enjoy[i]
}##END FOR i

# Prior Specification
#Priors on regression coefficients
v[1:3]  ~ dmnorm(v1[1:3], H0[1:3,1:3])
a1[1:6] ~ dmnorm(u1[1:6], H1[1:6,1:6])  
b1[1:2] ~ dmnorm(u2[1:2], H2[1:2,1:2])
c1[1:4] ~ dmnorm(u3[1:4], H3[1:4,1:4])

#Priors on Precisions
tau1 <- 1/var1
tau2 <- 1/var2
tau3 <- 1/var3
var1 ~ dunif(0, 10^(10))
var2 ~ dunif(0, 10^(10))
var3 ~ dunif(0, 10^(10))
PI[1:3, 1:3] ~ dwish(R[1:3,1:3], 4) 
                                    
}  ##END Model
"
writeLines(modelstring,con="model.bug")

# READ IN DATA AND PREPARE FOR JAGS 

semdata  <- read.csv(file.choose(),header=TRUE) #browse to select data "sem.csv" #colnames(semdata) 
y        <- semdata$mathscor
momedu   <- semdata$momeduc
dadedu   <- semdata$dadeduc
perteach <- semdata$perteach
importnt <- semdata$importnt
enjoy    <- semdata$enjoy
x <- as.matrix(cbind(momedu, dadedu, perteach))#newly added 
nData<-NROW(y)

#Specify variables
Semdata1 <- list(y=y, x=x, nData=nData, v1=rep(0,3),
                u1=rep(0,6),u2=rep(0,2),u3=rep(0,4),H0=diag(10^(-4),3,3),
                H1=diag(c(10^(-10), rep(0.1, 5)),6,6),H2=diag(c(10^(-10), 0.1),2,2), 
                H3=diag(c(10^(-10), rep(0.1, 3)),4,4),R=diag(10^(-4),3,3),
                momedu=momedu, dadedu=dadedu, perteach=perteach, 
                importnt= importnt, enjoy=enjoy)

# RUN CHAIN

# Initialize Model
adaptSteps = 5000
burnInSteps = 5000
nChains = 2
thinSteps = 50
nPerChain = 500000

semModel1 = jags.model("model.bug",data=Semdata1,n.chains=nChains, n.adapt=adaptSteps)

# Obtain the Posterior Sample of Factor Loadings:
parameter=c("a1","b1","c1","tau1","tau2","tau3","v","PI","var1","var2","var3") 
#Specify the Parameters to Be Estimated
cat("Burning in the MCMC chain ...\n")
update(semModel1, n.iter=burnInSteps)
cat("Sampling from the final MCMC chain ... \n")
codaSamples1 = coda.samples(semModel1, n.iter=nPerChain, variable.names=parameter,
                            thin=thinSteps,seed=5555)

summary(codaSamples1)    #Posterior Mean, posterior SD and posterior probablity interval (PPI) 
                         #Note that tau1, tau2 and tau 3 are error precisions instead of variances

plot(codaSamples1[[1]][,c(11,19,17,28)])               #Trace plots and Density plots

par(mfrow=c(2,2))
acf(codaSamples1[[1]][,11],main="MATHSCORE on MOMEDU")             #Auto-correlation plots 
acf(codaSamples1[[1]][,19],main="IMPORTNT on MOMEDU")              
acf(codaSamples1[[1]][,17],main="ENJOY on PERTEACH")               
acf(codaSamples1[[1]][,28],main="Error Variance of MATHSCORE")  

# Diagnostics
#Geweke Diagnostics and Plot
geweke.diag(codaSamples1[[1]][,c(11,19,17,28)])
geweke.plot(codaSamples1[[1]][,c(11,19,17,28)],frac1 = 0.1, frac2 = 0.5) 

#Heidelberg-Welsh Diagnostic
heidel.diag(codaSamples1[[1]][,c(11,19,17,28)])

#Raftery.diag
raftery.diag(codaSamples1[[1]][,c(11,19,17,28)])

#-----------------------------POSTERIOR PREDICTIVE CHECK----------------------------------#

#install.package(ggplot2)
#install.package(gridExtra)
#install.package(sem)
library(ggplot2)
library(gridExtra)
library(MASS)
library(sem)

semdf <- data.frame(momedu = semdata$momeduc, dadedu = semdata$dadeduc, 
  perteach = semdata$perteach, enjoy = semdata$enjoy, importnt = semdata$importnt, y = semdata$mathscor)

sim.y.reps <- function(x) {

#---grab the xth draw from the posterior---#
x <- codaSamples1[[1]][x,1:27]
names(x) <- c(paste("PI",1:9,sep = ""),"mu.y","beta.ym","beta.yd","beta.yp","beta.yi","beta.ye","mu.e",
"beta.ep","mu.i","beta.im","beta.ip","beta.ie",paste("tau",1:3,sep = ""),paste("v",1:3,sep=""))

# Model Implied Covariance Matrix--#
#--inverse of precision matrix obtains X covariance matrix--#
sigma.xx <- solve(matrix(x[1:9],nrow = 3, ncol = 3), byrow = TRUE)

#---Beta Matrix---#
beta.enjoy <- c(0,0,0)
beta.imp <- c(x["beta.ie"],0,0)
beta.math <- c(x["beta.ye"],x["beta.yi"],0)
beta.mat <- matrix(rbind(beta.enjoy,beta.imp,beta.math),nrow = 3, ncol = 3)

#---Gamma Matrix---#
gamma.enjoy <- c(0,0,x["beta.ep"])
gamma.imp <- c(x["beta.im"],0,x["beta.ip"])
gamma.math <- c(x["beta.ym"],x["beta.yd"],x["beta.yp"])
gamma.mat <- matrix(rbind(gamma.enjoy,gamma.imp,gamma.math), nrow = 3, ncol = 3)

#---Structural to Regression conversions---#
err.vect <- c(1/x["tau2"],1/x["tau3"],1/x["tau1"]) ###--Y structural form variances--#
inv.beta <- solve(diag(3) - beta.mat) ###--I minus Beta inverse--#
PI <- inv.beta %*% gamma.mat 
alpha <- inv.beta %*% c(x["mu.e"],x["mu.i"],x["mu.y"]) 
exp.x <- c(x["v1"],x["v2"],x["v3"]) ###--expected value of the X variables--#
exp.y <- alpha + PI %*% exp.x 
sigma.yy <- PI %*% sigma.xx %*% t(PI) + inv.beta %*% diag(err.vect) %*% t(inv.beta) 
sigma.xy <- sigma.xx %*% t(PI) 

#---Need the replicated X variables first---#
x.rep <- mvrnorm(nrow(semdf),exp.x,sigma.xx)

#---Obtain the replicated Y variables next---#
error.mat <- inv.beta %*% 
matrix(
c(rnorm(nrow(semdf),0,sqrt(err.vect[1])),
rnorm(nrow(semdf),0,sqrt(err.vect[2])),
rnorm(nrow(semdf),0,sqrt(err.vect[3]))),
nrow = 3, ncol = nrow(semdf), byrow = TRUE)
y.rep <- matrix(
c(rep(alpha[1],nrow(semdf)),
rep(alpha[2],nrow(semdf)),
rep(alpha[3],nrow(semdf))),
nrow = 3, ncol = nrow(semdf), byrow = TRUE) + 
(PI %*% t(x.rep)) + 
error.mat

#----Use these to double check matrix alebra, as they should give the same thing as above---##
enjoy <- x["mu.e"] + x["beta.ep"]*x.rep[,3] + rnorm(length(semdata$perteach),0,sqrt(1/x["tau2"]))
importnt <- x["mu.i"] + x["beta.im"]*x.rep[,1] + x["beta.ip"]*x.rep[,3] + x["beta.ie"]*enjoy + rnorm(length(semdata$perteach),0,sqrt(1/x["tau3"]))
y <- x["mu.y"] + x["beta.ym"]*x.rep[,1] + x["beta.yd"]*x.rep[,2] + x["beta.yp"]*x.rep[,3] + 
 x["beta.yi"]*importnt + x["beta.ye"]*enjoy + rnorm(length(semdata$perteach),0,sqrt(1/x["tau1"]))

#---Two replicates, one using matrix, one using single line equations (df1 used for double checking)---#
rep.sem.df1 <- data.frame(momedu = x.rep[,1], dadedu = x.rep[,2], perteach = x.rep[,3], enjoy = enjoy, importnt = importnt, y = y)
rep.sem.df <- data.frame(momedu = x.rep[,1], dadedu = x.rep[,2], perteach = x.rep[,3], enjoy = y.rep[1,], importnt = y.rep[2,], y = y.rep[3,])

#---we have three covariance matrices, MODEL, OBSERVED, REPLICATED---#
obs.cov.mat <- cov(semdf)
rep.cov.mat <- cov(rep.sem.df)
mod.cov.mat <- rbind(cbind(sigma.xx,sigma.xy),cbind(t(sigma.xy),sigma.yy))
dimnames(mod.cov.mat) <- dimnames(obs.cov.mat)

#---optional plots to check to make sure the replicates look correct---#
#q1 <- qplot(semdf[,4], main = "enjoy (Observed)"); q2 <- qplot(semdf[,5], main = "important (Observed)"); q3 <- qplot(semdf[,6], main = "math score (Observed)")
#q4 <- qplot(y.rep[1,], main = "enjoy (Replicated)"); q5 <- qplot(y.rep[2,], main = "important (Replicated)"); q6 <- qplot(y.rep[3,], main = "math score (Replicated)")
#pdf(file=' ') #####ADD FILE PATH######
#print(grid.arrange(q1,q2,q3,q4,q5,q6, ncol = 3, nrow = 2))
#dev.off()
#q1 <- qplot(x = rep.sem.df1$y, y = rep.sem.df1$enjoy) + stat_smooth(method = "lm")
#q2 <- qplot(x = rep.sem.df$y, y = rep.sem.df$enjoy) + stat_smooth(method = "lm")
#q3 <- qplot(x = rep.sem.df1$importnt, y = rep.sem.df$importnt)
#print(grid.arrange(q1,q2, ncol = 2))
#q1 <- qplot(semdf[,1], main = "momedu (Observed)"); q2 <- qplot(semdf[,2], main = "dadedu (Observed)"); q3 <- qplot(semdf[,3], main = "perteach (Observed)")
#q4 <- qplot(x.rep[,1], main = "momedu (Replicated)"); q5 <- qplot(x.rep[,2], main = "dadedu (Replicated)"); q6 <- qplot(x.rep[,3], main = "perteach (Replicated)")
#pdf(file=' ') #####ADD FILE PATH######
#print(grid.arrange(q1,q2,q3,q4,q5,q6, ncol = 3, nrow = 2))
#dev.off()

#---Optional for printing out the differing covariances matrices to visually check their accuracy---#
#print("OBS")
#print(obs.cov.mat)
#print("REP")
#print(rep.cov.mat)
#print("MODEL")
#print(mod.cov.mat)

# End Model Implied Covariance Matrix #

#----Calculate the observed and replicated CHI-SQUARE values using the model implied covariance matrix-----#
chi.sq.obs <- 0.5*nrow(semdf)*(determinant(mod.cov.mat)$modulus + 
## note: determinant(x)$modulus returns log of the determinant of x
sum(diag( ##gets the trace
solve(mod.cov.mat) %*% (obs.cov.mat +
(c(exp.x,exp.y) - colMeans(semdf)) %*% t(c(exp.x,exp.y) - colMeans(semdf)) ##mean matrix
)
)) - 
determinant(obs.cov.mat)$modulus - 6 ### p + q = 6?
)##END CHI.SQ calculation
chi.sq.rep <- 0.5*nrow(rep.sem.df)*(determinant(mod.cov.mat)$modulus + ## note: determinant(x)$modulus returns log of the determinant of x
sum(diag( ##gets the trace
solve(mod.cov.mat) %*% (rep.cov.mat +
(c(exp.x,exp.y) - colMeans(rep.sem.df)) %*% t(c(exp.x,exp.y) - colMeans(rep.sem.df)) ##mean matrix
)
)) - 
determinant(rep.cov.mat)$modulus - 6 ### p + q = 6?
)##END CHI.SQ calculation

return(list(rep.sem.df,c(chi.sq.obs,chi.sq.rep)))

}##END function sim.y.reps

##RUN THE POSTERIOR FUNCTIONS
#rand.seed <- round(runif(1,1,10000),0)
#set.seed(rand.seed)
set.seed(515)
y.rep <- lapply(1:1000, sim.y.reps) 
chi.sq <- cbind(unlist(lapply(y.rep,function(x) return(x[[2]][1]))),unlist(lapply(y.rep,function(x) return(x[[2]][2]))))
chi.rep <- chi.sq[,2]
chi.obs <- chi.sq[,1]

##SAVE PDF POSTERIOR CHECKING PLOTS
liks.dif <- chi.obs - chi.rep
range <- max(liks.dif) - min(liks.dif)
p.value <- round(length(which(liks.dif < 0))/length(liks.dif),3)
pdf(file=' ') #####ADD FILE PATH######
par(ask = FALSE)
hist(liks.dif, xlab = expression(Chi["obs"]^2 - Chi["rep"]^2), main = "")
abline(v = 0, lty = 2, lwd = 2)
text(x = min(liks.dif) + range/5, y = 200, label = paste("p-value = ",p.value,sep = ""))
dev.off()
pdf(file=' ') #ADD FILE PATH
qplot(chi.obs,chi.rep, shape = I(1)) + 
geom_abline(slope = 1, intercept = 0) +
theme_bw() +
geom_text(aes(x = -2 + min(chi.obs,chi.rep) , y = - 2 + max(chi.obs,chi.rep), label = paste("p-value =",p.value)), size = 4.5) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
ylim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) + 
xlim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) +
xlab("Observed") + ylab("Replicated")
dev.off()

##DISPLAY POSTERIOR CHECKING PLOTS & INFORMATION
hist(liks.dif, xlab = expression(Chi["obs"]^2 - Chi["rep"]^2), main = "")
abline(v = 0, lty = 2, lwd = 2)
text(x = min(liks.dif) + range/5, y = 200, label = paste("p-value = ",p.value,sep = ""))
par(ask = TRUE)
qplot(chi.obs,chi.rep, shape = I(1)) + 
geom_abline(slope = 1, intercept = 0) +
theme_bw() +
geom_text(aes(x = -2 + min(chi.obs,chi.rep) , y = - 2 + max(chi.obs,chi.rep), label = paste("p-value =",p.value)), size = 4.5) +
#theme(panel.grid.major = element_blank()) +
ylim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) + 
xlim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) +
xlab("Observed") + ylab("Replicated")


