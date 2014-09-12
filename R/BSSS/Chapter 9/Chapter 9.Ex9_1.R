#-------------------------------------------------------
# EXAMPLE 9.1 Bayesian CFA with non-informative and
# informative priors
#         Program Steps
# 1. Specify model.  Note that one loading in each column of the loading matrix
#     must be set to 1.0.
# 2. Place priors on all model parameters except the fixed factor loading
# 3. End modelstring and write lines to a model.bug connection
# 4. Read in data and specify parameters and initialization values for the chain
# 5. Summarize
#
#  Remaining code obtains posterior predictive plots
#-------------------------------------------------------

install.packages("rjags")
install.packages("MCMCpack")
require(rjags)    
require(MCMCpack) 

 
# Specify model
modelstring = "								
model {										 		
for (i in 1 : nData) {						
for (j in 1 : 13) {        # NCOL(y) = 13
y[i,j] ~ dnorm(mu[i,j],psi[j])
}
mu[i,1] <- w[1]+1*xi[i,1] 					# Factor 1, lam[1]=1
mu[i,3] <- w[3]+lam[3]*xi[i,1]
mu[i,5] <- w[5]+lam[5]*xi[i,1]
mu[i,6] <- w[6]+lam[6]*xi[i,1]
mu[i,9] <- w[9]+lam[9]*xi[i,1]
mu[i,11] <-w[11]+lam[11]*xi[i,1]
mu[i,13] <-w[13]+lam[13]*xi[i,1]

mu[i,2]  <-w[2]+1*xi[i,2]				  	# Factor 2, lam[2]=1
mu[i,4] <- w[4]+lam[4]*xi[i,2]
mu[i,7] <- w[7]+lam[7]*xi[i,2]
mu[i,8] <- w[8]+lam[8]*xi[i,2]
mu[i,10] <-w[10]+lam[10]*xi[i,2]
mu[i,12] <-w[12]+lam[12]*xi[i,2]

xi[i,1:2] ~ dmnorm(u[1:2],phi[1:2,1:2])
}

# Distributions and Priors
# Priors on Intercepts and Loadings
lam[1]<-1
lam[3] ~ dnorm(0, 10^(-6)) #prior mean and prior precision
lam[5] ~ dnorm(0, 10^(-6))
lam[6] ~ dnorm(0, 10^(-6))
lam[9] ~ dnorm(0, 10^(-6))
lam[11] ~ dnorm(0, 10^(-6))
lam[13] ~ dnorm(0, 10^(-6))

lam[2]<-1 
lam[4] ~ dnorm(0, 10^(-6))
lam[7] ~ dnorm(0, 10^(-6))
lam[8] ~ dnorm(0, 10^(-6))
lam[10] ~ dnorm(0, 10^(-6))
lam[12] ~ dnorm(0, 10^(-6))

#prior mean and prior precision
w[1] ~ dnorm(0, 10^(-6)) 
w[3] ~ dnorm(0, 10^(-6)) 
w[5] ~ dnorm(0, 10^(-6))
w[6] ~ dnorm(0, 10^(-6))
w[9] ~ dnorm(0, 10^(-6))
w[11] ~ dnorm(0, 10^(-6))
w[13] ~ dnorm(0, 10^(-6))

w[2] ~ dnorm(0, 10^(-6)) 
w[4] ~ dnorm(0, 10^(-6))
w[7] ~ dnorm(0, 10^(-6))
w[8] ~ dnorm(0, 10^(-6))
w[10] ~ dnorm(0, 10^(-6))
w[12] ~ dnorm(0, 10^(-6))

#---------------------------------------
#Priors on Precisions
for(j in 1:13) {
psi[j] ~ dgamma(1, 0.001) # Precision of residuals
sgm[j] <- 1/psi[j]        # Error variance
}
          
phi[1:2,1:2] ~ dwish(R[1:2,1:2], 3)      #phi is precision matrix

Sigma[1,1] <- phi[2,2]/(phi[1,1]*phi[2,2]-phi[2,1]^2) #Sigma is covariance matrix
Sigma[2,1] <- -phi[1,2]/(phi[1,1]*phi[2,2]-phi[2,1]^2)
Sigma[2,2] <- phi[1,1]/(phi[1,1]*phi[2,2]-phi[2,1]^2)
Sigma[1,2] <- Sigma[2,1]
}  
"
writeLines(modelstring,con="model.bug")

# READ IN DATA AND PREPARE FOR JAGS 
cfadata = read.csv(file.choose(),header=TRUE) #browse to select data "cfa_imputed.csv" (using regression imputation).
y = as.matrix(cfadata)
nData=NROW(y)

#Specify parameters
u=c(0,0)
R=matrix(c(1,0,0,1),nrow=2)
cfadata <- list(y=y, nData=nData,u=u,R=R)

# RUN CHAIN
# Initialize Model
parameters = c("lam","Sigma","sgm","phi","xi","w")  #Specify the Parameters to Be Estimated
adaptSteps = 5000
burnInSteps = 5000
nChains = 2
thinSteps = 500
nPerChain = 100000

cfaModel1 = jags.model("model.bug",data=cfadata,n.chains=nChains, n.adapt=adaptSteps)
				
# Obtain the Posterior Sample of Factor Loadings: #
cat("Burning in the MCMC chain ...\n")
update(cfaModel1, n.iter=burnInSteps)
cat("Sampling from the final MCMC chain ... \n")
codaSamples1 = coda.samples(cfaModel1, variable.names=parameters,
				n.iter=nPerChain, thin=thinSteps,seed=5555)

#Posterior Mean, posterior SD and posterior probablity interval (PPI) for the first chain #
summary(codaSamples1[[1]])     
#posterior mean of factor correlation #
mean(codaSamples1[[1]][,2]/sqrt(codaSamples1[[1]][,1]*codaSamples1[[1]][,4]))    
#posterior SD of factor correlation #
sd(codaSamples1[[1]][,2]/sqrt(codaSamples1[[1]][,1]*codaSamples1[[1]][,4]))      

# 95% PPI of factor corr.#
quantile(codaSamples1[[1]][,2]/sqrt(codaSamples1[[1]][,1]*codaSamples1[[1]][,4]),c(.025,.975))
#Selected Trace plots and Density plots for loadings of items 8, 9, 10 and 12 (1st chain). #
plot(codaSamples1[[1]][,c(12:14,16)])  
#Selected Trace plots and Density plots for variance terms #
plot(codaSamples1[[1]][,c(1,2,4,22)])  

#Auto-correlation plots for selected items (first chain) #
par(mfrow=c(2,2))
acf(codaSamples1[[1]][,12],main="Item 8")    
acf(codaSamples1[[1]][,13],main="Item 9")      
acf(codaSamples1[[1]][,14],main="Item 10")     
acf(codaSamples1[[1]][,16],main="Item 12")     

par(mfrow=c(2,2))
acf(codaSamples1[[1]][,1],main="Variance of Factor TEABEHA")     
acf(codaSamples1[[1]][,2],main="Covariance betw TEABEHA and STUDBEHA")     
acf(codaSamples1[[1]][,4],main="Variance of Factor STUDBEHA")     
acf(codaSamples1[[1]][,22],main="Error Variance of Item 1")   

# Diagnostics
#Geweke Diagnostic and Plot
geweke.diag(codaSamples1[[1]][,c(14,1,2,22)])
geweke.plot(codaSamples1[[1]][,c(14,1,2,22)]) 

#Heidelberger-Welch diagnostics
heidel.diag(codaSamples1[[1]][,c(14,1,2,22)])
                                    
#Raftery.diag
raftery.diag(codaSamples1[[1]][,c(14,1,2,22)])

#------------------------POSTERIOR PREDICTIVE CHECK ----------------------------------#
    NEEDS TO BE FIXED
####


#install.package(ggplot2)
#install.package(gridExtra)
#install.package(lavaan)
library(ggplot2)
library(gridExtra)
library(lavaan)

##CREATE THE REPLICATED DATA
draw.reps <- function(x) {
(x <- codaSamples1[[1]][x,])
(f.draw <- mvrnorm(n = 165, mu = c(0,0), Sigma = matrix(c(x[1],x[2],x[3],x[4]),nrow = 2, ncol = 2)))
#print(f.draw)

## create the replicated data set
y.rep <- matrix(nrow = 165, ncol = 13)
for(i in 1:13) {
if(is.na(match(i,c(1,3,5,6,9,11,13)))) { fact <- 2 } else { fact <- 1 } 
y.rep[,i] <- x[i+34] + x[i+4]*f.draw[,fact] + rnorm(165, 0, sqrt(x[i+21]))
}##END for

## create the model implied covariance mat
obs.cov.mat <- cov(cfadata$y)
rep.cov.mat <- cov(y.rep)
lambda <- matrix(rep(0,26),nrow = 13, ncol = 2)
lambda[c(1,3,5,6,9,11,13),1] <- x[4 + c(1,3,5,6,9,11,13)]
lambda[c(2,4,7,8,10,12),2] <- x[4 + c(2,4,7,8,10,12)]
F.cov <- matrix(c(x[1],x[2],x[3],x[4]),nrow = 2, ncol = 2)
phi <- diag(x[22:34])
exp.cov.mat <- lambda %*% F.cov %*% t(lambda) + phi

## Calculate the observed and replicated X^2 values using the model implied covariance matrix
chi.sq.obs <- nrow(cfadata$y)*(determinant(exp.cov.mat)$modulus + ## note: determinant(x)$modulus returns log of the determinant of x
sum(diag( ##gets the trace
(obs.cov.mat) %*% solve(exp.cov.mat)
)) - 
determinant(obs.cov.mat)$modulus - 13 ### p + q = 13
)##END CHI.SQ calculation
chi.sq.rep <- nrow(y.rep)*(determinant(exp.cov.mat)$modulus + ## note: determinant(x)$modulus returns log of the determinant of x
sum(diag( ##gets the trace
(rep.cov.mat) %*% solve(exp.cov.mat)
)) - 
determinant(rep.cov.mat)$modulus - 13) 

return(list(y.rep,c(chi.sq.obs,chi.sq.rep)))
}##END function draw.reps


##DISCREP STAT USING LAVAAN
discrep.stat.lav <- function(x) {
lav.path <- ' chi1 =~ SC17Q01 + SC17Q03 + SC17Q05 + SC17Q06 + SC17Q09 + SC17Q11 + SC17Q13
chi2 = ~ SC17Q02 + SC17Q04 + SC17Q07 + SC17Q08 + SC17Q10 + SC17Q12
'
x <- data.frame(x[[1]])
names(x) <- dimnames(cfadata$y)[[2]]
rep.fit <- cfa(lav.path, data = x)
return(rep.fit@Fit@test[[1]]$stat)
}##END function discrep.stat.lav

##EASY FUNCTION TO MAKE FUNCTION CALL
get.yreps <- function(z) {
y.reps <- lapply(1:z,draw.reps)
return(y.reps)
}##END function

##BEGIN POSTERIOR CHECK
set.seed(515)
y.reps <- get.yreps(1000)
chi.sq <- cbind(unlist(lapply(y.reps,function(x) return(x[[2]][1]))),unlist(lapply(y.reps,function(x) return(x[[2]][2]))))
chi.rep <- chi.sq[,2]
chi.obs <- chi.sq[,1]

##ChiSq Using lavaan 
#lav.path <- ' chi1 =~ SC17Q01 + SC17Q03 + SC17Q05 + SC17Q06 + SC17Q09 + SC17Q11 + SC17Q13
#chi2 = ~ SC17Q02 + SC17Q04 + SC17Q07 + SC17Q08 + SC17Q10 + SC17Q12
#'
#fit <- cfa(lav.path, data = data.frame(cfadata$y))
#chi.obs.lav <- fit@Fit@test[[1]]$stat
#chi.rep.lav <- unlist(lapply(y.reps,discrep.stat.lav)) ##grabs discrep stat from y.reps, will take a minute

##SAVE PDF POSTERIOR CHECKING PLOTS
liks.dif <- chi.obs - chi.rep
range <- max(liks.dif) - min(liks.dif)
p.value <- round(length(which(liks.dif < 0))/length(liks.dif),3)
pdf(file=' ') #####ADD FILE PATH######
par(ask = FALSE)
hist(liks.dif, xlab = expression(Chi["obs"]^2 - Chi["rep"]^2), main = "", xlim = c(-60,150))
abline(v = 0, lty = 2, lwd = 2)
text(x = -35, y = 200, label = paste("p-value = ",p.value,sep = ""))
dev.off()
pdf(file=' ') #####ADD MY FILE PATH######
qplot(chi.obs,chi.rep, shape = I(1)) + 
geom_abline(slope = 1, intercept = 0) +
theme_bw() +
geom_text(aes(x = 7 + min(chi.obs,chi.rep) , y = - 2 + max(chi.obs,chi.rep), label = paste("p-value =",p.value)), size = 4.5) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
ylim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) + 
xlim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) +
xlab("Observed") + ylab("Replicated")
dev.off()


##DISPLAY POSTERIOR CHECKING PLOTS & INFORMATION
hist(liks.dif, xlab = expression(Chi["obs"]^2 - Chi["rep"]^2), main = "", xlim = c(-60,150))
abline(v = 0, lty = 2, lwd = 2)
text(x = -35, y = 200, label = paste("p-value = ",p.value,sep = ""))
par(ask = TRUE)
qplot(chi.obs,chi.rep, shape = I(1)) + 
geom_abline(slope = 1, intercept = 0) +
theme_bw() +
geom_text(aes(x = 7 + min(chi.obs,chi.rep) , y = - 2 + max(chi.obs,chi.rep), label = paste("p-value =",p.value)), size = 4.5) +
#theme(panel.grid.major = element_blank()) +
ylim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) + 
xlim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) +
xlab("Observed") + ylab("Replicated")


####Informative Prior####
##Specify model
require(rjags)

modelstring = "								
model {										 		
for (i in 1 : nData) {						
for (j in 1 : 13) {        # NCOL(y) = 13
y[i,j] ~ dnorm(mu[i,j],psi[j])
}

mu[i,1] <- w[1]+1*xi[i,1] 					# Factor 1, lam[1]=1
mu[i,3] <- w[3]+lam[3]*xi[i,1]
mu[i,5] <- w[5]+lam[5]*xi[i,1]
mu[i,6] <- w[6]+lam[6]*xi[i,1]
mu[i,9] <- w[9]+lam[9]*xi[i,1]
mu[i,11] <-w[11]+lam[11]*xi[i,1]
mu[i,13] <-w[13]+lam[13]*xi[i,1]

mu[i,2]  <-w[2]+1*xi[i,2]				  	# Factor 2, lam[2]=1
mu[i,4] <- w[4]+lam[4]*xi[i,2]
mu[i,7] <- w[7]+lam[7]*xi[i,2]
mu[i,8] <- w[8]+lam[8]*xi[i,2]
mu[i,10] <-w[10]+lam[10]*xi[i,2]
mu[i,12] <-w[12]+lam[12]*xi[i,2]

xi[i,1:2] ~ dmnorm(u[1:2],phi[1:2,1:2])
}

# Distributions and Priors
#Priors on Loadings and Intercepts
lam[1]<-1
lam[3] ~ dnorm(0.729, 1/.018) #prior mean and prior precision 
lam[5] ~ dnorm(1.176, 1/0.028)
lam[6] ~ dnorm(0.887, 1/0.025)
lam[9] ~ dnorm(1.013, 1/0.029)
lam[11] ~ dnorm(0.645, 1/0.018)
lam[13] ~ dnorm(1.121, 1/0.031)

lam[2]<-1 
lam[4] ~ dnorm(0.810, 1/0.026)
lam[7] ~ dnorm(1.359, 1/0.055)
lam[8] ~ dnorm(0.931, 1/0.031)
lam[10] ~ dnorm(0.970, 1/0.035)
lam[12] ~ dnorm(0.942, 1/0.030)

w[1] ~ dnorm(0, 10^(-6)) 
w[3] ~ dnorm(0, 10^(-6)) #prior mean and prior precision
w[5] ~ dnorm(0, 10^(-6))
w[6] ~ dnorm(0, 10^(-6))
w[9] ~ dnorm(0, 10^(-6))
w[11] ~ dnorm(0, 10^(-6))
w[13] ~ dnorm(0, 10^(-6))

w[2] ~ dnorm(0, 10^(-6)) 
w[4] ~ dnorm(0, 10^(-6))
w[7] ~ dnorm(0, 10^(-6))
w[8] ~ dnorm(0, 10^(-6))
w[10] ~ dnorm(0, 10^(-6))
w[12] ~ dnorm(0, 10^(-6))

#---------------------------------------
#Priors on Precisions
for(j in 1:13) {
psi[j] ~ dgamma(1, 0.001) # Precision of residuals
sgm[j] <- 1/psi[j]        # Error variance
}

phi[1:2,1:2] ~ dwish(R[1:2,1:2], 3)                   #phi is precision matrix
Sigma[1,1] <- phi[2,2]/(phi[1,1]*phi[2,2]-phi[2,1]^2) #Sigma is covariance matrix
Sigma[2,1] <- -phi[1,2]/(phi[1,1]*phi[2,2]-phi[2,1]^2)
Sigma[2,2] <- phi[1,1]/(phi[1,1]*phi[2,2]-phi[2,1]^2)
Sigma[1,2] <- Sigma[2,1]

}  
"
writeLines(modelstring,con="model.bug")

# READ IN DATA AND PREPARE FOR JAGS 

cfadata = read.csv(file.choose(),header=TRUE) 
y = as.matrix(cfadata)
nData=NROW(y)

#Specify parameters
u=c(0,0)
R=matrix(c(1,0,0,1),nrow=2)
cfadata <- list(y=y, nData=nData,u=u,R=R)

# RUN CHAIN

# Initialize Model
parameters = c("lam","Sigma","sgm","phi","xi","w")  #Specify the Parameters to Be Estimated
adaptSteps = 5000
burnInSteps = 5000
nChains = 2
thinSteps = 500
nPerChain = 1000000

cfaModel2 = jags.model("model.bug",data=cfadata,n.chains=nChains, n.adapt=adaptSteps)
				
# Obtain the Posterior Sample of Factor Loadings:
cat("Burning in the MCMC chain ...\n")
update(cfaModel2, n.iter=burnInSteps)
cat("Sampling from the final MCMC chain ... \n")
codaSamples2 = coda.samples(cfaModel2, variable.names=parameters,
				n.iter=nPerChain, thin=thinSteps,seed=5555)
#Posterior Mean, posterior SD and posterior probablity interval (PPI) for the 1st chain
summary(codaSamples2[[1]])  
#Factor correlation  
mean(codaSamples2[[1]][,2]/sqrt(codaSamples2[[1]][,1]*codaSamples2[[1]][,4]))   
#Posterior SD of factor correlation
sd(codaSamples2[[1]][,2]/sqrt(codaSamples2[[1]][,1]*codaSamples2[[1]][,4]))       
quantile(codaSamples2[[1]][,2]/sqrt(codaSamples2[[1]][,1]*codaSamples2[[1]][,4]),c(.025,.975))#95% PPI of factor corr.

# Trace and Density Plots for Selected Loadings and Variance Terms (first chain) #
plot(codaSamples2[[1]][,c(12:14,16)])  
plot(codaSamples2[[1]][,c(1,2,4,22)])   

#Selected Auto-Correlation Plots (first chain) #
par(mfrow=c(2,2))
acf(codaSamples2[[1]][,12],main="Item 8")    
acf(codaSamples2[[1]][,13],main="Item 9")    
acf(codaSamples2[[1]][,14],main="Item 10")     
acf(codaSamples2[[1]][,16],main="Item 12")    

par(mfrow=c(2,2))
acf(codaSamples2[[1]][,1],main="Variance of Factor TEABEHA")    
acf(codaSamples2[[1]][,2],main="Covariance betw TEABEHA and STUDBEHA")     
acf(codaSamples2[[1]][,4],main="Variance of Factor STUDBEHA")     
acf(codaSamples2[[1]][,22],main="Error Variance of Item 1")   

# Diagnostics
#Geweke Diagnostics and Plot
geweke.diag(codaSamples2[[1]][,c(14,1,2,22)])
geweke.plot(codaSamples2[[1]][,c(14,1,2,22)],frac1 = 0.1, frac2 = 0.5) 

#Heidelberger-Welch diagnostics
heidel.diag(codaSamples2[[1]][,c(14,1,2,22)])
  
#Raftery.diag
raftery.diag(codaSamples2[[1]][,c(14,1,2,22)])


