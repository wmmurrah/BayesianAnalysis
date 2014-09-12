#------------------------------------------------------------
# Example 9.4 Bayesian GCM: With predictor
#
# Programming steps similar to Example 9.3.  Note
# models for random intercept and random slope as functions of
# gender
#------------------------------------------------------------

require(rjags)   #load the package "rjags"

####Noninformative Prior####
##Specify model
modelstring = "								
model {										 		
for (i in 1 : nData) {						
for (j in 1 : 5) {        # NCOL(y) = 5
y[i,j] ~ dnorm(mu[i,j],psi[j])
}

mu[i,1] <- 1*xi[i,1]+ 0*xi[i,2]	 # xi[i,1]: Random Intercept, Initial Status
mu[i,2] <- 1*xi[i,1]+ 1*xi[i,2]         # xi[i,2]: Random Slope, Growth Rate    
mu[i,3] <- 1*xi[i,1]+ 2*xi[i,2]
mu[i,4] <- 1*xi[i,1]+ 3*xi[i,2]
mu[i,5] <- 1*xi[i,1]+ 4*xi[i,2]

xi[i,1] <- alpha1+ beta1*gender[i]+ eta[i,1] #Model for random intercept
xi[i,2] <- alpha2+ beta2*gender[i]+ eta[i,2] # Model for random slope

eta[i,1:2] ~ dmnorm(v[1:2],phi[1:2,1:2])
}

# Distributions and Priors

#Priors on Precisions
for(j in 1:5) {
psi[j] ~ dgamma(1, 0.001) # Error Precision 
sgm[j] <- 1/psi[j]        # Error variance
}

alpha1 ~ dnorm(0, 10^(-6)) #Mean of intercept    
alpha2 ~ dnorm(0, 10^(-6)) #Mean of slope 
beta1 ~ dnorm(0, 10^(-6))  #Intercept on gender
beta2 ~ dnorm(0, 10^(-6))  #Slope on gender    

phi[1:2,1:2] ~ dwish(R[1:2,1:2], 5)                   #phi is precision matrix of factors
Sigma[1,1] <- phi[2,2]/(phi[1,1]*phi[2,2]-phi[2,1]^2) #Sigma is covariance matrix of factors
Sigma[2,1] <- -phi[1,2]/(phi[1,1]*phi[2,2]-phi[2,1]^2)
Sigma[2,2] <- phi[1,1]/(phi[1,1]*phi[2,2]-phi[2,1]^2)
Sigma[1,2] <- Sigma[2,1]

}  
"
writeLines(modelstring,con="model.bug")

# READ IN DATA AND PREPARE FOR JAGS 

gcmdata = read.csv(file.choose(),header=TRUE) #browse to select data "lsay_ex.csv".
y = as.matrix(gcmdata[,2:6])
gender = gcmdata$gender
nData=NROW(y)

#Specify parameters
v<-c(0,0)
R<-matrix(c(1,0,0,1),nrow=2)
gcmdat <- list(y=y, gender=gender, nData=nData, v=v, R=R)

# RUN CHAIN

# Initialize Model
parameters = c("alpha1", "alpha2","beta1","beta2","Sigma","sgm")  #Specify the Parameters to Be Estimated
adaptSteps = 5000
burnInSteps = 5000
nChains = 2
thinSteps = 100
nPerChain = 100000

gcmModel2 = jags.model("model.bug",data=gcmdat,n.chains=nChains, n.adapt=adaptSteps)
				
# Obtain the Posterior Sample #

cat("Burning in the MCMC chain ...\n")
update(gcmModel2, n.iter=burnInSteps)
cat("Sampling from the final MCMC chain ... \n")
codaSamples2 = coda.samples(gcmModel2, variable.names=parameters,
				n.iter=nPerChain, thin=thinSteps,seed=5555)

summary(codaSamples2[[1]])     #Posterior Mean, posterior SD and posterior probablity interval (PPI) for the first chain

plot(codaSamples2[[1]][,c(7,8,1,4)])  #Selected Trace plots and Density plots

par(mfrow=c(2,2))
acf(codaSamples2[[1]][,7],main="Intercept on Gender")     #Auto-correlation plots  
acf(codaSamples2[[1]][,8],main="Slope on Gender")     
acf(codaSamples2[[1]][,1],main="Variance of Random Intercept") 
acf(codaSamples2[[1]][,4],main="Variance of Random Slope") 

# Diagnostics
#Geweke Diagnostic and Plot
require(coda)

geweke.diag(codaSamples2[[1]])
geweke.plot(codaSamples2[[1]][,c(7,8,1,4)]) 

#Heidelberger-Welch diagnostics
heidel.diag(codaSamples2[[1]])

#Raftery.diag
raftery.diag(codaSamples2[[1]])

