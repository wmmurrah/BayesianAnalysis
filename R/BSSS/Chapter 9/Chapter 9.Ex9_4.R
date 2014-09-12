#------------------------------------------------------------
# Example 9.4 Bayesian GCM: Baseline model
#
# Programming steps similar to Example 9.3.  Note
# fixing of parameters to yield slope and intercept of growth curve
#------------------------------------------------------------
install.packages("rjags")
require(rjags)   #load the package "rjags"

####Noninformative Prior####
##Specify model
modelstring = "								
model {										 		
for (i in 1 : nData) {						
for (j in 1 : 5) {        # NCOL(y) = 5
y[i,j] ~ dnorm(mu[i,j],psi[j])
}

mu[i,1] <- 1*xi[i,1]+ 0*xi[i,2]	 # xi[i,1] Random Intercept, Initial Status
mu[i,2] <- 1*xi[i,1]+ 1*xi[i,2]         # xi[i,2]: Random Slope, Growth Rate    
mu[i,3] <- 1*xi[i,1]+ 2*xi[i,2]
mu[i,4] <- 1*xi[i,1]+ 3*xi[i,2]
mu[i,5] <- 1*xi[i,1]+ 4*xi[i,2]

xi[i,1:2] ~ dmnorm(v[1:2],phi[1:2,1:2])

}
# Distributions and Priors
#Priors on Precisions
for(j in 1:5) {
psi[j] ~ dgamma(1, 0.001) # Error Precision 
sgm[j] <- 1/psi[j]        # Error variance
}

v[1] ~ dnorm(0, 10^(-6)) #Factor mean, i.e., mean of intercept here    
v[2] ~ dnorm(0, 10^(-6)) #Factor mean, i.e., mean of slope here   

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
y = as.matrix(gcmdata)
gender = gcmdata$gender
nData=NROW(y)

#Specify parameters
R<-matrix(c(1,0,0,1),nrow=2)
gcmdat <- list(y=y, nData=nData, R=R)

# RUN CHAIN

# Initialize Model
parameters = c("v","Sigma","sgm")  #Specify the Parameters to Be Estimated
adaptSteps = 5000
burnInSteps = 5000
nChains = 2
thinSteps = 100
nPerChain = 100000

gcmModel1 = jags.model("model.bug",data=gcmdat,n.chains=nChains, n.adapt=adaptSteps)
				
# Obtain the Posterior Sample of Factor Loadings:

cat("Burning in the MCMC chain ...\n")
update(gcmModel1, n.iter=burnInSteps)
cat("Sampling from the final MCMC chain ... \n")
codaSamples1 = coda.samples(gcmModel1, variable.names=parameters,
				n.iter=nPerChain, thin=thinSteps,seed=5555)

summary(codaSamples1[[1]])     #Posterior Mean, posterior SD and posterior probablity interval (PPI) for the first chain

plot(codaSamples1[[1]][,c(10,11,1,2)])  #Selected Trace plots and Density plots

par(mfrow=c(2,2))
acf(codaSamples1[[1]][,10], main="Fixed Slope of Time")     #Auto-correlation plots  
acf(codaSamples1[[1]][,11], main="Fixed Intercept")     
acf(codaSamples1[[1]][,1], main="Variance of Random Intercept") 
acf(codaSamples1[[1]][,2], main="Covariance of Random Intercept & Random Slope") 

# Diagnostics
#Geweke Diagnostic and Plot

geweke.diag(codaSamples1[[1]])
geweke.plot(codaSamples1[[1]][,c(10,11,1,2)]) 

#Heidelberger-Welch diagnostics
heidel.diag(codaSamples1[[1]])

#Raftery.diag
raftery.diag(codaSamples1[[1]])

\end{spverbatim}
\end{scriptsize}
\newpage

\begin{scriptsize}
\begin{spverbatim}
