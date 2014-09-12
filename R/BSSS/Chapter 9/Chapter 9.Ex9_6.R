#-----------------------------------------------
# Example 9.5 Bayesian Mixture Factor Analysis
#   
# Programming steps are similar to Example 9.1.  Note
# the discrete categorical distribution "dcat" used to specify a
# categorical distribution with K categories (here 2). Note also
# that the parameter pi of the categorical distribution is given 
# a Dirichlet distribution with hyperparameter alpha.
#-----------------------------------------------

require(rjags) #load the package "rjags"

# Noninformative Prior
# Specify model
modelstring = "								
model {										 		
for (i in 1 : nData) {						
for (j in 1 : 13) {        # NCOL(y) = 13
y[i,j] ~ dnorm(mu[i,j],psi[L[i],j])
}

mu[i,1] <- 1*xi[i,1] 					# Factor 1, lam[1]=1
mu[i,3] <- lam[L[i],3]*xi[i,1]
mu[i,5] <- lam[L[i],5]*xi[i,1]
mu[i,6] <- lam[L[i],6]*xi[i,1]
mu[i,9] <- lam[L[i],9]*xi[i,1]
mu[i,11] <-lam[L[i],11]*xi[i,1]
mu[i,13] <-lam[L[i],13]*xi[i,1]

mu[i,2]  <-1*xi[i,2]				  	# Factor 2, lam[2]=1
mu[i,4] <- lam[L[i],4]*xi[i,2]
mu[i,7] <- lam[L[i],7]*xi[i,2]
mu[i,8] <- lam[L[i],8]*xi[i,2]
mu[i,10] <-lam[L[i],10]*xi[i,2]
mu[i,12] <-lam[L[i],12]*xi[i,2]

L[i]~dcat(pi[1:K])  # Define a categorical distribution

xi[i,1:2] ~ dmnorm(u[1:2],phi[1:2,1:2])
}

# Distributions and Priors

#Priors on Loadings
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

# Priors on Precisions
for(j in 1:13) {
psi[j] ~ dgamma(1, 0.001) # Error variances
sgm[j] <- 1/psi[j]
}
phi[1:2,1:2] ~ dwish(R[1:2,1:2], 5) # Precision matrix

Sigma[1,1] <- phi[2,2]/(phi[1,1]*phi[2,2]-phi[2,1]^2) #Sigma is covariance matrix
Sigma[2,1] <- -phi[1,2]/(phi[1,1]*phi[2,2]-phi[2,1]^2)
Sigma[2,2] <- phi[1,1]/(phi[1,1]*phi[2,2]-phi[2,1]^2)
Sigma[1,2] <- Sigma[2,1]

# Priors on mixture probability
pi[1:K] ~ ddirch(alpha[])
for (j in 1:K) {alpha[j]<-1}
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
K=2;
cfd <- list(y=y, nData=nData,u=u,R=R,K=K)

# RUN CHAIN

# Initialize Model
parameters = c("lam","sgm","Sigma","pi")  #Specify the Parameters to Be Estimated
adaptSteps = 5000
burnInSteps = 5000
nChains = 2
thinSteps = 500
nPerChain = 500000

cfaModel1 = jags.model("model.bug",data=cfd,n.chains=nChains, n.adapt=adaptSteps)
				
# Obtain the Posterior Samples #

cat("Burning in the MCMC chain ...\n")
update(cfaModel1, n.iter=burnInSteps)
cat("Sampling from the final MCMC chain ... \n")
codaSamples1 = coda.samples(cfaModel1, variable.names=parameters,
				n.iter=nPerChain, thin=thinSteps,seed=5555)

summary(codaSamples1)    #Posterior Mean, posterior SD and posterior probablity interval (PPI) 


plot(codaSamples1[[1]][,c(9,1,2,37)])               #Trace plots and Density plots

par(mfrow=c(2,2))
acf(codaSamples1[[1]][,9],main="Loadings of Item 3 in the First Latent Class")   #Auto-correlation plots 
acf(codaSamples1[[1]][,1],main="Variance of Factor TEABEHA")   
acf(codaSamples1[[1]][,2],main="Covariance betw TEABEHA & STUDBEHA")    
acf(codaSamples1[[1]][,37],main="Error Variance of Item 3 in the First Latent Class")    

# Diagnostics
# Geweke Diagnostics and Plot
require(coda)
geweke.diag(codaSamples1[[1]][,c(9,1,2,37)])
geweke.plot(codaSamples1[[1]][,c(9,1,2,37)]) 

# Heidelberger-Welch diagnostics
heidel.diag(codaSamples1[[1]][,c(9,1,2,37)])


# Raftery.diag
raftery.diag(codaSamples1[[1]][,c(9,1,2,37)])
