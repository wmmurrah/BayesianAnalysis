install.packages("pastecs")
install.packages("coda")
require(pastecs)
require(coda)

#-----------------------------------------------------------------------------------------
# EXAMPLE 4.1:  Metropolis-Hastings Algorithm with a normal proposal distribution 
#                             Program Steps:
# 1.  Define a numeric vector of length n and set the first value to zero 
# 2.  Draw one value from a proposal normal distribution with mean 0 and variance 1
# 3.  Define a candidate value as the value of the target plus the first value 
# 4.  Define the acceptance probability "aprob" as the minimum value of 1 or the 
#      ratio of the density values of candidate value and x.
# 5.  Decide to accept or reject the candidate value by comparing ``aprob" to a random  
#      draw from a uniform(0, 1) distribution
# 6.  Summarize and plot the results   
#-----------------------------------------------------------------------------------------#

MHNorm <- function (n,burnin,print=FALSE) 
{
        vec <- vector("numeric", n)
        x <- 0
        vec[1] <- x
        for (i in 2:n) {
                target <- rnorm(1,0,1)
                can <- x + target
                aprob <- min(1, dnorm(can)/dnorm(x))
                u <- runif(1)
                if (u < aprob) 
                        x <- can
                vec[i] <- x
        }
      		  summary <- stat.desc(vec[burnin:n], basic=F)
            print(summary)
        	       	
        	  if (print==TRUE) print(vec[burnin:n])
	        	
par(mfrow=c(2,1))
plot(ts(vec[burnin:n]))
hist(vec[burnin:n],30,main=NULL,xlab=NULL)
par(mfrow=c(1,1))
}
# End

#--Metropolis-Hastings Algorithm with a Uniform Proposal Distribution--#
require(pastecs)
require(coda)

MHUnif <- function (n,burnin,tau,print=FALSE) 
{
        vec <- vector("numeric",n)
        x <- 0
        vec[1] <- x
        for (i in 2:n) {
                target <- runif(1, -tau, tau)
                can <- x + target
                aprob <- min(1, dnorm(can)/dnorm(x))
                u <- runif(1)
                if (u < aprob) 
                        x <- can
                			  vec[i] <- x         
        }
            summary <- stat.desc(vec[burnin:n], basic=F)
            print(summary)
        	
        	    if (print==TRUE) print(vec[burnin:n])        		
        	
par(mfrow=c(2,1))
plot(ts(vec[burnin:n]))
hist(vec[burnin:n],30)
par(mfrow=c(1,1))
}
# End

#--------------------------------------------------------------------------
# EXAMPLE 4.2:  Gibbs Sampler
#         Program Steps
# 1. Define a sample of size N and generate data from a normal distribution with
       chosen mean and variance (here, 0 and 5, respectively)
# 2. Create a table to be filled in later
# 3. Write JAGS code beginning with "modelstring="
# 4. Define a probability distribution for x and for hyperparameters mu and tau
# 5. End JAGS code and return to rjags
# 6. Create a .bug file that contains the model string
# 7. Define rjags parameters
# 8. Run jags.model which reads in model.bug and model parameters
# 9. Summarize with coda and diagnostics
#--------------------------------------------------------------------------


install.packages("rjags")  # requires that "jags" be installed"
install.packages("coda")
require(rjags)
require(coda)

N <- 1000
x <- rnorm(N, 0, 5)

write.table(x,
            file = 'example1.data',
            row.names = FALSE,
            col.names = FALSE)

#--------------------------------#
# JAGS code starts here     
#--------------------------------#

modelstring="
model {
for (i in 1:N) {
x[i] ~ dnorm(mu, tau)
}
mu ~ dnorm(0, .0001)
tau ~ dgamma(0.001,0.001)
}
"
#----------------------------------#
# JAGS code ends here    
#---------------------------------#
writeLines(modelstring,con="model.bug")
parameters = c("mu","tau")  #Specify the Parameters to Be Estimated
adaptSteps = 500
burnInSteps = 1000
nChains = 2
thinSteps = 10
nPerChain = 100000

foo <- jags.model("model.bug",
                   data = list('x' = x,
                               'N' = N),
                   n.chains = nChains,
                   n.adapt = adaptSteps)

cat("Burning in the MCMC chain ...\n")
update(foo, n.iter=burnInSteps)
cat("Sampling from the final MCMC chain ... \n")

codaSamples1 = coda.samples(foo, variable.names=parameters,
                            n.iter=nPerChain, thin=thinSteps,seed=2847)

summary(codaSamples1[[1]])     #Posterior Mean, posterior SD and posterior probablity interval (PPI) for the first chain

plot(codaSamples1,trace=F)
plot(codaSamples1,density=F,col="black")

par(mfrow=c(2,2))
autocorr.plot(codaSamples1[[1]])
autocorr.plot(codaSamples1[[2]])

geweke.plot(codaSamples1)
geweke.diag(codaSamples1)

#Gelman Plot
gelman.diag(codaSamples1)
gelman.plot(codaSamples1)

#Heidelberger-Welch diagnostics
heidel.diag(codaSamples1[[1]])

#Raftery.diag
raftery.diag(codaSamples1[[1]])

# End