
#---------------------------------------------------------------------------------------#
#  The following R code generates artificial data and plots the data distribution, 
#  the prior distribution under varying degrees of precision of the prior,
#  and the resulting posterior distribution.  
#
#               Program steps common for each plot
# 1. Create a sequence of values for the x-axis and artifical data for the y-axis.
# 2. Call a variety of available density functions in R.  
# 3. Define the posterior as the product of the prior and the likelihood.
# 4. Generate plots.
#---------------------------------------------------------------------------------------#

install.packages("scatterplot3d")
install.packages("pscl")
require(scatterplot3d)
require(pscl)

#-----------------------------------------------------------------------------#
#  Example 3.1.  NORMAL DISTRIBUTION: MEAN KNOWN, VARIANCE UNKNOWN 
#  PRIOR FOR SIGMA^2 IS INVERSE-GAMMA
#------------------------------------------------------------------------------#

x <- seq(0.1,10,by=.001)  
y <- c(-1,-0.7,-0.5,-0.2,0.1,0.3,0.6,1.2)  
n=8
par(mfrow=c(2,2),lwd=2,mar=c(3,3,4,4),cex.axis=.6)

# INVERSE-GAMMA(5,3) -- HIGH PRECISION

prior = densigamma(x,5,3)  
like = densigamma(x,sum(y^2),n/2-1)  
and known mean parameter 
post_propor = prior*like 
#Approx. posterior density; sum(post)*0.001 = 1
post = post_propor / (0.001*sum(post_propor)) 
plot(x,post,type="l",ylab="Density",lty=1,xlim=c(0,8),ylim=c(0,2.5),
     lwd=2,col="gray",
     main=paste("Prior is Inverse-Gamma(5,3)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend("topright",c("prior","Likelihood","Posterior"), 
	lty=c(3,2,1),col=c("black","black","gray"))

# INVERSE GAMMA (3,3)
prior = densigamma(x,3,3) #scale parameter 3, shape parameter 3
like = densigamma(x,sum(y^2),n/2-1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) 

plot(x,post,type="l",ylab="Density",lty=1,xlim=c(0,8),
     ylim=c(0,2.5),lwd=2,col="gray",
     main=paste("Prior is Inverse-Gamma(3,3)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend("topright",c("prior","Likelihood","Posterior"), 
    lty=c(3,2,1),col=c("black","black","gray"))

# INVERSE GAMMA (1,3)
prior = densigamma(x,1,3)  
like = densigamma(x,sum(y^2),n/2-1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) 

plot(x,post,type="l",ylab="Density",lty=1,xlim=c(0,8),
     ylim=c(0,2.5),lwd=2,col="gray",
     main=paste("Prior is Inverse-Gamma(1,3)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend("topright",c("prior","Likelihood","Posterior"), 
     lty=c(3,2,1),col=c("black","black","gray"))

#--------- INVERSE GAMMA (0,3) -- VERY LOW PRECISION-----------#
prior = densigamma(x,0.1,3) 
like = densigamma(x,sum(y^2),n/2-1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) 

plot(x,post,type="l",ylab="Density",lty=1,xlim=c(0,8),
     ylim=c(0,2.5),lwd=2,col="gray",
     main=paste("Prior is Inverse-Gamma(0.1,3)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend("topright",c("prior","Likelihood","Posterior"), 
     lty=c(3,2,1),col=c("black","black","gray"))

#---------NORMAL DISTRIBUTION: MEAN UNKNOWN, VARIANCE KNOWN----------#
x<- seq(-5,5,by=.001) #length(x)=10001
par(mfrow=c(2,2),lwd=2,mar=c(3,3,4,4),cex.axis=.6)

# Normal--High precision
prior = dnorm(x,0,0.3) 
like = dnorm(x,0,1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,1.5),
     main=paste("Prior is Normal (0,0.3)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend(1.6,1.58,c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Normal--Intermediate precision
prior = dnorm(x,0,0.5) 
like = dnorm(x,0,1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) 

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,1.5),
     main=paste("Prior is Normal (0,0.5)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend(1.6,1.58,c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Normal--low precision
prior = dnorm(x,0,1.2) 
like = dnorm(x,0,1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) 

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,1.5),
     main=paste("Prior is Normal (0,1.2)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend(1.6,1.58,c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Normal--very low precision
prior = dnorm(x,0,3) 
like = dnorm(x,0,1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) 

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,1.5),
     main=paste("Prior is Normal (0,3)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend(1.6,1.58,c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))
# End

#----------------UNIFORM PRIOR--------------#
x<- seq(-5,5,by=.001) 

par(mfrow=c(2,2),lwd=2,mar=c(3,3,4,4),cex.axis=.6)

# Prior: Uniform --High precision
prior = dunif(x,min=-1,max=1)
like = dnorm(x,0,1) #given parameter, the distribution of data
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: sum(post)*0.001 = 1

plot(x,post,type="l",ylab="",lty=1,lwd=2,col="gray",ylim=c(0,0.8),main=paste("Prior is Uniform(-1,1)"))
 lines(x,like,lty=2,lwd=2)
 lines(x,prior,lty=3,lwd=2)
 legend("topright",c("prior","Likelihood","Posterior"),lty=c(3,2,1),col=c("black","black","gray"))

# Prior: Uniform --Intermediate precision
prior = dunif(x,min=-3,max=1)
like = dnorm(x,0,1) #given parameter, the distribution of data
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: sum(post)*0.001 = 1

plot(x,post,type="l",ylab="",lty=1,lwd=2,col="gray",ylim=c(0,0.8),main=paste("Prior is Uniform(-3,1)"))
 lines(x,like,lty=2,lwd=2)
 lines(x,prior,lty=3,lwd=2)
 legend("topright",c("prior","Likelihood","Posterior"),lty=c(3,2,1),col=c("black","black","gray"))

# Prior: Uniform --Low precision
prior = dunif(x,min=-5,max=2)
like = dnorm(x,0,1) #given parameter, the distribution of data
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: sum(post)*0.001 = 1

plot(x,post,type="l",ylab="",lty=1,lwd=2,col="gray",ylim=c(0,0.8),main=paste("Prior is Uniform(-5,2)"))
 lines(x,like,lty=2,lwd=2)
 lines(x,prior,lty=3,lwd=2)
 legend("topright",c("prior","Likelihood","Posterior"),lty=c(3,2,1),col=c("black","black","gray"))
 
# Prior: Uniform --Very low precision
prior = dunif(x,min=-8,max=10)
like = dnorm(x,0,1) #given parameter, the distribution of data
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: sum(post)*0.001 = 1

plot(x,post,type="l",ylab="",lty=1,lwd=2,col="gray",ylim=c(0,0.8),main=paste("Prior is Uniform(-8,10)"))
 lines(x,like,lty=2,lwd=2)
 lines(x,prior,lty=3,lwd=2)
 legend("topright",c("prior","Likelihood","Posterior"),lty=c(3,2,1),col=c("black","black","gray"))
 # End

#--------------POISSON DISTRIBUTION-------------#	

x<- seq(0.1,10,by=.001) #length(x)=9901;
y=c(1,3,4,7)
n=4

par(mfrow=c(2,2),lwd=2,mar=c(3,3,4,4),cex.axis=.6)

# Prior: Gamma(10,0.2)--High precision
prior = dgamma(x,shape=10,scale=0.2)
like=dgamma(x,shape=sum(y)+1,scale=1/n)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,1.1),
     main=paste("Prior is Gamma(10,0.2)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend(5.5,1.15,c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Prior: Gamma(8,0.5)--Intermediate precision
prior = dgamma(x,shape=8,scale=0.5)
like=dgamma(x,shape=sum(y)+1,scale=1/n)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,1.1),
     main=paste("Prior is Gamma(8,0.5)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend(5.5,1.15,c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Prior: Gamma(3,1)--Low precision
prior = dgamma(x,shape=3,scale=1)
like=dgamma(x,shape=sum(y)+1,scale=1/n)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,1.1),
     main=paste("Prior is Gamma(3,1)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend(5.5,1.15,c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Prior: Gamma(2.1,3)--Very low precision
prior = dgamma(x,shape=2.1,scale=3)
like=dgamma(x,shape=sum(y)+1,scale=1/n)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,1.1),
     main=paste("Prior is Gamma(2.1,3)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend(5.5,1.15,c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))
# End

#-------BINOMIAL DISTRIBUTION-------------#
x<- seq(0.001,1,by=.001)    #length(x)=1000; parameter
y=9
n=20

par(mfrow=c(2,2),lwd=2,mar=c(3,3,4,4),cex.axis=.6)

# Prior: Beta--High Precision
prior = dbeta(x,17,19)
like=dbeta(x,y+1,n-y+1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,6),
     main=paste("Prior is Beta(17,19)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend("topright",c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Prior: Beta--Intermediate Precision
prior = dbeta(x,7,8)
like=dbeta(x,y+1,n-y+1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,6),
     main=paste("Prior is Beta(7,8)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend("topright",c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Prior: Beta--Low Precision
prior = dbeta(x,3,4)
like=dbeta(x,y+1,n-y+1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,6),
     main=paste("Prior is Beta(3,4)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend("topright",c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))

# Prior: Beta--Very low Precision
prior = dbeta(x,1.2,1.5)
like=dbeta(x,y+1,n-y+1)
post_propor = prior*like 
post = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

plot(x,post,type="l",ylab="Density",lty=1,lwd=2,col="gray",ylim=c(0,6),
     main=paste("Prior is Beta(1.2,1.5)"))
lines(x,like,lty=2,lwd=2)
lines(x,prior,lty=3,lwd=2)
legend("topright",c("prior","Likelihood","Posterior"), lty=c(3,2,1),col=c("black","black","gray"))
# End

#---------MULTINOMIAL DISTRIBUTION--------------#
library(MCMCpack)    	# For dirichlet distribution
library(MVA)     			# For 3D plot

x<-matrix(rep(0,1000*3),nrow=1000)
x[,1]<- sample(seq(0.001,1,by=.001)) #length(x[,1])=1000; parameter theta1
x[,2]<- runif(1000,0,1-x[,1])                 #length(x[,2])=1000; parameter theta2
x[,3]<-1-x[,1]-x[,2]                                #length(x[,3])=1000; parameter theta3
y1=3
y2=5 
y3=6
n=14

# Prior: Dirichlet (High Precision)
par(mfrow=c(3,3),lwd=2,mar=c(3,3,4,4),cex.axis=.6)

prior1 =  ddirichlet(x, alpha=c(5,10,14))
like1=ddirichlet(x,alpha=c(y1+1,y2+1,y3+1))
post_propor = prior1*like1 
post1 = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

scatterplot3d(x[,1],x[,2],prior1,type="p",angle=55,color="black", xlab=expression(pi[1]), 
             ylab=expression(pi[2]),zlab="Prior Density",main="Dirichlet Prior for Multinomial Distribution--High Precision")
scatterplot3d(x[,1],x[,2],like1,type="p",angle=55,color="black", xlab=expression(pi[1]), 
             ylab=expression(pi[2]),  zlab="Likelihood",main="Likelihood")
scatterplot3d(x[,1],x[,2],post1,type="p",angle=55,color="black", xlab=expression(pi[1]), 
             ylab=expression(pi[2]), zlab="Posterior Density",main="Posterior")

# Prior: Dirichlet (Intermediate Precision)

prior2 =  ddirichlet(x, alpha=c(3,4,3))
like2=ddirichlet(x,alpha=c(y1+1,y2+1,y3+1))
post_propor = prior2*like2 
post2 = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

scatterplot3d(x[,1],x[,2],prior2,type="p",angle=55,color="black", xlab=expression(pi[1]), 
             ylab=expression(pi[2]),zlab="Prior Density",main="Dirichlet Prior for Multinomial Distribution--Moderate Precision")
scatterplot3d(x[,1],x[,2],like2,type="p",angle=55,color="black", xlab=expression(pi[1]), 
             ylab=expression(pi[2]),
             zlab="Likelihood",main="Likelihood")
scatterplot3d(x[,1],x[,2],post2,type="p",angle=55,color="black",xlab=expression(pi[1]), 
             ylab=expression(pi[2]), zlab="Posterior Density",main="Posterior")

# Prior: Dirichlet  (Low Precision)

prior3 =  ddirichlet(x, alpha=c(2,2,2))
like3=ddirichlet(x,alpha=c(y1+1,y2+1,y3+1))
post_propor = prior3*like3 
post3 = post_propor / (0.001*sum(post_propor)) #Approx. posterior density: #sum(post)*0.001 = 1

scatterplot3d(x[,1],x[,2],prior3,type="p",angle=55,color="black", xlab=expression(pi[1]), 
             ylab=expression(pi[2]),zlab="Prior Density",main="Dirichlet Prior for Multinomial Distribution--Low Precision")
scatterplot3d(x[,1],x[,2],like3,type="p",angle=55,color="black", xlab=expression(pi[1]), ylab=expression(pi[2]),
             zlab="Likelihood",main="Likelihood")
scatterplot3d(x[,1],x[,2],post3,type="p",angle=55,color="black", xlab=expression(pi[1]), ylab=expression(pi[2]),
             zlab="Posterior Density",main="Posterior")
# End