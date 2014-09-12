#-----------EXAMPLE 7.1: A COMPARATIVE EXAMPLE OF MISSING DATA METHODS----------#
# Read in data file 
datafile <- read.csv(file.choose(),header=T)
datafile9 <- subset(datafile, select=c(rcomb1, gender, native,  slang,  ESCS,
       JOYREAD, DIVREAD, MEMOR, ELAB, CSTRAT))
head(datafile9)
nrow(datafile9)

#################################################
# Norm (Data Augmentation)
#################################################
install.packages("norm")
require("norm")
library(MCMCpack)

#Transfer data into a matrix 
a <- as.matrix(datafile9)

#Preliminary manipuation
s  <-  prelim.norm(a)

#find the MLE for a starting value
thetahat <- em.norm(s)

#set random number generator seed
rngseed(1234567) 
  
#Run the NORM program, specifiying the data object, start value, and number of steps
theta <- da.norm(s,thetahat,steps=20,showits=TRUE)

#Generated Imputed data set using theta
IMP.norm<-imp.norm(s,theta,a )
IMP.norm<-as.data.frame(IMP.norm)

#Conduct Bayesian Regression using imputed data
noninf.norm<- MCMCregress(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, data=IMP.norm, prior.m=0, prior.p=0)

summary(noninf.norm)

#Conduct Bayesian Regression with informative prior
inf.norm<- MCMCregress(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, data=IMP.norm, marginal.likelihood="Chib95",mcmc=10000,
	 b0=c( 477.6,15.3, 4.8, 20.4, 36.2, 21.1, 1.7, -17.8,  -13.3,  25.5  ),
	 B0=c( 0.0435, 0.1089, 0.0247, 0.0298, 0.3824, 0.3202,  
			0.4207, 0.2397, 0.2306, 0.1950  ))

summary(inf.norm)

#Conduct frequentist regression 
freq.norm<- lm(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, data=IMP.norm)

summary(freq.norm)

##################################################
# Chained equations
##################################################
install.packages("mi")
library(mi) 
#Impute data
IMP <- mi(datafile9, n.imp=5,n.iter=300, add.noise=FALSE,seed=12345 )

#Conduct Bayesian Regression with non informative prior
noninf<- Bayesglm.mi.ce(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP, prior.m=0, prior.p=0)

#Conduct Bayesian Regression with informative prior
inf<- Bayesglm.mi.ce(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP, 
	 prior.m=c( 477.6,15.3, 4.8, 20.4, 36.2, 21.1, 1.7, -17.8,  -13.3,  25.5  ),
	 prior.p=c( 0.0435, 0.1089, 0.0247, 0.0298, 0.3824, 0.3202,  
			0.4207, 0.2397, 0.2306, 0.1950  ))

#Conduct frequentist regression 
b<-lm.mi(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP)
summary(b)
								
#################################################
# Amelia (EM Boostrap)
#################################################
install.packages("Amelia")
require(Amelia)

#Set bounds on variables to be imputed. These should be determined by the actual distributions of each variable.
bds <- matrix(c(1,-3.5,3.5, 2,1,4, 3,-3.5,3.5, 4,-3.5,3.5), nrow = 4, ncol=3, byrow=TRUE)

#Run the AMELIA program, specifiying the data object, the number of imputed data sets desired, and the bounds for imputed variables
IMP.amelia <- amelia(x=datafile9,m=5, bounds=bds)

#Conduct Bayesian Regression using imputed data
noninf.amelia<- Bayesglm.mi.amelia(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP.amelia, prior.m=0, prior.p=0)

#Conduct Bayesian Regression with informative prior
inf.amelia<- Bayesglm.mi.amelia(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP.amelia, 
	 prior.m=c( 477.6,15.3, 4.8, 20.4, 36.2, 21.1, 1.7, -17.8,  -13.3,  25.5  ),
	 prior.p=c( 0.0435, 0.1089, 0.0247, 0.0298, 0.3824, 0.3202,  
			0.4207, 0.2397, 0.2306, 0.1950  ))

#Conduct frequentist regression 
freq.amelia<-freq.mi.amelia(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP.amelia)
       
##################################################
# BaBooN (Bayesian Predictive Mean Matching)           
##################################################
install.packages("BaBooN")
require(BaBooN)

#Run BaBooN program, specifying the data object, number of iterations desired, and the number of imputed data sets desired. 
IMP.bbpmm <- BBPMM(datafile9, nIter=5, M=5)

#Conduct Bayesian Regression using imputed data
noninf.bbpmm<- Bayesglm.mi.bbpmm(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP.bbpmm, prior.m=0, prior.p=0)

#Conduct Bayesian Regression with informative prior
inf.bbpmm<- Bayesglm.mi.bbpmm(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP.bbpmm, 
	 prior.m=c( 477.6,15.3, 4.8, 20.4, 36.2, 21.1, 1.7, -17.8,  -13.3,  25.5  ),
	 prior.p=c( 0.0435, 0.1089, 0.0247, 0.0298, 0.3824, 0.3202,  
			0.4207, 0.2397, 0.2306, 0.1950  ))

#Conduct frequentist regression 
freq.bbpmm<-freq.mi.bbpmm(formula=rcomb1~gender+ native+  slang+  ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT, mi.object=IMP.bbpmm)

####### FUNCTIONS FOR COMBINING FILES ##############      

#  MCMCreg.mi.ce   For chained equations

MCMCreg.mi.ce <-function(formula, mi.object, prior.m, prior.p)
{
    call   <-match.call()	
    library(MCMCpack)
    m      <- m(mi.object)
    result <- vector( "list", m )
    names( result ) <- as.character(paste( "Chain", seq( m ), sep = "" ))
    mi.data <- mi.completed(mi.object)

    for ( i in 1:m ) {
	as.data.frame(do.call(rbind, mi.data[[i]]))
      result[[i]] <- MCMCregress(formula,data=mi.data[[i]],marginal.likelihood="none"
,mcmc=10000,b0=prior.m,B0=prior.p)
    }
    coef   <- vector( "list", m )
    se     <- vector( "list", m )
    for( j in 1:m ) {
      coef[[j]]<- lapply( result, summary )[[ j ]]$statistics[ ,1]
      se[[j]]  <- lapply( result, summary )[[ j ]]$statistics[ ,2]
    }
    pooled <- mi.pooled(coef, se)
    print(pooled)
}

# MCMCreg.mi.bbpmm  For Bayesian bootstrap pred. mean matching

MCMCreg.mi.bbpmm <-function(formula, mi.object, prior.m, prior.p)
{
    call   <-match.call()	
    library(MCMCpack)
    m      <- mi.object$M
    result <- vector( "list", m )
    names( result ) <- as.character(paste( "Chain", seq( m ), sep = "" ))
    mi.data <- mi.object$impdata


    for ( i in 1:m ) {
	as.data.frame(do.call(rbind, mi.data[[i]]))
      result[[i]] <- MCMCregress(formula,data=mi.data[[i]],marginal.likelihood="none"
,mcmc=10000,b0=prior.m,B0=prior.p)
    }
    coef   <- vector( "list", m )
    se     <- vector( "list", m )
    for( j in 1:m ) {
      coef[[j]]<- lapply( result, summary )[[ j ]]$statistics[ ,1]
      se[[j]]  <- lapply( result, summary )[[ j ]]$statistics[ ,2]
    }
    pooled <- mi.pooled(coef, se)
    print(pooled)
}

# MCMCreg.mi.amelia  For the EM boostrap

MCMCreg.mi.amelia <-function(formula, mi.object, prior.m, prior.p)
{
    call   <-match.call()	
    library(MCMCpack)
    m      <- mi.object$m
    result <- vector( "list", m )
    names( result ) <- as.character(paste( "Chain", seq( m ), sep = "" ))
    mi.data <- mi.object$imputations


    for ( i in 1:m ) {
	as.data.frame(do.call(rbind, mi.data[[i]]))
      result[[i]] <- MCMCregress(formula,data=mi.data[[i]],marginal.likelihood="none"
,mcmc=10000,b0=prior.m,B0=prior.p)
    }
    coef   <- vector( "list", m )
    se     <- vector( "list", m )
    for( j in 1:m ) {
      coef[[j]]<- lapply( result, summary )[[ j ]]$statistics[ ,1]
      se[[j]]  <- lapply( result, summary )[[ j ]]$statistics[ ,2]
    }
    pooled <- mi.pooled(coef, se)
    print(pooled)
}

# freq.mi.bbpmm # 

freq.mi.bbpmm <-function(formula, mi.object)
{
    call   <-match.call()	
     m      <- mi.object$M
    result <- vector( "list", m )
    names( result ) <- as.character(paste( "Chain", seq( m ), sep = "" ))
    mi.data <- mi.object$impdata

    for ( i in 1:m ) {
	as.data.frame(do.call(rbind, mi.data[[i]]))
      result[[i]] <- lm(formula,data=mi.data[[i]])
    }
    coef   <- vector( "list", m )
    se     <- vector( "list", m )
    for( j in 1:m ) {
      coef[[j]]<- lapply( result, summary )[[ j ]]$coefficients[ ,1]
      se[[j]]  <- lapply( result, summary )[[ j ]]$coefficients[ ,2]
    }
    pooled <- mi.pooled(coef, se)
    print(pooled)
}

# freq.mi.amelia# 

freq.mi.amelia <-function(formula, mi.object)
{
    call   <-match.call()	
    m      <- mi.object$m
    result <- vector( "list", m )
    names( result ) <- as.character(paste( "Chain", seq( m ), sep = "" ))
    mi.data <- mi.object$imputations

    for ( i in 1:m ) {
	as.data.frame(do.call(rbind, mi.data[[i]]))
      result[[i]] <- lm(formula,data=mi.data[[i]])

    }
    coef   <- vector( "list", m )
    se     <- vector( "list", m )
    for( j in 1:m ) {
      coef[[j]]<- lapply( result, summary )[[ j ]]$coefficients[ ,1]
      se[[j]]  <- lapply( result, summary )[[ j ]]$coefficients[ ,2]
    }
    pooled <- mi.pooled(coef, se)
    print(pooled)
}
