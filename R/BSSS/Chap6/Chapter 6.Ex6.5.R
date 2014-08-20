#------------- EXAMPLE 6.5 BAYESIAN LOGISTIC REGRESSION---------#

install.packages("MCMCpack")
install.packages("BMA")

#Multiple Regression Analysis :
library(MCMCpack)
datafile <- read.csv(file.choose(),header=T)
datafile9 <- subset (datafile, select=c(private, gender, native, slang, ESCS))

head(datafile9)
nrow(datafile9)
datafile9<-na.omit(datafile9)
nrow(datafile9)

#FullModel
FullModel <- MCMClogit(private~gender+ native+  slang+  ESCS
       ,data=datafile9,marginal.likelihood="Laplace",mcmc=100000,b0=0,
	B0=(.01))
plot(FullModel) # Produces the convergence plots and the posterior densities
dev.off()
summary(FullModel)

FullModel_inf <- MCMClogit(private~gender+ native+  slang+  ESCS
	,data=datafile9,marginal.likelihood="Laplace",mcmc=100000,
	b0=c(-3.4118, 0.2513,  0.6171,  -0.3057, 0.2446 ),
	B0=c( 3.454347, 27.272662, 3.018396, 6.184967, 93.873638))
plot(FullModel_inf) # Produces the convergence plots and the posterior densities
dev.off()
summary(FullModel_inf)

#Convergence Diagnostics
library(coda)
geweke.diag(FullModel_inf, frac1=0.1, frac2=0.5)  
heidel.diag(FullModel_inf,eps=0.1,pvalue=0.05) 
raftery.diag(FullModel_inf,q=0.5,r=0.05,s=0.95,converge.eps=0.001) 

geweke.diag(FullModel_inf, frac1=0.1, frac2=0.5) 
heidel.diag(FullModel_inf,eps=0.1,pvalue=0.05) 
raftery.diag(FullModel_inf,q=0.5,r=0.05,s=0.95,converge.eps=0.001) 

#Bayesian Model Averaging
library(BMA)
attach(datafile9)
bma=bic.glm(cbind(gender, native,  slang,  ESCS),private,
glm.family=binomial("logit"), strict=FALSE,OR=20)
summary(bma)
plot(bma) # Plots of BMA posterior distributions
imageplot.bma(bma) # The image plot shows which predictors are included in each model

