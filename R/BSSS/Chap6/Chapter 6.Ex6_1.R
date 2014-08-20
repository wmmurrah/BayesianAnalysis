#-----------------------------------------------------------------------------------------
# EXAMPLE 6:1: BAYESIAN MULTIPLE REGRESSION WITH NON-INFORMATIVE PRIORS
#                Program Steps
# 1. Read in data and apply tranformations or handle missing data as needed
# 2. Call MCMCregress
# 3. Default prior mean b0=0 precision B0=0 tells MCMCregress to use to use an improper
#     uniform prior
# 4. Summarize and obtain diagnostic plots
#-----------------------------------------------------------------------------------------

# install.packages("MCMCpack")
# install.packages("BMA")
# install.packages("coda")
require(MCMCpack)
require(coda)
require (BMA)

#  Read in data PISA2009.csv

datafile <- read.csv(file.choose(),header=T)
datafile9 <- subset(datafile, select=c(rcomb1, gender, native,  slang,  ESCS,
       JOYREAD, DIVREAD, MEMOR, ELAB, CSTRAT))
head(datafile9)
nrow(datafile9)
datafile9<-na.omit(datafile9)
nrow(datafile9)


# lm ----------------------------------------------------------------------

lm.fit <- lm(rcomb1~gender+native+ slang+ESCS+
               JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT,
             data=datafile9)
summary(lm.fit)

# MCMCpack ----------------------------------------------------------------

FullModel <- MCMCregress(rcomb1~gender+native+ slang+ESCS+
       JOYREAD+ DIVREAD+ MEMOR+ ELAB+ CSTRAT,
       data=datafile9,burnin=5000,mcmc=100000,thin=10,b0=0,B0=0)
plot(FullModel)
autocorr.plot(FullModel)
dev.off()
summary(FullModel)

# Diagnostics
geweke.diag(FullModel, frac1=0.1, frac2=0.5)  
heidel.diag(FullModel,eps=0.1,pvalue=0.05)  
raftery.diag(FullModel,q=0.5,r=0.05,s=0.95,converge.eps=0.001)  


# STAN --------------------------------------------------------------------
library(rstan)
df <- datafile9
N <- as.numeric(nrow(df))
K <- as.numeric(ncol(df))-1
rcomb1 <- df$rcomb1
x  <- as.matrix(df[ ,2:10])
datlist <- c('N', 'K','rcomb1', 'x')

stan.fullmod <- '
data{
  int<lower=0> N;                 // number of cases
  int<lower=0> K;                 // number of predictors
  vector[N] rcomb1;               // outcome
  matrix[N,K] x;                  // matrix of predictors
}
parameters {
  real b0;                        // intercept
  vector[K] beta;                 // vector of regression coefficients
  real<lower=0, upper=100> sigma; // residual standard error
}
model {
  b0 ~ uniform(-1000, 1000);
  beta[K] ~ uniform(-1000, 1000);
  rcomb1 ~ normal(b0 + x * beta, sigma); // likelihood
}
'
stan.fullfit <- stan(model_code = stan.fullmod, data=datlist, iter=1000,
                     chains=2)

traceplot(stan.fullfit)
stan.fullfit
# run.jags ----------------------------------------------------------------




# end ---------------------------------------------------------------------

