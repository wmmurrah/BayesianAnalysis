#----------------------------------------------------------------------------
# EXAMPLE 8.1: BAYESIAN RANDOM EFFECTS ANOVA 
#        Program Steps
# 1. Read in data
# 2. Call MCMChregress from MCMCpack and specify model
# 3. For random effects ANOVA regress outcome on the intercept
#     and define the intercept as random
# 4. Set starting values and priors
# 5. Summarize the values in the model object (here model_inf) 
# 6. Obtain plots
#
# Remaining code shows MCMChregress for informative priors
# and code to obtain posterior predictive checks.
#---------------------------------------------------
# install.packages("MCMCpack")
# install.packages("coda")
require(MCMCpack)
require(coda)

hlmdata <- read.csv(file.choose(),header=T)     

model_inf <- MCMChregress(fixed=rcomb1~1, random=~1, 				          
	        group="SCHOOLID", data=hlmdata, burnin=1000, mcmc=10000, thin=10, verbose=1,
              seed=2012, beta.start=0, sigma2.start=1,
              Vb.start=1, mubeta=500, Vbeta=100,
              r=1, R=diag(100,1,1), nu=0.001, delta=0.001)
              
summary(model_inf$mcmc[,1])  #Posterior Mean and SD for the Fixed Effect

summary(model_inf$mcmc[,167:169]) #Posterior Mean and SD for variance components

plot(model_inf$mcmc[,c(1,167)])     #Trace Plots and Density Plots for Selected Fixed Effects

geweke.plot(model_inf$mcmc[,c(1,167)])
geweke.diag(model_inf$mcmc[,c(1,167)])
heidel.diag(model_inf$mcmc[,c(1,167)])
raftery.diag(model_inf$mcmc[,c(1,167)])


# STAN --------------------------------------------------------------------
library(rstan)
model <- '
data{ 
   int<lower=0> G;
   int<lower=0> N;
   vector[N] y;
   
}
parameters {
   real mu;
   real<lower=0> tau;
   real[G] eta;
}
model {
   mu ~ inv_gamma(500, 100);
   

   y ~ normal(mu + eta, sigma)
}
'

 #----------EXAMPLE 8.2: Multilevel Model: Informative Priors from PISA 2000---------#
#install.packages("MCMCpack")
require(MCMCpack)

hlmdata <- read.csv(file.choose(),header=T)     

model_inf <- MCMChregress(fixed=rcomb1~JOYREAD+gender+MEMOR+DISCLIMA+SCHSIZE+TCSHORT+JOYREAD:TCSHORT, 
              random=~JOYREAD+gender+MEMOR+DISCLIMA, 				          
	        group="SCHOOLID",data=hlmdata, burnin=1000, mcmc=100000, thin=100,verbose=1,
              seed=2012, beta.start=0, sigma2.start=1,
              Vb.start=1, mubeta=c(484.4,27.84, 14.75,-1.17,-6.63, 0.59,-11.24, -1.45), Vbeta=2*c(66.85,
              3.47, 11.17, 2.58, 3.57, 0.33, 26.23,3.43),
              r=5, R=diag(c(2000,50,1,1,100)), nu=0.001, delta=0.001)


summary(model_inf$mcmc[,1:8])  #Posterior Mean and SD for Selected Fixed Effects

# varaince components
ncol(model_inf$mcmc)
summary(model_inf$mcmc[,834:860])  

#  #Trace Plots and Density Plots for Selected Fixed Effects
plot(model_inf$mcmc[,c(1,2,6)])    

#------------------------POSTERIOR PREDICTIVE CHECK ----------------------------------#

library(ggplot2)
library(gridExtra)
library(coda)
require(MCMCpack)

set.seed(515)
###These grab the data from the MCMC object
df.hlm <- data.frame(as.matrix(model_inf$mcmc)) 
var.list <- c(".Intercept.","JOYREAD","gender","MEMOR","DISCLIMA","SCHSIZE","TCSHORT","JOYREAD.TCSHORT")
var.list <- paste("beta.",var.list, sep = "")
eff.list <- c(".Intercept.","JOYREAD","gender","MEMOR","DISCLIMA")
eff.list <- sapply(eff.list,function(x) paste(x,".",x,sep = ""))
eff.list <- paste("VCV.",eff.list, sep = "")
df.hlm <- df.hlm[,c(var.list,eff.list,"sigma2")]

##function to calculate the chi-square discrepancy statistic
discrep.stat <- function(x) {

y.rep <- y.posterior[[3*x-2]] ##creates all the y.reps for simulation number x
nas <- which(is.na(y.rep)) ##remove the NA's due to missing values in the X matrix
y.rep <- y.rep[-nas]
y.rep.ev <- y.posterior[[3*x-1]] ##creates the expected values based upon the theta draws and the X matrix
y.rep.ev <- y.rep.ev[-nas]
y.rep.var <- y.posterior[[3*x]] ##creates the variance based on the theta draws
obs <- hlmdata$rcomb1 ##grabs the observed values
obs <- obs[-nas]

val1 <- ((y.rep - y.rep.ev)^2)/y.rep.var ##y.rep chi.square
val2 <- ((obs - y.rep.ev)^2)/y.rep.var ##y.obs chi.square

return(c(sum(val1),sum(val2)))

}##END function lik.ratio.chi.square

##Function to calculate the posterior y.reps
simulate.posterior <- function(x, obs = as.matrix(hlmdata), betas = as.matrix(df.hlm)) {
obs.cols <- match(c("JOYREAD","gender","MEMOR","DISCLIMA","SCHSIZE","TCSHORT"),attr(obs,"dimnames")[[2]] )
obs.x <- obs[,obs.cols]
obs.x <- cbind(obs.x,as.matrix(obs[,obs.cols[1]] * obs[,obs.cols[6]]))
pred.x <- cbind(1,obs.x)  %*% as.numeric(betas[x,1:8])
size <- nrow(obs)
pred.x.rand <- pred.x + rnorm(size,0,sqrt(betas[x,9])) + ## random error from Intercept
rnorm(size,0,sqrt(betas[x,10])) + ## random error from JOYREAD
rnorm(size,0,sqrt(betas[x,11])) + ## random error from gender
rnorm(size,0,sqrt(betas[x,12])) + ## random error from MEMOR
rnorm(size,0,sqrt(betas[x,13])) + ## random error from DISCLIMA
rnorm(size,0,sqrt(betas[x,14])) ## overall random error
var.x = sum(betas[x,9:14])
return(list(pred.x.rand, pred.x, var.x))
}##END function simulate.posterior

##These are the lines that actually run the above functions
y.posterior <- sapply(1:1000,simulate.posterior) ##create the 1000 posterior simulations of sample size n = 5322
posterior.check <- t(sapply(1:1000,discrep.stat)) ##create the 1000 discrepancy statistics for both y.rep and y.obs
chi.obs <- posterior.check[,2]
chi.rep <- posterior.check[,1]
chi.discrepancy <- posterior.check[,2] - posterior.check[,1] ##get the difference between y.obs and y.rep
p.value <- round(length(which(chi.discrepancy < 0))/length(chi.discrepancy),3)
range <- max(chi.discrepancy) - min(chi.discrepancy) ##range
x.breaks <- c(round(median(chi.discrepancy),1),0,seq(round(range(chi.discrepancy)[1],0),round(range(chi.discrepancy)[2],0),round(range/4,0)))
x.breaks <- x.breaks[order(x.breaks)]

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
pdf(file=' ') #####ADD FILE PATH######
qplot(chi.obs,chi.rep, shape = I(1)) + 
geom_abline(slope = 1, intercept = 0) +
theme_bw() +
geom_text(aes(x = 100 + min(chi.obs,chi.rep) , y = - 2 + max(chi.obs,chi.rep), label = paste("p-value =",p.value)), size = 4.5) +
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
geom_text(aes(x = 100 + min(chi.obs,chi.rep) , y = - 2 + max(chi.obs,chi.rep), label = paste("p-value =",p.value)), size = 4.5) +
#theme(panel.grid.major = element_blank()) +
ylim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) + 
xlim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) +
xlab("Observed") + ylab("Replicated")

##PRINT OUT MEDIAN DIFFERENCE and P-VALUE
print(data.frame(Median.Difference = median(chi.discrepancy), p.value = p.value))