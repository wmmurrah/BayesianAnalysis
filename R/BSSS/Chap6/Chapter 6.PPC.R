#------------------------POSTERIOR PREDICTIVE CHECK ----------------------------------#

library(gridExtra)
library(coda)
require(MCMCpack)

#set.seed(515)
###These grab the data from the MCMC object
beta.df <- data.frame(as.matrix(FullModel_inf))
beta.mat <- as.matrix(beta.df[sample(8000:10000,1000),])?beta.non.df <- data.frame(as.matrix(FullModel))
beta.non.mat <- as.matrix(beta.non.df[sample(8000:10000,1000),])

ppc <- function(x, df1) {

obs <- as.matrix(cbind(Intercept = 1,datafile9[,2:10])) ## Independent variables in regression model
betas <- as.matrix(df1[x,1:10]) ## betas drawn from the posterior distribution
y.rep.ev <- obs %*% betas       ## expected value y
y.rep.var <- df1[x,11] ## variance of y

#print(c(y.rep.ev,y.rep.var))
y.rep <- obs %*% betas + rnorm(nrow(datafile9),0,sqrt(df1[x,11])) ## replicated data
val1 <- ((y.rep - y.rep.ev)^2)/y.rep.var ##y.rep chi.square ## Replicated Chi Square Statistic
val2 <- ((datafile9$rcomb1 - y.rep.ev)^2)/y.rep.var ## Observed Chi Square Statistic
return(list(y.rep,c(chi.rep = sum(val1),chi.obs = sum(val2))))

}##END function ppc

## SIMULATED THE REPLICATED DATA AND OBTAIN THE DISCREPANCY STATISTIC
ppc.data <- lapply(1:1000,ppc, df1 = beta.mat)
ppc.non.data <- lapply(1:1000,ppc, df1 = beta.non.mat)
posterior.check <- t(sapply(ppc.data, function(x) return(x[[2]])))
posterior.non.check <- t(sapply(ppc.non.data, function(x) return(x[[2]])))
chi.obs <- posterior.non.check[,2]
chi.rep <- posterior.non.check[,1]
chi.discrepancy <- posterior.non.check[,2] - posterior.non.check[,1] ##get the difference between y.obs and y.rep
p.value <- round(length(which(chi.discrepancy < 0))/length(chi.discrepancy),3)
range <- max(chi.discrepancy) - min(chi.discrepancy) ##range

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
geom_text(aes(x = 65 + min(chi.obs,chi.rep) , y = - 2 + max(chi.obs,chi.rep), label = paste("p-value =",p.value)), size = 4.5) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
ylim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) + 
xlim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) +
xlab("Observed") + ylab("Replicated")
#dev.off()

##DISPLAY POSTERIOR CHECKING PLOTS & INFORMATION
hist(liks.dif, xlab = expression(Chi["obs"]^2 - Chi["rep"]^2), main = "")
abline(v = 0, lty = 2, lwd = 2)
text(x = min(liks.dif) + range/5, y = 200, label = paste("p-value = ",p.value,sep = ""))
par(ask = TRUE)
qplot(chi.obs,chi.rep, shape = I(1)) + 
geom_abline(slope = 1, intercept = 0) +
theme_bw() +
geom_text(aes(x = 65 + min(chi.obs,chi.rep) , y = - 2 + max(chi.obs,chi.rep), label = paste("p-value =",p.value)), size = 4.5) +
#theme(panel.grid.major = element_blank()) +
ylim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) + 
xlim(c(-7 + min(chi.obs,chi.rep),max(chi.obs,chi.rep))) +
xlab("Observed") + ylab("Replicated")

###
pdf(file=' ') #####ADD FILE PATH######
hist(FullModel_inf[,11], main = "Histogram of posterior draws for Sigma", xlab = "Sigma")
#dev.off()

