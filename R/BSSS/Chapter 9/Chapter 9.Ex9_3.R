#--------------------------------------------------------------
# Example 9.3 Multilevel Path Analysis--Noninformative Prior
# 
# Programming steps the same as in Example 9.2 except that we loop
# over students i and schools g
#--------------------------------------------------------------

install.packages("rjags")
require(rjags)

##Specify model
modelstring = "								
model {										 		
#subject-level models
for (i in 1:nData)
{
for(j in (sch_order[i]):(sch_order[i]))
{
y[i]~dnorm(mu1[i],tau1)
mu1[i]<- a0[j,1] + a1[1]*momedu[i] + a1[2]*dadedu[i]+a1[3]*perteach[i]+
            a2*importnt[i]+a0[j,2]*enjoy[i]

enjoy[i]~dnorm(mu2[i], tau2)
mu2[i]<-b0[j]+b1*perteach[i]

importnt[i]~dnorm(mu3[i], tau3)
mu3[i]<-c0[j]+c1[1]*momedu[i]+c1[2]*perteach[i]+c2*enjoy[i]
}
}

# school level models
for(j in 1:G){
a0[j,1:2]~dmnorm(z[j,1:2],phi[1:2,1:2])

z[j,1]<- m[1]+ A[1]*newmetho[j] + A[2]*enthusia[j] + A[3]*cnsensus[j]
          + A[4]*cndition[j]+ A[5]*encourag[j]
z[j,2]<- m[2]+ B[1]*newmetho[j] + B[2]*enthusia[j] + B[3]*cnsensus[j]
          + B[4]*cndition[j]+ B[5]*encourag[j]

b0[j]~ dnorm(v[j],psi)
v[j]<- m[3]+ C[1]*newmetho[j] + C[2]*enthusia[j] + C[3]*cnsensus[j]
          + C[4]*cndition[j]+ C[5]*encourag[j]

c0[j]~ dnorm(w[j],xi)
w[j]<- m[4]+D[1]*newmetho[j] + D[2]*enthusia[j] + D[3]*cnsensus[j]
          + D[4]*cndition[j]+ D[5]*encourag[j]

encourag[j]~dnorm(r[j],pi)
r[j]<-m[5]+E*enthusia[j]

enthusia[j]~dnorm(p[j],delta)
p[j]<-m[6]+F[1]*newmetho[j]+F[2]*cnsensus[j]+ F[3]*cndition[j]

}

# Prior Specification

# Priors on regression coefficients
m[1:6] ~ dmnorm(u1[1:6], H1[1:6,1:6]) # H1 is prior precision matrix.
A[1:5] ~ dmnorm(u2[1:5], H2[1:5,1:5])
B[1:5] ~ dmnorm(u3[1:5], H3[1:5,1:5])
C[1:5] ~ dmnorm(u4[1:5], H4[1:5,1:5])
D[1:5] ~ dmnorm(u5[1:5], H5[1:5,1:5])
E ~ dnorm(u6,H6)
F[1:3] ~ dmnorm(u7[1:3], H7[1:3,1:3])
a1[1:3]  ~ dmnorm(u8[1:3], H8[1:3,1:3])
a2~dnorm(u9,H9)
b1 ~ dnorm(u10, H10)
c1[1:2]  ~ dmnorm(u11[1:2], H11[1:2,1:2])
c2~dnorm(u12,H12)

#Priors on Precisions
tau1 ~ dgamma(0.001,0.001)
tau2 ~ dgamma(0.001,0.001)
tau3 ~ dgamma(0.001,0.001)

phi[1:2,1:2] ~ dwish(R1[1:2,1:2], 2) # Precision matrix
psi ~ dgamma(0.001,0.001)
xi  ~ dgamma(0.001,0.001)
pi  ~ dgamma(0.001,0.001)
delta  ~ dgamma(0.001,0.001)
}  
"
writeLines(modelstring,con="model.bug")

# READ IN DATA AND PREPARE FOR JAGS 

semdata  <- read.csv(file.choose(),header=TRUE) #browse to select data "PISA2003.semmodel.csv" (using regression imputation).
#colnames(semdata) 
y        <- semdata$mathscor
momedu   <- semdata$momeduc
dadedu   <- semdata$dadeduc
perteach <- semdata$perteach
importnt <- semdata$importnt
enjoy    <- semdata$enjoy
newmetho <- semdata$newmetho
enthusia <- semdata$enthusia
cnsensus <- semdata$cnsensus
cndition <- semdata$cndition
encourag <- semdata$encourag
schoolid <- semdata$schoolid
nData<-NROW(y)

# Calculate the order of the school that each student belongs to
# n[g] is the number of student within the g_th school (unequal school sizes)

sch_id<-unique(schoolid) 
rk<-rank(sch_id)
G<-length(sch_id) #number of schools = G
sch_order<-rep(0,nrow(semdata))
n<-rep(0,G) #initialize n
for(g in 1:G)
{
  n[g]<-sum(schoolid==sch_id[g]) #n[g] is the sample size of school g
  sch_order[schoolid==sch_id[g]]=rep(rk[g],n[g])
}
# Specify parameters
semdata <- list(y=y, G=G, sch_order=sch_order,nData=nData, R1=diag(1,2,2), 
                u1=rep(0,6),u2=rep(0,5),u3=rep(0,5),u4=rep(0,5),u5=rep(0,5),
                u6=rep(0,1),u7=rep(0,3), u8=rep(0,3),u9=rep(0,1),u10=rep(0,1),
                u11=rep(0,2),u12=rep(0,1),
                H1=diag(10^(-1),6,6),H2=diag(10^(-1),5,5),H3=diag(10^(-1),5,5),
                H4=diag(10^(-1),5,5),H5=diag(10^(-1),5,5), 
                H6=diag(10^(-1),1,1),H7=diag(10^(-1),3,3),H8=diag(10^(-1),3,3),
                H9=diag(10^(-1),1,1),H10=diag(10^(-1),1,1),
                H11=diag(10^(-1),2,2),H12=diag(10^(-1),1,1),
                momedu=momedu, dadedu=dadedu, perteach=perteach, importnt= importnt,
                enjoy=enjoy, newmetho=newmetho, enthusia=enthusia, cnsensus=cnsensus,
                cndition=cndition, encourag=encourag)

# RUN CHAIN

# Initialize Model
adaptSteps = 5000
burnInSteps = 5000
nChains = 2
thinSteps = 1000
nPerChain = 1000000

semModel1 = jags.model("model.bug",data=semdata,n.chains=nChains, n.adapt=adaptSteps)
				
# Obtain the Posterior Sample of Factor Loadings:
parameter=c("m","A","B","C","D","E","F","a1","a2","b1","c1","c2","tau1","tau2","tau3",
            "phi","psi","xi","pi","delta") #Specify the Parameters to Be Estimated
cat("Burning in the MCMC chain ...\n")
update(semModel1, n.iter=burnInSteps)
cat("Sampling from the final MCMC chain ... \n")
codaSamples1 = coda.samples(semModel1, n.iter=nPerChain, variable.names=parameter,
                            thin=thinSteps,seed=5555)

summary(codaSamples1[[1]])    #Posterior Mean, posterior SD and posterior probablity interval (PPI) 

#Transform precision matrix of random effects to variance-covariance matrix
phi11<-codaSamples1[[1]][,40]
phi21<-codaSamples1[[1]][,41]
phi12<-codaSamples1[[1]][,42]
phi22<-codaSamples1[[1]][,43]

Sigma11 <- phi22/(phi11*phi22-phi21^2) #Sigma is covariance matrix
Sigma21 <- phi12/(phi11*phi22-phi21^2)
Sigma22 <- phi11/(phi11*phi22-phi21^2)
Sigma12 <- Sigma21

summary(Sigma11)#Posterior mean, SD and PPI of variance of random intercept
summary(Sigma22)#Posterior mean, SD and PPI of variance of random slope

plot((codaSamples1[,c(25, 6, 43, 46)]))               #Trace plots and Density plots

par(mfrow=c(2,2))
acf(codaSamples1[[1]][,25],main="MATHSCORE on MOMEDU")   #Auto-correlation plots
acf(codaSamples1[[1]][,6],main="Slope of ENJOY on NEWMETHOD") 
acf(Sigma22,main="Variance of Random Slope of ENJOY") 
acf(codaSamples1[[1]][,46],main="Error Precison of MATHSCORE") 

require(coda)
geweke.plot(codaSamples1[[1]][,c(25, 6, 43, 46)])
geweke.diag(codaSamples1[[1]][,c(25, 6, 43, 46)])

#Heidelberger-Welch diagnostics
heidel.diag(codaSamples1[[1]][,c(25, 6, 43, 46)])

#Raftery.diag
raftery.diag(codaSamples1[[1]][,c(25, 6, 43, 46)])

