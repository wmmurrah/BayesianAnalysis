#R program for Gibbs sampling using composition method in OLS

#number of iterations
m=20000

x=as.matrix(read.table("c:\\ols_examp.dat")[1:2313,2:10]
y=as.matrix(read.table("c:\\ols_examp.dat")[1:2313,11]

#establish parameter vectors and constant quantities
s2=matrix(1,m); b=matrix(0,m,9)
xtxi=solve(t(x)%*%x)
pars=coefficients(lm(y[,1] ~ x-1))

#simulate sigma from its inverse gamma marginal
s2=1/rgamma(m,(2313-9)/2,.5*t(residuals(lm(y[,1] ~ x-1)))%*%residuals(lm(y[,1] ~ x-1)))

#simulate beta vector from appropriate mvn
for(i in 1:m)
{
b[i,]=pars+t(rnorm(9,mean=0,sd=1))%*%chol(s2[i]*xtxi)

#write output to file and screen
write(c(b[i,],s2[i]), 
      file="c:\\ols_examp.out", append=T, ncolumns=10)
if(i%%50==0){print(c(i,b[i,1],s2[i]))}
}
