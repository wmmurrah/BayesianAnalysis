#R program for Gibbs sampling from full conditionals in OLS example

#number of iterations
m=5000

#read only observations with complete information, n=2313
x=as.matrix(read.table("c:\\ols_examp.dat")[1:2313,2:10]
y=as.matrix(read.table("c:\\ols_examp.dat")[1:2313,11]

#establish parameter vectors and constant quantities
s2=matrix(1,m); b=matrix(0,m,9)

xtxi=solve(t(x)%*%x)
pars=coefficients(lm(y ~ x-1))

#Gibbs sampling begins
for(i in 2:m){
#simulate beta from its multivariate normal conditional
b[i,]=pars+t(rnorm(9,mean=0,sd=1))%*%chol(s2[i-1]*xtxi)

#simulate sigma from its inverse gamma distribution
s2[i]=1/rgamma(1,2313/2,.5*t(y-x%*%(b[i,]))%*%(y-x%*%(b[i,])))

#write output to file and screen
write(c(b[i,],s2[i]),file="c:\\ols_examp.out", append=T, ncol=10)
if(i%%50==0){print(c(i,b[i,1],s2[i]))}
}
