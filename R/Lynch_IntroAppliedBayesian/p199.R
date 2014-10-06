#R program for dichotomous probit model
x=as.matrix(read.table("c:\\bookdead.dat")[,3:9])
y=as.matrix(read.table("c:\\bookdead.dat")[,10])

#create variables, set values, and write out starting values
b=matrix(0,7); vb<- solve(t(x)%*%x); ch<-chol(vb)

write(c(i,t(b)), file="c:\\dprob_gibbs.out", append=T, ncolumns=8)

#begin MCMC simulation
for(i in 2:25000){

#simulate latent data from truncated normal distributions
u=as.matrix(runif(length(y),min=0,max=1))
xb=as.matrix(x%*%b)

ystar=qnorm(y*u + u*(-1)^y*pnorm(0,mean=xb,sd=x[,1]) + 
            y*pnorm(0,mean=xb,sd=1), mean=xb, sd=x[,1])

#simulate beta vector from appropriate mvn
b=vb%*%(t(x)%*%ystar) + t((rnorm(7,mean=0,sd=1))%*%ch)

write(c(i,t(b)), file="c:\\dprob_gibbs.out", append=T, ncolumns=8)
if(i%%10==0){print(c(i,t(b))),digits=2)}
}



