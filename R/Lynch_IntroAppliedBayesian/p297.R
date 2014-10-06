#R program for multivariate probit model

x=as.matrix(read.table("c:\\mvnprob.dat1")[,1:7])
z=as.matrix(read.table("c:\\mvnprob.dat1")[,8:9])

#create variables and starting values
zstar=matrix(0,nrow(z),2)
d=2;k=7
b=matrix(0,(d*k)) 
s=cs=diag(d)

tz=matrix(0,d,4);ctz=matrix(0,d,4)
tz[,1]=-Inf; tz[,2]=0; tz[,4]=Inf
tz[1,3]=qnorm(sum(z[,1]<=2)/nrow(z), mean=-qnorm(sum(z[,1]==1)/nrow(z),mean=0,sd=1), sd=1)
tz[2,3]=qnorm(sum(z[,2]<=2)/nrow(z), mean=-qnorm(sum(z[,2]==1)/nrow(z),mean=0,sd=1), sd=1)

ctz=tz

acc1=acc2=acctot=0

write(c(0,t(b),t(s),tz[1,3],tz[2,3]),file="c:\\mvnprob.out",ncol=(d*k+k*k +3),append=T)

#begin Gibbs sampling
for(i in 2:6000){

#draw latent data: one-iteration gibbs sampler for tmvn simulation
bb=matrix(b,k,2)
m=x%*%bb

for(j in 1:d)
{ 
mm=m[,j] + t(s[j,-j])%*%solve(s[-j,-j])%*%(zstar[,-j]-m[,-j])
ss=s[j,j]-t(s[j,-j])%*%(solve(s[-j,-j]))%*%s[j,-j]

zstar[,j]=qnorm(runif(nrow(z),
min=pnorm(tz[j,z[,j]],mm,sqrt(ss)),
max=pnorm(tz[j,(z[,j]+1)],mm,sqrt(ss))),mean=mm,sd=sqrt(ss))
}

#draw thresholds using Cowles' algorithm

ctz[1,3]=qnorm(runif(1,min=pnorm(0,mean=tz[1,3],sd=.01),max=1),tz[1,3],sd=.01)
r=as.matrix((pnorm(ctz[1,z[,1]+1]-m[,1],0,1)-pnorm(ctz[1,z[,1]]-m[,1],0,1))/(pnorm(tz[1,z[,1]+1]-m[,1],0,1)-pnorm(tz[1,z[,1]]-m[,1],0,1)))
r=t(log(r))%*%matrix(1,nrow(z))+log((1-pnorm(-tz[1,3]/.01,0,1))/(1-pnorm(-ctz[1,3]/.01,0,1)))
if(r>log(runif(1,0,1))){tz[1,3]=ctz[1,3]; acc1=acc1+1}

ctz[2,3]=qnorm(runif(1,min=pnorm(0,mean=tz[2,3],sd=.01),max=1),tz[2,3],sd=.01)
r=as.matrix((pnorm(ctz[2,z[,2]+1]-m[,2],0,1)-pnorm(ctz[2,z[,2]]-m[,2],0,1))/(pnorm(tz[2,z[,2]+1]-m[,2],0,1)-pnorm(tz[2,z[,2]]-m[,2],0,1)))
r=t(log(r))%*%matrix(1,nrow(z))+log((1-pnorm(-tz[2,3]/.01,0,1))/(1-pnorm(-ctz[2,3]/.01,0,1)))
if(r>log(runif(1,0,1))){tz[2,3]=ctz[2,3]; acc2=acc2+1}


#draw b from mvn
vb=solve(solve(s)%x%(t(x)%*%x))
mn=vb%*%(as.vector(t(x)%*%zstar%*%t(solve(s))))
b=mn+t(rnorm((d*k),0,1)%*%chol(vb))

e=matrix((as.vector(zstar)-(diag(d)%x%x%*%b)),nrow(z),d)
v=t(e)%*%e

like=-.5*(d+nrow(z)+1)*log(det(s))-.5*sum(diag(v%*%solve(s)))
cs[1,2]=cs[2,1]=s[1,2]+rnorm(1,mean=0,sd=.01)

if(abs(cs[1,2])<1){
cslike=-.5*(d+nrow(z)+1)*log(det(cs))-.5*sum(diag(v%*%solve(cs)))

if((cslike-like)>log(runif(1,0,1)))
{s[1,2]=s[2,1]=cs[1,2]; acctot=acctot+1} 
}

if(i%%10==0){print(i)}
if(i%%5==0){
write(c(i,t(b),t(s),tz[1,3],tz[2,3]),file="c:\\mvnprob.out",ncol=(d*k+k*k+3),append=T)}
}


