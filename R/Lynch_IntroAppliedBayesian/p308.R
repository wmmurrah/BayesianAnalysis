#read data
x=as.matrix(read.table("c:\\mvnprob.dat2")[,1:9])
z=as.matrix(read.table("c:\\mvnprob.dat2")[,10:11])

#establish variables
zstar=matrix(0,nrow(z),2)
b=matrix(0,(18)) 
s=diag(2); cs=diag(2)

acctot=0

#define thresholds--note 'trick' for t3 and t4
tz=matrix(0,4);ctz=matrix(0,4)
tz[1]=-Inf; tz[2]=0; tz[3]=tz[4]=Inf

write(c(0,t(b),t(s)),file="c:\\mshaz.out",ncol=23,append=T)

#begin Gibbs sampler
for(i in 2:10000){

#draw latent data
bb=matrix(b,9,2)
m=x%*%bb

#mini mvn gibbs sampler
mm=m[,2] + s[1,2]*(zstar[,1]-m[,1])
ss=1-s[1,2]^2

zstar[,2]=qnorm(runif(nrow(z),
min=pnorm(tz[z[,2]+1],mm,sqrt(ss)),
max=pnorm(tz[z[,2]+2],mm,sqrt(ss))),mean=mm,sd=sqrt(ss))

mm=m[,1] + s[1,2]*(zstar[,2]-m[,2])
ss=1-s[1,2]^2

zstar[,1]=qnorm(runif(nrow(z),
min=pnorm(tz[z[,1]-z[,2]+1],mm,sqrt(ss)),
max=pnorm(tz[z[,1]+z[,2]+2],mm,sqrt(ss))),mean=mm,sd=sqrt(ss))

#draw b from mvn
vb=solve(solve(s)%x%(t(x)%*%x))
mn=vb%*%(as.vector(t(x)%*%zstar%*%t(solve(s))))
b=mn+t(rnorm((d*k),0,1)%*%chol(vb))

#simulate s using MH sampling
e=matrix((as.vector(zstar)-(diag(d)%x%x%*%b)),nrow(z),d)
v=t(e)%*%e

like=-.5*(d+nrow(z)+1)*log(det(s))-.5*sum(diag(v%*%solve(s)))
cs[1,2]=cs[2,1]=s[1,2]+rnorm(1,mean=0,sd=.03)

if(abs(cs[1,2])<1){
cslike=-.5*(d+nrow(z)+1)*log(det(cs))-.5*sum(diag(v%*%solve(cs)))

if((cslike-like)>log(runif(1,0,1)))
{s[1,2]=s[2,1]=cs[1,2]; acctot=acctot+1} 
}

if(i%%10==0){print(c(i,b[1],b[1+k],s[1,2],acctot/i),digits=5)}
if(i%%5==0){write(c(i,t(b),t(s)),file="c:\\mshaz.out",ncol=23,append=T)}
}
print(Sys.time()-starttime)

