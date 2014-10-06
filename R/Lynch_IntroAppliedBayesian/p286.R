covmat=diag(2)
covmat[1,2]=covmat[2,1]=.5
z=matrix(0,2000,2)
q=matrix(0,2000,2)
count=0

for(i in 1:2000){
#naive simulation
z[i,]=0
while(z[i,1]<=0 | z[i,2]<=0)
{count=count+1; 
 z[i,]=rnorm(2,0,1)%*%(chol(covmat))
}

#conditional simulation based on decomposition
q[i,1]=qnorm(runif(1,min=.5,max=1),0,1)
mm=covmat[1,2]*q[i,1]
ss=1-.5^2
q[i,2]=qnorm(runif(1,min=pnorm(0,mm,sqrt(ss)),max=1),mm,sqrt(ss))
}
