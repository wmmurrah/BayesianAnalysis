x=as.matrix(read.table("c:\\bookheal.dat")[,1:7])
y=as.matrix(read.table("c:\\bookheal.dat")[,8])

t=matrix(0,6)
t[1]=-Inf; t[2]<-0; t[6]<-Inf
t[3]=qnorm(sum(y<=2)/nrow(y),-qnorm(sum(y==1)/nrow(y),0,1),1)
t[4]=qnorm(sum(y<=3)/nrow(y),-qnorm(sum(y==1)/nrow(y),0,1),1)
t[5]=qnorm(sum(y<=4)/nrow(y),-qnorm(sum(y==1)/nrow(y),0,1),1)

b=matrix(0,7); vb=solve(t(x)%*%x); ch=chol(vb)

write(c(1,0,0,0,0,0,0,0,0,0,0),file="c:\\oprob_gibbs.out", append=T, ncolumns=11)

for(i in 2:100000){
#simulate latent data from truncated normal distributions
xb=as.matrix(x%*%b)
ystar=qnorm(runif(nrow(y),min=pnorm(t[y],xb,1),
                          max=pnorm(t[y+1],xb,1)),mean=xb,1)

#simulate thresholds
for(k in 3:5){t[k]=runif(1,min=max(ystar[y==k-1]), max=min(ystar[y==k]))}

#simulate beta vector from appropriate mvn
b<-vb%*%(t(x)%*%ystar) + t((rnorm(7,mean=0,sd=1))%*%ch)

write(c(i,t(b),t[3],t[4],t[5]), file="c:\\oprob_gibbs.out", 
      append=T, ncolumns=11)
if(i%%10==0){print(c(i,t(b),t[3],t[4],t[5]),digits=2)}
}



