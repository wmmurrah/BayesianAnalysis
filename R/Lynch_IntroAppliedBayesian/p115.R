#R program for Random Walk Metropolis algorithm

m=matrix(0,5000); b=matrix(.2,5000); z=matrix(0,1377)
z[1:17]=0;    z[18:32]=1;   z[33:103]=2
z[104:425]=3; z[426:826]=4; z[827:1337]=5
acctot=0

for(i in 2:5000)
 {
  m[i]=m[i-1]+rnorm(1,mean=0,sd=.002)
  b[i]=(2-25*m[i])/10
  acc=tot=1
  for(j in 1:1377)
   {tot=tot*((z[j]*m[i]+b[i])/(z[j]*m[i-1]+b[i-1]))}
  if(runif(1,min=0,max=1)>tot || b[i]<0)
   {m[i]=m[i-1]; b[i]=b[i-1]; acc=0}
  acctot=acctot+acc
 
  if(i%%10==0){print(c(i,m[i],b[i],acctot/i))}
 }

