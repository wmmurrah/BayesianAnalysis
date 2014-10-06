#R program for MH sampling of parameters in linear regression
#number of iterations
m=20000

#read in data, establish x and y matrices 
x=as.matrix(read.table("c:\\olsexamp.dat")[1:2313,2:10])
y=as.matrix(read.table("c:\\olsexamp.dat")[1:2313,11])

#establish parameter vectors, proposal scales and acceptance rates
s2=matrix(1,m); b=matrix(0,m,9)
bscale=sqrt(diag(vcov(lm(y~x-1))))*.5
s2scale=sqrt(var(residuals(lm(y~x-1))*(2313-1)/(2313-9)))*.5
accrate=matrix(0,m,9); s2accrate=matrix(0,m)

#unnormalized posterior distribution function
post<-function(x,y,b,s2)
 {return((-1157.5*log(s2) + (-.5/s2) * (t(y- x%*%b)%*%(y- x%*%b))))}

#Begin MH Sampling
for(i in 2:m){
#temporarily set ‘new’ values of b
b[i,]=b[i-1,] 

#update regression parameters
for(j in 1:9){
#generate candidate and assume it will be accepted...
b[i,j]=b[i-1,j]+rnorm(1,mean=0, sd=bscale[j]); acc=1

#...until it is evaluated for rejection
if((post(x,y,b[i,],s2[i-1]) - post(x,y,b[i-1,],s2[i-1]))
    <log(runif(1,min=0,max=1)))
 {b[i,j]=b[i-1,j]; acc=0}	
accrate[i,j]=(accrate[i-1,j]*(i-1)+acc)/i
}

#update s2.  generate candidate and assume accepted
s2[i]=s2[i-1]+rnorm(1,mean=0, sd=s2scale); acc=1

#...until it is evaluated for rejection
if(s2[i]<0 || 
  (post(x,y,b[i,],s2[i]) - post(x,y,b[i,],s2[i-1]))
   <log(runif(1,min=0,max=1)))
 {s2[i]=s2[i-1]; acc=0}
s2accrate[i]=(s2accrate[i-1]*(i-1)+acc)/i

#write output to file and screen
write(c(b[i,],s2[i],accrate[i,],s2accrate[i]), 
      file="c:\\ols_examp.out", append=T, ncolumns=20)
if(i%%10==0){print(c(i,b[i,1],s2[i],accrate[i,1],s2accrate[i]))}
}
