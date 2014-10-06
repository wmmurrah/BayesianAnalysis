a=matrix(10,100000);b=matrix(10,100000); acca=0; accb=0
y=matrix(c(556,346,312,284),4); n=matrix(c(1067,685,637,628),4)
k=matrix((y)/n,100000,4,byrow=T)

apost=function(f,g,k)
{
post=4*(lgamma(f+g)-lgamma(f)-lgamma(g)) + f * sum(log(k))
post=post+(1-1)*log(f)-(f*.1)
return(post)
}

bpost=function(f,g,k)
{
post=4*(lgamma(f+g)-lgamma(f)-lgamma(g)) + g * sum(log(1-k))
post=post+(1-1)*log(g)-(g*.1)
return(post)
}

for(i in 2:100000){

#draw a
a[i]=a[i-1]+rnorm(1,0,20)
if(a[i]>0)
{
acca=acca+1
newpost=apost(a[i],b[i-1],k[i-1,])
oldpost=apost(a[i-1],b[i-1],k[i-1,])

if(log(runif(1,min=0,max=1))>(newpost-oldpost))
 {a[i]=a[i-1]; acca=acca-1}
}
if(a[i]<0){a[i]=a[i-1]}

#draw b
b[i]=b[i-1]+rnorm(1,0,20)
if(b[i]>0)
{
accb=accb+1
newpost=bpost(a[i],b[i],k[i-1,])
oldpost=bpost(a[i],b[i-1],k[i-1,])

if(log(runif(1,min=0,max=1))>(newpost-oldpost))
 {b[i]=b[i-1]; accb=accb-1}
}
if(b[i]<0){b[i]=b[i-1]}

#draw k from beta
k[i,]=rbeta(4,(y+a[i]),(n-y+b[i]))

if(i%%10==0){print(c(i,a[i],b[i],acca/i,accb/i))}
}
