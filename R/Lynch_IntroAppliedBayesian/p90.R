#R program for Gibbs sampling using inversion method
x=matrix(-5,2000); y=matrix(-5,2000)
for(i in 2:2000)
{
 #sample from x | y
 u=runif(1,min=0,max=1)
 x[i]=sqrt(u*(6*y[i-1]+8)+(1.5*y[i-1]+1)*(1.5*y[i-1]+1))-(1.5*y[i-1]+1)
 
 #sample from y | x
 u=runif(1,min=0,max=1)
 y[i]=sqrt((2*u*(4*x[i]+10))/3 + ((2*x[i]+2)/3)*((2*x[i]+2)/3))-((2*x[i]+2)/3)
}


