#R program for Gibbs sampling using rejection sampling
x=matrix(-1,2000); y=matrix(-1,2000)
for(i in 2:2000)
{
 #sample from y | x using rejection sampling
 z=0
 while(z==0)
  {
   u=runif(1,min=0, max=2)
   if( ((2*u)+(3*y[i-1])+2) > (25*runif(1,min=0,max=1)*.5))
    {x[i]=u; z=1}
  }
 #sample from y | x using rejection sampling
 z=0
 while(z==0)
  {
   u=runif(1,min=0,max=2)
   if( ((2*x[i])+(3*u)+2) > (25*runif(1,min=0,max=1)*.5))
    {y[i]=u; z=1}
  }
}
