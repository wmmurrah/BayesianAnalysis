#R program for Gibbs sampling from a bivariate normal pdf
x=matrix(-10,2000); y=matrix(-10,2000)
for(j in 2:2000)
 {
  #sampling from x|y
  x[j]=rnorm(1,mean=(.5*y[j-1]),sd=sqrt(1-.5*.5))
  #sampling from y|x
  y[j]=rnorm(1,mean=(.5*x[j]),sd=sqrt(1-.5*.5))
 }
