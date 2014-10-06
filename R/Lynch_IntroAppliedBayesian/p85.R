#R program for rejection method of sampling
count=0; k=1; f=matrix(NA,1000)
while(k<1001)
 {
  z=runif(1,min=0,max=5)
  r=((1/40)*(2*z+3))/(2*.2)
  if(r>runif(1,min=0,max=1))
   {f[k]=z; k=k+1}
  count=count+1
 }


