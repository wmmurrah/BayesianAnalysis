
#R: sampling from marginal for variance and conditional for mean 
x=as.matrix(read.table("c:\\education.dat",header=F)[,1])
sig=rgamma(2000,(length(x)-1)/2 , rate=((length(x)-1)*var(x)/2))
sig=1/sig
mu=rnorm(2000,mean=mean(x),sd=(sqrt(sig/length(x))))







