#R: sampling from conditionals for both variance and mean

x=as.matrix(read.table("c:\\education.dat",header=F)[,1])
mu=matrix(0,2000); sig=matrix(1,2000)
for(i in 2:2000)
 {
  sig[i]=rgamma(1,(length(x)/2),rate=sum((x-mu2[i-1])^2)/2)
  sig[i]=1/sig2[i]
  mu[i]=rnorm(1,mean=mean(x),sd=(sqrt(sig2[i]/length(x))))
 }
