#R program for simple random effects model

#read data
y=as.matrix(read.table("c:\\internet_examp.dat")[,3:4]

m=0; s2=10000; a=c=.01; b=d=.01; tau2=1; sigma2=1; malpha=0
n=nrow(y)

for(i in 1:20000){

#draw alpha_i
alpha= rnorm(n,
mean=(((tau2*(y[,1]+y[,2]))+sigma2*malpha)/(2*tau2+sigma2)),
sd=sqrt((tau2*sigma2)/(2*tau2+sigma2)))

#draw malpha
malpha=rnorm(1,
mean=(tau2*m+s2*sum(alpha))/((tau2+n*s2)),
sd=sqrt((tau2*s2)/((tau2+n*s2))))

#draw tau2
tau2=rgamma(1,
shape=(n/2+a),
rate=(sum((alpha-malpha)^2)+2*b)/2)

tau2=1/tau2

#draw sigma2
sigma2=rgamma(1,
shape=n+c,
rate=(sum((y-alpha)^2) +2*d)/2)
sigma2=1/sigma2

#write results to file
if(i%%10==0 | i==1){print(c(i,alpha[1],malpha,tau2,sigma2))
write(c(i,alpha[1],malpha,tau2,sigma2),file="c:\\internet_examp.out",
      append=T,ncol=5)}
}
