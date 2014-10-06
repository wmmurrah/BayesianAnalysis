d=5
x=matrix(NA,10000,d)

mu=matrix(c(1,2,3,4,5),d)
sig=matrix(.5,d,d)
sig[1,1]=5; sig[2,2]=4; sig[3,3]=3; sig[4,4]=2; sig[5,5]=1

for(i in 1:10000){
#simulate 5th element
x[i,d]=rnorm(1,mu[d],sqrt(sig[d,d]))

#simulate 4th, 3rd ...
for(j in d:2)
 {
  m=mu[j-1] + 
    sig[(j-1),(j:d)]%*%solve(sig[j:d,j:d])%*%(x[i,j:d]-mu[j:d])
  s=sig[j-1,j-1]-
    sig[(j-1),(j:d)]%*%(solve(sig[j:d,j:d]))%*%sig[(j:d),(j-1)]

  x[i,j-1]=rnorm(1,m,sqrt(s))
 }
print(c(i,x[i,]))
}
