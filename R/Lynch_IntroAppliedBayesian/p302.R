g=as.matrix(read.table("c:\\mvprob1.out"))
summ=matrix(0,31,4)

for(m in 1:31){

x=matrix(c(1,m+73,30,1,1,1,12),7)
cellsum=matrix(0,1000)

#loop over 1000 post-burnin samples of parameters
for(i in 1:1000){
s=matrix(c(1,g[i,17],g[i,17],1),2,2)
b=matrix(g[i,2:15],7,2)
t1=g[i,20]; t2=g[i,21]
xb=t(b)%*%x

#generate 1000 ppd samples to compute probabilities
#this step is equivalent to integration over desired cell
for(j in 1:1000){
zz=xb+t(rnorm(2,0,1)%*%chol(s))
if(zz[1]>t1 & zz[2]>t2){cellsum[i-201]=cellsum[i-201]+1}
}
if(i%%5==0){print(c(m,i))}
}
summ[m,1]=mean(cellsum/1000)
summ[m,2]=sd(cellsum/1000)
summ[m,3]=sort(cellsum)[25]
summ[m,4]=sort(cellsum)[975]
}
