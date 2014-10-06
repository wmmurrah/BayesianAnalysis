ageints=9; n=5
cv=matrix(c(0,1,1,1,16),ageints,5,byrow=T)
x=matrix(1,ageints,9); x[,2]=seq(0,ageints); x[,5:9]=cv
g=as.matrix(read.table("c:\\mshaz.out"))
b=matrix(0,9,2)

radix=c(.848,.152,0)

mpower=function(mat,power)
 {ma<-diag(3);for(i in 1:power){ma=ma%*%mat};return(ma)}

for(m in 1001:2000){
#read in parameter sample
b[(1:9),1]=g[m,(2:10)]; b[(1:9),2]=g[m,(11:19)]
rho=g[m,21]
sig=matrix(c(1,rho,rho,1),2,2)

#compute predicted values for transitions: start h
x[,3]=0; x[,4]=0; hfb=x%*%b

#compute predicted values for transitions: start u
x[,3]=1; x[,4]=x[,2]; ufb=x%*%b

#establish life table variables
l=array(0,c(ageintervals,3,3)); l[1,,]=diag(3)*radix
bl=matrix(0,ageintervals,2); tl=matrix(0,ageintervals,2)

#compute transition probabilities matrices
for(a in 1:ageints){
pmat[1,3]=pmvnorm(lower=c(-Inf,-Inf),upper=c(+Inf,hfb[a,2]),mean=c(0,0),corr=sig)
pmat[2,3]=pmvnorm(lower=c(-Inf,-Inf),upper=c(+Inf,ufb[a,2]),mean=c(0,0),corr=sig)
pmat[1,2]=pmvnorm(lower=c(-Inf,hfb[a,2]),upper=c(hfb[a,1],+Inf),mean=c(0,0),corr=sig)
pmat[2,2]=pmvnorm(lower=c(-Inf,ufb[a,2]),upper=c(ufb[a,1],+Inf),mean=c(0,0),corr=sig)
pmat[1,1]=1-(pmat[1,2]+pmat[1,3])
pmat[2,1]=1-(pmat[2,2]+pmat[2,3])

#convert tp to m via Sylvester's formula
mmat=0
lam2=(pmat[2,2]+pmat[1,1]+sqrt((pmat[2,2]+pmat[1,1])^2-4*(pmat[1,1]*pmat[2,2]-pmat[1,2]*pmat[2,1])))/2
lam3=(pmat[2,2]+pmat[1,1]-sqrt((pmat[2,2]+pmat[1,1])^2-4*(pmat[1,1]*pmat[2,2]-pmat[1,2]*pmat[2,1])))/2
mmat= (log(lam2)/((lam2-1)*(lam2-lam3)))*
((pmat-diag(3))%*%(pmat-lam3*diag(3)))+(log(lam3)/((lam3-1)*(lam3-lam2)))*
((pmat-diag(3))%*%(pmat-lam2*diag(3)))

mmat=-mmat

#compute lx and Lx for next age group
if(a<ageints){
expm=diag(3);pyr=diag(3)
for(j in 1:20)
{
expm=expm + ((-1)^j)*mpower(mmat,j)/factorial(j);
pyr=pyr  +  ((-1)^j)*mpower(mmat,j)/factorial(j+1);
}

lx=l[a,,]%*%(expm)
blx=n*(l[a,,]%*%pyr)

l[a+1,1,1]=sum(lx[,1]); l[a+1,2,2]=sum(lx[,2]); l[a+1,3,3]=0;
bl[a,1]=sum(blx[,1]); bl[a,2]=sum(blx[,2])
}

if(a==ageints){
blx=l[a,1:2,1:2]%*%solve(mmat[1:2,1:2])
bl[a,1]=sum(blx[,1]); bl[a,2]=sum(blx[,2])
}
}

le=matrix(NA,ageints,2)
for(a in 1:ageints){
tl[a,1]=sum(bl[a:ageintervals,1]); tl[a,2]=sum(bl[a:ageintervals,2])
le[a,]=tl[a,]/sum(l[a,,])
}
write(c(t(le)),file="c:\\lifetab.out",append=T,ncolumns=(2*ageints))
print(c(m,le[1,1],le[1,2]))
}
print(Sys.time()-starttime)