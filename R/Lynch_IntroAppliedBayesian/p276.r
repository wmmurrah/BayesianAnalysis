#R program for multivariate regression
x<-as.matrix(read.table("c:\\mvn_examp.dat")[,2:10])
y<-as.matrix(read.table("c:\\mvn_examp.dat")[,11:14])

d=4;k=9
b=matrix(0,(d*k)); s=diag(d)

for(i in 2:10000){
#draw b from mvn
vb=solve(solve(s)%x%(t(x)%*%x))
mn=vb%*%(as.vector(t(x)%*%y%*%t(solve(s))))

b=mn+t(rnorm((d*k),0,1)%*%chol(vb))

#draw s from inverse wishart
e=matrix((as.vector(y)-(diag(d)%x%x%*%b)),nrow(y),d)
v=t(e)%*%e

s=riwish(nrow(y)-1,v)

print(c(i,b[1],b[10],b[19],b[28],s[1,1],s[2,2],s[3,3],s[4,4]))
write(c(t(b),t(s)),file="c:\\mvn_examp.out",ncolumns=52,append=T)
}
