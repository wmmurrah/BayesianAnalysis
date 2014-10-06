#R program for simulating parameters from BVN distribution

xx<-matrix(NA,2,1377)
xx[,1:361]=c(1,1); xx[,362:448]=c(1,2); xx[,449:487]=c(1,3);
xx[,488:495]=c(1,4); xx[,496:497]=c(1,5); xx[,498:499]=c(1,6);
xx[,500:608]=c(2,1); xx[,609:801]=c(2,2); xx[,802:852]=c(2,3);
xx[,853:865]=c(2,4); xx[,866:867]=c(2,5); xx[,868:870]=c(2,6);
xx[,871:915]=c(3,1); xx[,916:1006]=c(3,2); xx[,1007:1190]=c(3,3);
xx[,1191:1215]=c(3,4); xx[,1216:1219]=c(3,5); xx[,1220:1224]=c(3,6);
xx[,1225:1239]=c(4,1); xx[,1240:1256]=c(4,2); xx[,1257:1291]=c(4,3);
xx[,1292:1308]=c(4,4); xx[,1309:1312]=c(4,5); xx[,1313:1314]=c(4,6);
xx[,1315:1324]=c(5,1); xx[,1325:1328]=c(5,2); xx[,1329:1337]=c(5,3);
xx[,1338:1342]=c(5,4); xx[,1343:1344]=c(5,5); xx[,1345:1355]=c(6,1);
xx[,1356:1364]=c(6,2); xx[,1365:1368]=c(6,3); xx[,1369:1371]=c(6,4);
xx[,1372:1372]=c(6,5); xx[,1373:1377]=c(6,6);
x=t(xx)[,1]; y=t(xx)[,2]

lnpost=function(ar,amx,amy,asx,asy,ax,ay,axy)
{return(-690*log((1-ar^2)*asx*asy) 
        +(-.5/(1-ar^2))*(ax/asx - 2*ar*axy/sqrt(asx*asy) 
        + ay/asy))}

mnx=mean(x); mny=mean(y); accr=0; accx=0; accy=0
mx=matrix(0,10000); my=matrix(0,10000); 
sx=matrix(1,10000); sy=matrix(1,10000); r=matrix(0,10000)

for(i in 2:10000){
#sample mx from normal
mx[i]=rnorm(1,mean=mnx+(r[i-1]*sx[i-1]*(my[i-1]-mny))/sy[i-1]
              ,sd=sqrt(sx[i-1]*(1-r[i-1]^2)/1377))

#sample my from normal
my[i]=rnorm(1,mean=mny+(r[i-1]*sy[i-1]*(mx[i]-mnx))/sx[i-1]
              ,sd=sqrt(sy[i-1]*(1-r[i-1]^2)/1377))

#update sums of squares
sx2=sum((x-mx[i])^2); sy2=sum((y-my[i])^2); 
sxy=sum((x-mx[i])*(y-my[i]))

#sample sx
sx[i]=sx[i-1]+runif(1,min=-.1,max=.1); acc=1
if(sx[i]<0){acc=0; sx[i]=sx[i-1]}
if((lnpost(r[i-1],mx[i],my[i],sx[i],sy[i-1],sx2,sy2,sxy)
   -lnpost(r[i-1],mx[i],my[i],sx[i-1],sy[i-1],sx2,sy2,sxy))
   <log(runif(1,min=0,max=1)))
{acc=0; sx[i]=sx[i-1]}
accx=accx+acc

#sample sy 
sy[i]=sy[i-1]+runif(1,min=-.1,max=.1); acc=1
if(sy[i]<0){acc=0; sy[i]=sy[i-1]}
if((lnpost(r[i-1],mx[i],my[i],sx[i],sy[i],sx2,sy2,sxy)
   -lnpost(r[i-1],mx[i],my[i],sx[i],sy[i-1],sx2,sy2,sxy))
   <log(runif(1,min=0,max=1)))
{acc=0; sy[i]=sy[i-1]}
accy=accy+acc

#sample r from full posterior using MH step
r[i]=r[i-1]+runif(1,min=-.05,max=.05); acc=1
if(abs(r[i])>1){acc=0; r[i]=r[i-1]}
if((lnpost(r[i],mx[i],my[i],sx[i],sy[i],sx2,sy2,sxy)
   -lnpost(r[i-1],mx[i],my[i],sx[i],sy[i],sx2,sy2,sxy))
   <log(runif(1,min=0,max=1)))
{acc=0; r[i]=r[i-1]}
accr=accr+acc

if(i%%100==0){print(c(i,accr/i,accx/i,accy/i,
                    mx[i],my[i],sx[i],sy[i],r[i]),digits=4)}

}