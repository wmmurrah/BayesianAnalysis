#---------EXAMPLE 6.4: BAYESIAN MODEL AVERAGING-------------#
install.packages("BMA")
require(BMA)
datafile <- read.csv(file.choose(),header=T)
datafile9 <- subset(datafile, select=c(rcomb1, gender, native,  slang,  ESCS,
       JOYREAD, DIVREAD, MEMOR, ELAB, CSTRAT))
attach(datafile9)
bma <- bicreg(cbind(gender, native,  slang,  ESCS,
       JOYREAD, DIVREAD, MEMOR, ELAB, CSTRAT),rcomb1,strict=FALSE,OR=20)
summary(bma)
plot(bma,include.intercept=FALSE)