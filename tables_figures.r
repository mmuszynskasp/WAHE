#####################################################################################################
#################### Aim: produce tables and graphs from the paper
#####################################################################################################
rm(list=ls())
dev.off()

library(sjmisc)
library(tidyr)
library(dplyr)
library(fBasics)
library(corrplot)
library(ggcorrplot)
library(psych)

data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC"
setwd(data.dir)

###################################################################
#####################data for mean prevalence and weights
chrt <- read.table(file="chrprev.csv", sep=",", header=TRUE)
othert <- read.table(file="otherprev.csv",sep=",", header=TRUE)
multipre <- read.table(file="multiprevtab.csv", sep=",", header=TRUE)

weightsall <- read.table(file="allweightsaver.csv", sep=",", header=TRUE)
weightsmult <- read.table(file="bycountrymul.csv", sep=",", header=TRUE)


#####################################################################################################
################SMPH
#####################################################################################################
males15 <- read.table(file="malesmph15.csv", sep=",", header=TRUE)
fem15 <- read.table(file="femsmph15.csv", sep=",", header=TRUE)
############################################################################################
####################### Table 1 - mean and sd of SMPH across countries
malest <- rbind(colMeans(males15[,-1]), colStdevs(males15[,-1]), 100*colStdevs(males15[,-1])/colMeans(males15[,-1]))
femst <- rbind(colMeans(fem15[,-1]), colStdevs(fem15[,-1]), 100*colStdevs(fem15[,-1])/colMeans(fem15[,-1]))

bothst <- rbind(malest,femst)
write.table(bothst, file="bothst.csv", sep=",")

#############################################################################################
#################### Table A9 and A10 from supplementary material
maleorder <- males15[order(males15$LE),]
femorder <- fem15[order(fem15$LE),]

write.table(maleorder, file="TableA9.csv", sep=",", row.names = FALSE)
write.table(femorder, file="TableA10.csv", sep=",", row.names = FALSE)

##############################################################################################


############################################################################################
######### Spearman correlation coefficients plots

femcor <- cor(fem15[,-1], method="spearman")
mencor <- cor(males15[,-1], method="spearman")

###significance
femcorp <- cor_pmat(fem15[,-1],method="spearman")
mencorp <- cor_pmat(males15[,-1],method="spearman")


dev.off()
pdf(file = "C:/Users/Magdalena/demography/happyindex/figures/corplots.pdf",  width = 10, height = 5) 
par(mfrow=c(1,2))
corrplot(mencor,method="color", type = "full", p.mat = mencorp, sig.level = 0.1, cl.lim=c(0,1), 
         col=colorRampPalette(c("grey100","grey70","grey26"))(20),tl.col="black", title="Men",addCoef.col = "black", mar=c(0,0,1,0))
corrplot(femcor,method="color", type = "full", p.mat = femcorp, sig.level = 0.1, cl.lim=c(0,1), 
         col=colorRampPalette(c("grey100","grey70","grey26"))(20),tl.col="black", title="Women", addCoef.col = "black", mar=c(0,0,1,0))
dev.off()


#########################################################################################################
####Intraclass correlation coefficient
###########################################################################################################
iccfem <- alpha(fem15[,3:11])
femdrop <- c(iccfem$total$raw_alpha, iccfem$alpha.drop$raw_alpha)

iccmale <- alpha(males15[,3:11])
maledrop <- c(iccmale$total$raw_alpha,iccmale$alpha.drop$raw_alpha)

tabicc <- rbind(maledrop,femdrop)
colnames(tabicc) <- c("all",colnames(fem15[,-c(1,2)]))
rownames(tabicc)<- c("male","fem")  

setwd(data.dir)
write.table(tabicc, file="tabicc.csv", sep=",")

##################################################################################
################# Bland-Altman plots
##################################################################################
myBA <- function(x,y){
  D1 <- abs(x-y)
  A1 <- (x+y)/2
  myf <- lm(D1 ~ A1)
  myfit <- fitted(myf)
  resf <- resid(myf)
  mylow <- myfit-1.96*sd(resf)
  myhigh <- myfit+1.96*sd(resf)
  myBA <- cbind(A1,D1,myfit,mylow,myhigh)
  myBA2 <- myBA[order(myBA[,1]),]
  return(myBA2)
}

r <- function(x,y){
  D1 <- abs(x-y)
  A1 <- (x+y)/2
  r2 <- cor.test(D1,A1)
  return(r2)
}


dev.off()
pdf(file = "C:/Users/Magdalena/demography/happyindex/figures/bothDALEBA.pdf",  width = 10, height = 10) 
par(mfrow=c(3,2))

myBAm <- myBA(x=males15$DALE, y=males15$HE.chr)
myBAf <- myBA(x=fem15$DALE, y=fem15$HE.chr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE Chronic Dis.")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=1)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=1)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)
legend("topright",legend = c("Men", "Women"), pch=c(1,1), lwd=c(1,1), col=c("black", "red"))

myBAm <- myBA(x=males15$DALE, y=males15$WAHE.chr)
myBAf <- myBA(x=fem15$DALE, y=fem15$WAHE.chr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="WAHE Chronic Dis.")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$DALE, y=males15$HE.GALI)
myBAf <- myBA(x=fem15$DALE, y=fem15$HE.GALI)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE GALI")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)


myBAm <- myBA(x=males15$DALE, y=males15$WAHE.GALI)
myBAf <- myBA(x=fem15$DALE, y=fem15$WAHE.GALI)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="WAHE GALI")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$DALE, y=males15$HE.sr)
myBAf <- myBA(x=fem15$DALE, y=fem15$HE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE Self-rated")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$DALE, y=males15$LE)
myBAf <- myBA(x=fem15$DALE, y=fem15$LE)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="LE")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)


dev.off()



dev.off()
pdf(file = "C:/Users/Magdalena/demography/happyindex/figures/bothWAHESR.pdf",  width = 10, height = 10) 
par(mfrow=c(3,2))


myBAm <- myBA(x=males15$HE.chr, y=males15$WAHE.sr)
myBAf <- myBA(x=fem15$HE.chr, y=fem15$WAHE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE Chronic Dis.")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=1)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=1)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)
legend("topright",legend = c("Men", "Women"), pch=c(1,1), lwd=c(1,1), col=c("black", "red"))

myBAm <- myBA(x=males15$WAHE.chr, y=males15$WAHE.sr)
myBAf <- myBA(x=fem15$WAHE.chr, y=fem15$WAHE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="WAHE Chronic Dis.")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$HE.GALI, y=males15$WAHE.sr)
myBAf <- myBA(x=fem15$HE.GALI, y=fem15$WAHE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE GALI")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)


myBAm <- myBA(x=males15$WAHE.GALI, y=males15$WAHE.sr)
myBAf <- myBA(x=fem15$WAHE.GALI, y=fem15$WAHE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="WAHE GALI")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$HE.sr, y=males15$WAHE.sr)
myBAf <- myBA(x=fem15$HE.sr, y=fem15$WAHE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE Self-rated")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$DALE, y=males15$WAHE.sr)
myBAf <- myBA(x=fem15$DALE, y=fem15$WAHE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="DALE")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)


dev.off()




dev.off()
pdf(file = "C:/Users/Magdalena/demography/happyindex/figures/bothWAHEchr.pdf",  width = 10, height = 10) 
par(mfrow=c(3,2))


myBAm <- myBA(x=males15$WAHE.chr, y=males15$WAHE.GALI)
myBAf <- myBA(x=fem15$WAHE.chr, y=fem15$WAHE.GALI)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="WAHE GALI")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$WAHE.chr, y=males15$HE.chr)
myBAf <- myBA(x=fem15$WAHE.chr, y=fem15$HE.chr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE Chronic Dis.")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=1)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=1)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)
legend("topright",legend = c("Men", "Women"), pch=c(1,1), lwd=c(1,1), col=c("black", "red"))


myBAm <- myBA(x=males15$WAHE.chr, y=males15$WAHE.sr)
myBAf <- myBA(x=fem15$WAHE.chr, y=fem15$WAHE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="WAHE self-rated")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)


myBAm <- myBA(x=males15$WAHE.chr, y=males15$HE.GALI)
myBAf <- myBA(x=fem15$WAHE.chr, y=fem15$HE.GALI)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE GALI")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$WAHE.chr, y=males15$LE)
myBAf <- myBA(x=fem15$WAHE.chr, y=fem15$LE)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="LE")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)



myBAm <- myBA(x=males15$WAHE.chr, y=males15$HE.sr)
myBAf <- myBA(x=fem15$WAHE.chr, y=fem15$HE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE Self-rated")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)



dev.off()





dev.off()
pdf(file = "C:/Users/Magdalena/demography/happyindex/figures/bothWAHEGALI.pdf",  width = 10, height = 10) 
par(mfrow=c(3,2))


myBAm <- myBA(x=males15$WAHE.GALI, y=males15$WAHE.chr)
myBAf <- myBA(x=fem15$WAHE.GALI, y=fem15$WAHE.chr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="WAHE chronic")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$WAHE.GALI, y=males15$HE.chr)
myBAf <- myBA(x=fem15$WAHE.GALI, y=fem15$HE.chr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE Chronic Dis.")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=1)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=1)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)
legend("topright",legend = c("Men", "Women"), pch=c(1,1), lwd=c(1,1), col=c("black", "red"))


myBAm <- myBA(x=males15$WAHE.GALI, y=males15$WAHE.sr)
myBAf <- myBA(x=fem15$WAHE.GALI, y=fem15$WAHE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="WAHE self-rated")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)


myBAm <- myBA(x=males15$WAHE.GALI, y=males15$HE.GALI)
myBAf <- myBA(x=fem15$WAHE.GALI, y=fem15$HE.GALI)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE GALI")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$WAHE.GALI, y=males15$LE)
myBAf <- myBA(x=fem15$WAHE.GALI, y=fem15$LE)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="LE")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)

myBAm <- myBA(x=males15$WAHE.GALI, y=males15$HE.sr)
myBAf <- myBA(x=fem15$WAHE.GALI, y=fem15$HE.sr)
rangex <- range(cbind(myBAm[,1],myBAf[,1]))
rangey <- range(cbind(myBAm[,-1],myBAf[,-1]))
plot(x=myBAm[,1], y=myBAm[,2], xlim=rangex, ylim=rangey, xlab="Mean of SMPHs",ylab="Difference between SMPHs", main="HE Self-rated")
points(x=myBAf[,1], y=myBAf[,2], col="red")
lines(x=myBAm[,1], y=myBAm[,3], lty=1, lwd=2)
lines(x=myBAm[,1], y=myBAm[,4], lty=4, lwd=2)
lines(x=myBAm[,1], y=myBAm[,5], lty=4, lwd=2)
lines(x=myBAf[,1], y=myBAf[,3], lty=1, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,4], lty=4, col="red", lwd=2)
lines(x=myBAf[,1], y=myBAf[,5], lty=4, col="red", lwd=2)



dev.off()




















