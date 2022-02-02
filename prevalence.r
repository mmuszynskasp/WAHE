#####################################################################################################
#################### Aim: estimate prevalence of health states from SILC
#####################################################################################################
rm(list=ls())
dev.off()

library(eurostat)
library(MASS)
library(sjmisc)
library(tidyr)
library(dplyr)
data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC"
setwd(data.dir)

mydata<- read.table(file="silc2018new.csv", sep=",", header=TRUE)

cname <- c( "Austria", "Belgium", "Bulgaria", "Switzerland", "Cyprus", "Czechia", "Denmark", "Estonia", "Greece", 
             "Spain", "Finland", "France", "Croatia","Hungary","Ireland","Italy", "Lithuania", "Luxembourg", "Latvia", "Netherlands","Norway", "Poland",
             "Portugal", "Romania", "Serbia",  "Sweden", "Slovenia", "Slovakia", "United Kingdom")

#countriesnames <- read.table(file="countriesnames.csv", sep=",")

###help variable for large age-groups
agegr <- rep(0, nrow=nrow(mydata))
mydata <- cbind(mydata,agegr)
ages <- c(15,30,40,seq(50,80,5))

for (i in 1:length(ages)){
  mydata$agegr[mydata$age>=ages[i]] <- ages[i] 
}

#########recode as in the weights
mydata$chronic[mydata$chronic==2] <- 0
mydata$country <- as.factor(as.numeric(as.character(mydata$country)))
mydata$GALI2 <- revalue(as.factor(mydata$GALI), c("1"="3", "3"="1"))
mydata$GALI2 <- as.factor(as.numeric(as.character(mydata$GALI2)))


#################################################################################################################
#######chronic
#################################################################################################################
chronic <- tapply(X=mydata$weight, INDEX=list(chronic=mydata$chronic, agegr=mydata$agegr,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)
chronic[is.na(chronic)] <- 0
pop<- chronic[1,,,]+chronic[2,,,]


#correct ages 15, 16 as half 17
chronic17 <- tapply(X=mydata$weight, INDEX=list(chronic=mydata$chronic, agegr=mydata$age,country=mydata$country, sex=mydata$sex), FUN=sum,na.rm=TRUE)[,1,,]
chronic17[is.na(chronic17)] <- 0
pop17 <- tapply(X=mydata$weight, INDEX=list(agegr=mydata$age,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)[1,,]

chronic[2,1,,] <- chronic[2,1,,]+chronic17[2,,]
chronic[1,1,,] <- chronic[1,1,,]-chronic17[2,,]+2*pop17
pop[1,,] <- pop[1,,] + 2*pop17
#####
chronicperc <- array(0, dim=dim(chronic), dimnames=dimnames(chronic))
for (i in 1:2){
  chronicperc[i,,,] <- chronic[i,,,]/pop
}
chronicperc[is.na(chronicperc)] <- 0

chronicprev <- as.data.frame.table(chronicperc)
write.table(chronicprev, file="chronicprev.csv", sep=",")



###########################################################################################
#GALI
###########################################################################################
GALI <- tapply(X=mydata$weight, INDEX=list(GALI=mydata$GALI2, agegr=mydata$agegr,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)
GALI[is.na(GALI)]<- 0
pop <- GALI[1,,,]+GALI[2,,,]+GALI[3,,,]

#correct ages 15, 16 as half 17
GALI17 <- tapply(X=mydata$weight, INDEX=list(GALI=mydata$GALI2, agegr=mydata$age,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)[,1,,]
GALI17[is.na(GALI17)] <- 0
pop17 <- tapply(X=mydata$weight, INDEX=list(agegr=mydata$age,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)[1,,]

GALI[2,1,,] <- GALI[2,1,,]+GALI17[2,,]
GALI[3,1,,] <- GALI[3,1,,]+GALI17[3,,]
GALI[1,1,,] <- GALI[1,1,,]-GALI17[2,,]-GALI17[3,,]+ 2*pop17


pop[1,,] <- pop[1,,] + 2*pop17
#####
GALIperc <- array(0, dim=dim(GALI), dimnames=dimnames(GALI))
for (i in 1:3){
  GALIperc[i,,,] <- GALI[i,,,]/pop
}

GALIperc[is.na(GALIperc)] <- 0

GALIprev <- as.data.frame.table(GALIperc)
write.table(GALIprev, file="GALIprev.csv", sep=",")


##############################################################################################################
#SR
##############################################################################################################

SR <- tapply(X=mydata$weight, INDEX=list(SR=mydata$selfrated, agegr=mydata$agegr,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)
SR[is.na(SR)] <- 0
pop<- SR[1,,,]+SR[2,,,]+SR[3,,,]+SR[4,,,]+SR[5,,,]

#correct ages 15, 16 as half 17
SR17 <- tapply(X=mydata$weight, INDEX=list(SR=mydata$selfrated, agegr=mydata$age,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)[,1,,]
SR17[is.na(SR17)] <- 0
pop17 <- tapply(X=mydata$weight, INDEX=list(agegr=mydata$age,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)[1,,]

SR[2:5,1,,] <- SR[2:5,1,,]+SR17[2:5,,]
SR[1,1,,] <- SR[1,1,,]-SR17[2,,]-SR17[3,,]-SR17[4,,]-SR17[5,,]+2*pop17
pop[1,,] <- pop[1,,] + 2*pop17
#####
SRperc <- array(0, dim=dim(SR), dimnames=dimnames(SR))
for (i in 1:5){
  SRperc[i,,,] <- SR[i,,,]/pop
}
SRperc[is.na(SRperc)] <- 0


SRprev <- as.data.frame.table(SRperc)
write.table(SRprev, file="SRprev.csv", sep=",")



#####################################################################################################################
#multi
####################################################################################################################
mydata$multi <- as.numeric(paste(mydata$chronic, mydata$GALI2, mydata$selfrated, sep=""))


multi <- tapply(X=mydata$weight, INDEX=list(multi=mydata$multi, agegr=mydata$agegr,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)
pop <- tapply(X=mydata$weight, INDEX=list(agegr=mydata$agegr,country=mydata$country, sex=mydata$sex), na.rm=TRUE, FUN=sum)

#correct ages 15, 16 as half 17
multi17 <- tapply(X=mydata$weight, INDEX=list(multi=mydata$multi, agegr=mydata$age,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)[,1,,]
multi17[is.na(multi17)] <- 0
pop17 <- tapply(X=mydata$weight, INDEX=list(agegr=mydata$age,country=mydata$country, sex=mydata$sex), FUN=sum, na.rm=TRUE)[1,,]

multi[2:27,1,,] <- multi[2:27,1,,]+multi17[2:27,,]
multi[1,1,,] <- multi[1,1,,]-colSums(multi17[2:27,,])+2*pop17


pop[1,,] <- pop[1,,] + 2*pop17
#####
multiperc <- array(0, dim=dim(multi), dimnames=dimnames(multi))
for (i in 1:30){
  multiperc[i,,,] <- multi[i,,,]/pop
}
multiperc[is.na(multiperc)] <- 0

multiprev <- as.data.frame.table(multiperc)
write.table(multiprev, file="multiprev.csv", sep=",")

