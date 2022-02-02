#####################################################################################################
#################### Aim: estimate well-being weights for health states for WAHE
#####################################################################################################
rm(list=ls())
dev.off()

library(MASS)
library(sjmisc)
library(tidyr)
library(dplyr)

data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC"
setwd(data.dir)

mydata <- read.table(file="silc2018new.csv", sep=",", header=TRUE)



# 
# ###data clean-up to be included later in the data preparation files 
# mydata <- read.table(file="silc2018.csv", sep=",", header=TRUE)
# #### if year of interview is missing we assume 2018
# ### if month of interview is missing we assume 06
# mydata$yearint[is.na(mydata$yearint)] <- 2018
# mydata$monthint[is.na(mydata$monthint)] <- 6
# mydata$monthbir[is.na(mydata$monthbir)] <- 6
# cases <- rep(1, nrow(mydata))
# age <- mydata$yearint-mydata$yearbir
# age[(mydata$monthint>mydata$monthbir)] <-0
# age2 <- mydata$yearint-mydata$yearbir+1
# mydata2 <- cbind(mydata,cases,age)
# mydata2$age[age==0] <- age2[age==0]
# 
# ## drop all missing values
# chronicdata <- mydata2[(!is.na(mydata2$chronic) & !is.na(mydata2$sex)),]
# mydata <- chronicdata[(!is.na(chronicdata$happy) & chronicdata$happy<99),]
# #do not remove countries DE,MT
# Ghappy <- mydata[(mydata$country!="MT"& mydata$country!="DE"),]
# 
# #we take only age 17+, 5-year age groups and 80+
# # remove DE (nobody older than 74), MT 
# mydata <- Ghappy[(Ghappy$age>16),]
# 
# agegr <- rep(0, nrow=nrow(mydata))
# mydata <- cbind(mydata,agegr)
# ages <- c(15,30,40,seq(50,80,5))
# 
# for (i in 1:length(ages)){
#   mydata$agegr[mydata$age>=ages[i]] <- ages[i] 
# }
# 
# ### replace country names with numbers
# replco <- cbind(sort(unique(mydata$country)),(1:length(unique(mydata$country))))
# 
# for (i in 1:nrow(replco)){
#   mydata$country[mydata$country==replco[i,1]] <- replco[i,2]
# }
# ###add 1 to well-being variable to remove 0, it might mess up the estimations
# mydata$happy <- mydata$happy+1
# 
# write.table(mydata, file="silc2018new.csv", sep=",")
# write.table(replco, file="countriesnames.csv", sep=",")


######## extra data preparation steps for the  ordered probit models
mydata$happy <- ordered(as.factor(mydata$happy))
mydata$country <- as.factor(as.numeric(as.character(mydata$country)))

###### help variables
sex <- sort(unique(mydata$sex))
country <- sort(unique(mydata$country))


#####################################################################################################
############## GALI
#####################################################################################################
mydata$GALI2 <- plyr::revalue(as.factor(mydata$GALI), c("1"="3", "3"="1"))
mydata$GALI2 <- as.factor(as.numeric(as.character(mydata$GALI2)))
mydata$age2 <- (mydata$age)^2/100


GALIweights <- function(thedata){
  GALImodel <- polr(happy ~ GALI2+age+age2, data= thedata, method = "probit")
  ourmodel <- coef(summary(GALImodel))

  coef0 <- ourmodel[substr(rownames(ourmodel),1,1)==1][1]
  coef9 <- ourmodel[substr(rownames(ourmodel),1,2)==10][1]
  divi <- coef9-coef0
  
  GALIw <- ourmodel[substr(rownames(ourmodel),1,1)=="G"][1:2]
  GALIweights <- cbind(i,j,t(GALIw), coef0, coef9, t(1+ GALIw/divi))
  return(GALIweights)
  }

##set the file first and then append for other country,sex,age combinations
i=1
j=1

thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
GALIw111 <- GALIweights(thedata)
colnames(GALIw111) <- c("country","sex", "coeffG1","coeffG2","thr1","thr10","GALI1", "GALI2")
write.table(GALIw111, file="GALIweights.csv", sep=",", row.names=FALSE)


j=2
thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
GALIw12k <- GALIweights(thedata)
write.table(GALIw12k, file="GALIweights.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)

for (i in 2:length(country)){
  for (j in 1:2){
      thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
      GALIwijk <- GALIweights(thedata)
      write.table(GALIwijk, file="GALIweights.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    }
  }




###############################################################################################################
##################### chronic
###############################################################################################################
mydata$chronic[mydata$chronic==2] <- 0
mydata$chronic <- as.factor(mydata$chronic)

chrweights <- function(thedata){
  chrmodel <- polr(happy ~ chronic+age+age2, data= thedata, method = "probit")
  ourmodel <- coef(summary(chrmodel))
  
  coef0 <- ourmodel[substr(rownames(ourmodel),1,1)==1][1]
  coef9 <- ourmodel[substr(rownames(ourmodel),1,2)==10][1]
  divi <- coef9-coef0
  
  chrw <- ourmodel[substr(rownames(ourmodel),1,1)=="c"][1]
  chrw <- cbind(i,j, chrw, coef0,coef9, t(1+ chrw/divi))
  return(chrw)
}

##set the file first and then append for other country,sex,age combinations
i=1
j=1

thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
chronicw <- chrweights(thedata)
colnames(chronicw) <- c("country","sex","coeff","thr1","thr10","chronic")
write.table(chronicw, file="chronicweights.csv", sep=",", row.names=FALSE)

j=2
thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
chronicw12k <- chrweights(thedata)
write.table(chronicw12k, file="chronicweights.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)



for (i in 2:length(country)){
  for (j in 1:2){
      thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
      chronicwijk <- chrweights(thedata)
      write.table(chronicwijk, file="chronicweights.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    }
  }



###################################################################################################################
############Self-rated health
###################################################################################################################
mydata$selfrated <- as.factor(mydata$selfrated)

SRweights <- function(thedata){
  SRmodel <- polr(happy ~ selfrated+age+age2, data= thedata, method = "probit")
  ourmodel <- coef(summary(SRmodel))
  
  coef0 <- ourmodel[substr(rownames(ourmodel),1,1)==1][1]
  coef9 <- ourmodel[substr(rownames(ourmodel),1,2)==10][1]
  divi <- coef9-coef0
  
  SRw <- ourmodel[substr(rownames(ourmodel),1,1)=="s"][1:4]
  SRwe <- cbind(i,j, t(SRw), coef0,coef9, t(1+ SRw/divi))
  return(SRwe)
}

##set the file first and then append for other country,sex,age combinations
i=1
j=1

thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
SRw111 <- SRweights(thedata)
colnames(SRw111) <- c("country","sex", "cSR2", "cSR3", "cSR4","cSR5", "coef1","coef10", "SR2","SR3", "SR4","SR5")
write.table(SRw111, file="selfweights.csv", sep=",", row.names=FALSE)

j=2
thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
SRw12k <- SRweights(thedata)
write.table(SRw12k, file="selfweights.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)




for (i in 2:length(country)){
  for (j in 1:2){
      thedata <- mydata[(mydata$country==country[i]& mydata$sex==sex[j]),]
      SRwijk <- SRweights(thedata)
      write.table(SRwijk, file="selfweights.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
    }
  }

  
  


####################################################################################################################
################# multidimensional
####################################################################################################################
Ghappy4 <- na.omit(mydata)
Ghappy4$multi <- as.numeric(paste(Ghappy4$chronic, Ghappy4$GALI2,Ghappy4$selfrated, sep=""))
Ghappy4$multi <- as.factor(Ghappy4$multi)


multiweights <- function(thedata){
  multimodel <- polr(happy ~ multi+age+age2, data= thedata, method = "probit")
  ourmodel <- coef(summary(multimodel))
  
  coef0 <- ourmodel[substr(rownames(ourmodel),1,1)==1][1]
  coef9 <- ourmodel[substr(rownames(ourmodel),1,2)==10][1]
  divi <- coef9-coef0
  
  multiw <- ourmodel[substr(rownames(ourmodel),6,9)!=""][1:(length(ourmodel[substr(rownames(ourmodel),6,9)!=""])/3)]
  multiwi <- 1+multiw/divi
  multiwnames <- substr(rownames(ourmodel)[substr(rownames(ourmodel),6,9)!=""],6,9)
  multiout <- cbind(i,j,multiwnames,multiw, coef0,coef9,multiwi)
  return(multiout)
}


i=1
j=1
thedata <- Ghappy4[(Ghappy4$country==country[i]& Ghappy4$sex==sex[j]),]
thedata <- droplevels(thedata)
multiw <- multiweights(thedata)

colnames(multiw) <- c("country","sex", "level", "coeff","thr1","thr10","weight")
write.table(multiw, file="multiweights.csv", sep=",", row.names=FALSE)

j=2
thedata <- Ghappy4[(Ghappy4$country==country[i]& Ghappy4$sex==sex[j]),]
thedata <- droplevels(thedata)
multiw12k <- multiweights(thedata)
write.table(multiw12k, file="multiweights.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)


for (i in 2:29){
  for (j in 1:2){
      thedata <- Ghappy4[(Ghappy4$country==country[i]& Ghappy4$sex==sex[j]),]
      thedata <- droplevels(thedata)
      multiwijk <- multiweights(thedata)
      write.table(multiwijk, file="multiweights.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}

