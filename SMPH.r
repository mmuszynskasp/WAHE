#####################################################################################################
#################### Aim: estimate summary measures
#####################################################################################################
rm(list=ls())
dev.off()

library(eurostat)
library(MASS)
library(sjmisc)
library(tidyr)
library(dplyr)
library(fBasics)

data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC"
lt.data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\EUROSTAT"
GBD.data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\GBD"

##prevalence data
setwd(data.dir)
chronicprev <- read.table(file="chronicprev.csv", sep=",", header=TRUE)
chronicperc <- tapply(X=chronicprev$Freq, INDEX=list(chronic=chronicprev$chronic, 
                                                     age=chronicprev$agegr, country=chronicprev$country, sex=chronicprev$sex), FUN=sum)

GALIprev <- read.table(file="GALIprev.csv", sep=",", header=TRUE)
GALIperc <- tapply(X=GALIprev$Freq, INDEX=list(GALI=GALIprev$GALI, 
                                                     age=GALIprev$agegr, country=GALIprev$country, sex=GALIprev$sex), FUN=sum)

SRprev <- read.table(file="SRprev.csv", sep=",", header=TRUE)
SRperc <- tapply(X=SRprev$Freq, INDEX=list(SR=SRprev$SR, 
                                               age=SRprev$agegr, country=SRprev$country, sex=SRprev$sex), FUN=sum)

multiprev <- read.table(file="multiprev.csv", sep=",",  header=TRUE)
multiperc <- tapply(X=multiprev$Freq, INDEX=list(multi=multiprev$multi, 
                                           age=multiprev$agegr, country=multiprev$country, sex=multiprev$sex), FUN=sum)


#cname <- c( "Austria", "Belgium", "Bulgaria", "Switzerland", "Cyprus", "Czechia", "Denmark", "Estonia", "Greece", 
#             "Spain", "Finland", "France", "Croatia","Hungary","Ireland","Italy", "Lithuania", "Luxembourg", "Latvia", "Netherlands","Norway", "Poland",
#             "Portugal", "Romania", "Serbia",  "Sweden", "Slovenia", "Slovakia", "United Kingdom")

#lt data
setwd(lt.data.dir)

Lx <- read.table(file="Lx.csv", sep=",", header=TRUE)
l15 <- read.table(file="l15.csv", sep=",", header = TRUE)


Lx <- tapply(Lx$values, INDEX=list(age=Lx$agegr, country=Lx$cno, sex=Lx$sex), FUN=sum)
lx15 <-  tapply(l15$values, INDEX=list(country=l15$cno, sex=l15$sex), FUN=sum)

##########################################################################################
######## Lx redistributed by health states
##########################################################################################
Lxchronic <- array(0, dim=dim(chronicperc), dimnames=dimnames(chronicperc))
for (i in 1:2){
  Lxchronic[i,,,] <- Lx[,,]*chronicperc[i,,,]
}
#
LxGALI <- array(0, dim=dim(GALIperc), dimnames=dimnames(GALIperc))
for (i in 1:3){
  LxGALI[i,,,]  <- GALIperc[i,,,]*Lx[,,]
}
#
LxSR <- array(0, dim=dim(SRperc), dimnames=dimnames(SRperc))
for (i in 1:5){
  LxSR[i,,,]  <- SRperc[i,,,]*Lx[,,]
}
#
Lxmulti <- array(0, dim=dim(multiperc), dimnames=dimnames(multiperc))
for (i in 1:30){
  Lxmulti[i,,,]  <- multiperc[i,,,]*Lx[,,]
}


##############################################################################
###Lx reditributed with well-being weights
#weights do not differentiated between age groups
#############################################################################
setwd(data.dir)
chronicallw <- read.table(file="chronicweights.csv", sep=",", header=TRUE)
selfweights <- read.table(file="selfweights.csv", sep=",", header=TRUE)
GALIallw <- read.table(file="GALIweights.csv", sep=",", header=TRUE)
multiweights <- read.table(file="multiweights.csv", sep=",", header=TRUE)

###chronic
#prepare weights
chronicw2 <- tapply(X=chronicallw$chronic, INDEX=list(country=chronicallw$country,sex=chronicallw$sex), FUN=sum)
Lxchronicw <- Lxchronic
for (i in 1:10){
  Lxchronicw[2,i,,] <- Lxchronic[2,i,,]*chronicw2
}
 
####GALI
#prepare weights
GALIallw2 <- GALIallw %>%
  pivot_longer(cols=GALI1:GALI2,
               names_to="GALI", values_to="weights")

GALI2w <- tapply(X=GALIallw2$weights, 
                 INDEX=list(GALI=GALIallw2$GALI, country=GALIallw2$country,sex=GALIallw2$sex), FUN=sum)

LxGALIw <- LxGALI
for (j in 1:2){
  for (i in 1:10){
    LxGALIw[j+1,i,,]  <- LxGALI[j+1,i,,]*GALI2w[j,,]
  }
}

####self-rated health
#prepare weights
selfweights2 <- selfweights %>%
  pivot_longer(cols = SR2:SR5,
               names_to = "sr", values_to = "weights")

selfw2 <- tapply(X=selfweights2$weights, 
                 INDEX=list(SR=selfweights2$sr, country=selfweights2$country,sex=selfweights2$sex), FUN=sum)

LxSRw <- LxSR
for (j in 1:4){
  for (i in 1:10){
    LxSRw[j+1,i,,]  <- LxSR[j+1,i,,]*selfw2[j,,]
  }
}


####multiple
#prepare weights
multiw2 <- tapply(X=multiweights$weight, 
                 INDEX=list(multi=multiweights$level,country=multiweights$country,sex=multiweights$sex), FUN=sum)

Lxmultiw <- Lxmulti
for (j in 1:29){
  for (i in 1:10){
    Lxmultiw[j+1,i,,]  <- Lxmultiw[j+1,i,,]*multiw2[j,,]
  }
}


##########################################################################################
########## WAHE indicator
##########################################################################################
WAHEchr15 <- array(0, dim=dim(Lxchronicw[1,1,,]), dimnames=(dimnames(Lxchronicw[1,1,,])))
for (i in 1:29){
  for (j in 1:2){
    WAHEchr15[i,j] <- sum(Lxchronicw[,,i,j])/lx15[i,j]
  }
}

WAHEGALI15 <- array(0, dim=dim(LxGALIw[1,1,,]), dimnames=(dimnames(LxGALIw[1,1,,])))
for (i in 1:29){
  for (j in 1:2){
    WAHEGALI15[i,j] <- sum(LxGALIw[,,i,j])/lx15[i,j]
  }
}

WAHESR15 <- array(0, dim=dim(LxSRw[1,1,,]), dimnames=(dimnames(LxSRw[1,1,,])))
for (i in 1:29){
  for (j in 1:2){
    WAHESR15[i,j] <- sum(LxSRw[,,i,j])/lx15[i,j]
  }
}

WAHEmulti15 <- array(0, dim=dim(Lxmultiw[1,1,,]), dimnames=(dimnames(Lxmultiw[1,1,,])))
for (i in 1:29){
  for (j in 1:2){
    WAHEmulti15[i,j] <- sum(Lxmultiw[,,i,j], na.rm=TRUE)/lx15[i,j]
  }
}

###########################################################################################
###################HE indicator
###########################################################################################
HEchr15 <- array(0, dim=dim(Lxchronic[1,1,,]), dimnames=(dimnames(Lxchronic[1,1,,])))
for (i in 1:29){
  for (j in 1:2){
    HEchr15[i,j] <- sum(Lxchronic[1,,i,j])/lx15[i,j]
  }
}

HEGALI15 <- array(0, dim=dim(LxGALIw[1,1,,]), dimnames=(dimnames(LxGALIw[1,1,,])))
for (i in 1:29){
  for (j in 1:2){
    HEGALI15[i,j] <- sum(LxGALIw[1,,i,j])/lx15[i,j]
  }
}

HESR15 <- array(0, dim=dim(LxSRw[1,1,,]), dimnames=(dimnames(LxSRw[1,1,,])))
for (i in 1:29){
  for (j in 1:2){
    HESR15[i,j] <- sum(LxSRw[1,,i,j])/lx15[i,j]
  }
}

HEmulti15 <- array(0, dim=dim(Lxmultiw[1,1,,]), dimnames=(dimnames(Lxmultiw[1,1,,])))
for (i in 1:29){
  for (j in 1:2){
    HEmulti15[i,j] <- sum(Lxmultiw[1,,i,j],na.rm=TRUE)/lx15[i,j]
  }
}

########################################################################################################
################# e15
########################################################################################################
e15 <- array(0, dim=dim(lx15), dimnames=dimnames(lx15))

for (j in 1:2){
  for (i in 1:29){
    e15[i,j]  <- sum(Lx[,i,j])/lx15[i,j]
  }
}

########################################################################################################
### DALE
########################################################################################################
setwd(GBD.data.dir)
gball <- read.table(file="gbd2.csv", sep=",", header=TRUE)

###order by country as in other indicators
cname <- c( "Austria", "Belgium", "Bulgaria", "Switzerland", "Cyprus", "Czechia", "Denmark", "Estonia", "Greece", 
             "Spain", "Finland", "France", "Croatia","Hungary","Ireland","Italy", "Lithuania", "Luxembourg", "Latvia", "Netherlands","Norway", "Poland",
             "Portugal", "Romania", "Serbia",  "Sweden", "Slovenia", "Slovakia", "United Kingdom")
cno <- 1:29
gball$cno <- 0
for (i in 1:length(cno)){
  gball$cno[rownames(gball)==cname[i]] <- i
}

#order by cno
DALE <- gball[order(gball$cno),]


##########################################################################################################
#########################  write out results, also appenix data table
##########################################################################################################

males15 <- cbind(cname,e15[,1], DALE[,2], HEchr15[,1], HEGALI15[,1], HESR15[,1], HEmulti15[,1],
                 WAHEchr15[,1], WAHEGALI15[,1], WAHESR15[,1],  WAHEmulti15[,1])
fem15 <- cbind(cname,e15[,2], DALE[,1], HEchr15[,2], HEGALI15[,2], HESR15[,2], HEmulti15[,2],
                 WAHEchr15[,2], WAHEGALI15[,2], WAHESR15[,2],  WAHEmulti15[,2])

colind <- c("country","LE", "DALE", "HE.chr", "HE.GALI",  "HE.sr", "HE.multi", "WAHE.chr", "WAHE.GALI",  "WAHE.sr", "WAHE.multi")
colnames(males15) <- colind
colnames(fem15) <- colind

setwd(data.dir)
write.table(males15,file="malesmph15.csv", sep=",")
write.table(fem15, file="femsmph15.csv", sep=",")
