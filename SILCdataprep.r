#####################################################################################################
#################### Aim: prepare data for WAHE based on 2018 SILC
#####################################################################################################

dev.off()
rm(list=ls())

########## I have downloaded original data structure that has folders by countries´ names and then folders with a year of survey 
setwd("C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC\\countries")
alldirs <- list.dirs("C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC\\countries")
myfold <- sort(unique(substr(alldirs,62,63)))[-1]

years <- 2007:2018

i=1
j=10

mypath <- paste("C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC\\countries\\",myfold[i], "\\", years[j],sep="")
setwd(mypath)

mydata1 <- read.table(file=paste("UDB_cAT",substr(years[j],3,4),"P.csv", sep=""), sep=",", header=TRUE)
mydata2 <- cbind(mydata1$PB020, mydata1$PB040, mydata1$PB100, mydata1$PB110,mydata1$PB130, mydata1$PB140,  mydata1$PB150, mydata1$PH010,mydata1$PH020, mydata1$PH030, mydata1$PW010T)

mydata1 <- read.table(file=paste("UDB_cAT",substr(years[j],3,4),"D.csv", sep=""), sep=",", header=TRUE)

mycolnames <- c("country", "weight", "monthint", "yearint", "monthbir", "yearbir", "sex", "selfrated", "chronic","GALI", "happy")
colnames(mydata2) <- mycolnames

setwd("C:\\Users\\mmuszynska\\happy\\happyindex\\data\\SILC\\all")
write.table(mydata2, file="silc2018.csv", sep=",", row.names = FALSE)

for (i in 2: length(myfold)){
  mypath <- paste("C:\\Users\\mmuszynska\\happy\\happyindex\\data\\SILC\\countries\\",myfold[i], "\\2018",sep="")
  setwd(mypath)
  myfile <- paste("UDB_c", myfold[i], "18P.csv", sep="")
  mydata1 <- read.table(file=myfile, sep=",", header=TRUE)
  mydata2 <- cbind(mydata1$PB020, mydata1$PB040, mydata1$PB100, mydata1$PB110,mydata1$PB130, mydata1$PB140,  mydata1$PB150, mydata1$PH010,mydata1$PH020, mydata1$PH030, mydata1$PW010T)
  setwd("C:\\Users\\mmuszynska\\happy\\happyindex\\data\\SILC\\all")
  write.table(mydata2, file="silc2018.csv", sep=",", row.names = FALSE, col.names=FALSE, append=TRUE)
  }

######################################################################################################
####### cosmetic changes in the variables
######################################################################################################
rm(list=ls())

setwd("C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC")
mydata <- read.table(file="silc2018.csv", sep=",", header=TRUE)

#### if year of interview is missing we assume 2018
### if month of interview is missing we assume 06
mydata$yearint[is.na(mydata$yearint)] <- 2018
mydata$monthint[is.na(mydata$monthint)] <- 6
mydata$monthbir[is.na(mydata$monthbir)] <- 6
#might be useful for tapply later
cases <- rep(1, nrow(mydata))

#age variables
age <- mydata$yearint-mydata$yearbir
age[(mydata$monthint>mydata$monthbir)] <-0
age2 <- mydata$yearint-mydata$yearbir+1
mydata <- cbind(mydata,cases,age)
mydata$age[age==0] <- age2[age==0]

agegr <- rep(0, nrow=nrow(mydata))
mydata <- cbind(mydata,agegr)
ages <- c(15,30,40,seq(50,80,5))
for (i in 1:length(ages)){
  mydata$agegr[mydata$age>=ages[i]] <- ages[i]
}
#we take only age 17+, 5-year age groups and 80+ (some error, ages 81 and 82 in the data)
mydata <- Ghappy[(Ghappy$age>16),]

## drop all missing values
chronicdata <- mydata2[(!is.na(mydata2$chronic) & !is.na(mydata2$sex)),]
mydata <- chronicdata[(!is.na(chronicdata$happy) & chronicdata$happy<99),]

#remove countries DE,MT
Ghappy <- mydata[(mydata$country!="MT"& mydata$country!="DE"),]


### replace country names with numbers
replco <- cbind(sort(unique(mydata$country)),(1:length(unique(mydata$country))))

for (i in 1:nrow(replco)){
  mydata$country[mydata$country==replco[i,1]] <- replco[i,2]
}
###add 1 to well-being variable to remove 0, it might mess up the estimations later
mydata$happy <- mydata$happy+1

write.table(mydata, file="silc2018new.csv", sep=",")
write.table(replco, file="countriesnames.csv", sep=",")





