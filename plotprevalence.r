#####################################################################################################
#################### Aim: plot the scatter plots of the age-standardized prevalence of decreased health
#####################################################################################################

rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Rcpp)
library(ggpubr)
library(eurostat)

data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC"
setwd(data.dir)

chronicpre <- read.table(file="chronicprev.csv", sep=",", header=TRUE)
selfpre <- read.table(file="SRprev.csv", sep=",", header=TRUE)
GALIpre <- read.table(file="GALIprev.csv", sep=",", header=TRUE)
multipre <- read.table(file="multiprev.csv", sep=",", header=TRUE)

######age standardize
###read-in total EU population from Eurosat, by sex and age, 28 EU countries, 2018
pop <- get_eurostat("demo_pjan", time_format = "num") 


standpop <- pop %>%
  dplyr::filter(geo == "EU28", time==2018, sex=="T")%>%
  mutate(age=replace(age,age=="Y_GE85", "Y85")) %>%
  mutate(age2=as.numeric(substr(age,2,3))) %>%
  mutate(agegr=0)%>%
  dplyr::filter(age2>=15) %>%
  select(age2,agegr,values)

#replace age groups as in prevalence
ages <- c(15,30,40,seq(50,80,5))

for (i in 1:length(ages)){
  standpop$agegr[standpop$age2>=ages[i]] <- ages[i] 
}

#standard population
stand <- standpop %>%
  group_by(agegr)%>%
  summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt))

###############################################################################################
###### standardize prevalence by country, sex and each dimension
##############################################################################################
chronicpre <- chronicpre[chronicpre$chronic==1,]
chronic <- tapply(X=chronicpre$Freq, INDEX=list(country=chronicpre$country,sex=chronicpre$sex, agegr=chronicpre$agegr), FUN=sum)
chronicst <- chronic[,,1]
for (i in 1:29){
  for (j in 1:2){
    chronicst[i,j] <- sum(chronic[i,j,]*stand$freq)
  }
}
chronicst <- as.data.frame.table(chronicst)

GALI <- tapply(X=GALIpre$Freq, INDEX=list(country=GALIpre$country,sex=GALIpre$sex, agegr=GALIpre$agegr, GALI=GALIpre$GALI), FUN=sum)
GALIst <- GALI[,,1,-1]
for(k in 1:2){
  for (i in 1:29){
    for (j in 1:2){
      GALIst[i,j,k] <- sum(GALI[i,j,,k+1]*stand$freq)
    }
  }
}


self <- tapply(X=selfpre$Freq, INDEX=list(country=selfpre$country,sex=selfpre$sex, agegr=selfpre$agegr, self=selfpre$SR), FUN=sum)
selfst <- self[,,1,]
for(k in 1:4){
  for (i in 1:29){
    for (j in 1:2){
      selfst[i,j,k] <- sum(self[i,j,,k]*stand$freq)
    }
  }
}



multipre2 <- multipre[(multipre$multi==11 |multipre$multi==12 | multipre$multi==123),]
multi <- tapply(X=multipre2$Freq, INDEX=list(country=multipre2$country,sex=multipre2$sex, agegr=multipre2$agegr, multi=multipre2$multi), FUN=sum)
multist <- multi[,,1,]

for(k in 1:3){
  for (i in 1:29){
    for (j in 1:2){
      multist[i,j,k] <- sum(multi[i,j,,k]*stand$freq)
    }
  }
}

multim <- as.data.frame(multist[,1,])
multif <- as.data.frame(multist[,2,])
bycountrymm <- multim %>%
  pivot_longer(cols=2:3, names_to="health", values_to="prev")
bycountrymf <- multif %>%
  pivot_longer(cols=2:3, names_to="health", values_to="prev")

####################################
##for summary tables only
tabmulti <- c("11","12","13","112","113","123")
multipre2 <- multipre[multipre$multi %in%tabmulti,]
multi <- tapply(X=multipre2$Freq, INDEX=list(country=multipre2$country,sex=multipre2$sex, agegr=multipre2$agegr, multi=multipre2$multi), FUN=sum)
multist <- multi[,,1,]

for(k in 1:4){
  for (i in 1:29){
    for (j in 1:2){
      multist[i,j,k] <- sum(multi[i,j,,k]*stand$freq)
    }
  }
}

multim <- as.data.frame(multist[,1,])
multif <- as.data.frame(multist[,2,])

bystatem <- multim %>%
  pivot_longer(cols=1:6, names_to="health", values_to="prev")%>%
  group_by(health) %>%
  summarise(avg=mean(prev))
bystatef <- multif %>%
  pivot_longer(cols=1:6, names_to="health", values_to="prev")%>%
  group_by(health) %>%
  summarise(avg=mean(prev))

tablemultiprev <- cbind(bystatem,bystatef)
write.table(tablemultiprev,file="multiprevtab.csv",sep=",")

##########countinue for the plot
bycountrym <- cbind(chronicst[chronicst$sex==1,],GALIst[,1,],selfst[,1,])
bycountryf <- cbind(chronicst[chronicst$sex==2,],GALIst[,2,],selfst[,2,])
namescol <- c("country","sex","chr","GALI1","GALI2", paste("level",1:5,sep=""))
colnames(bycountryf) <- namescol
colnames(bycountrym) <- namescol

bycountrym <- bycountrym %>%
  pivot_longer(cols=c(4:5,7:10), names_to="health", values_to="prev")

bycountryf <- bycountryf %>%
  pivot_longer(cols=c(4:5,7:10), names_to="health", values_to="prev")


####write out tables
chrt <- chronicst %>%
  group_by(sex) %>%
  summarise(avg = mean(Freq))

bycountry <- rbind(bycountryf,bycountrym)
othert <- bycountry %>%
  group_by(sex,health) %>%
  summarise(avg = mean(prev))

write.table(chrt, file="chrprev.csv", sep=",")
write.table(othert, file="otherprev.csv",sep=",")

##############################################################################################
######## plots
##############################################################################################
##male
p1 <- ggplot(subset(bycountrym, (health=="GALI1"| health=="GALI2")) , aes(chr, prev,col=health, label=country))  
p2 <- ggplot(subset(bycountrym, (health=="level2" | health=="level3")) , aes(level1, prev,col=health, label=country))  
p3 <- ggplot(subset(bycountrym, (health=="level4" | health=="level5")) , aes(chr, prev,col=health, label=country))  
p4 <- ggplot(bycountrymm, aes(`11`, prev,col=health))  



dev.off()
pdf(file = "C:/Users/Magdalena/demography/happyindex/figures/prevalencemale.pdf", width = 10, height = 9)

ggarrange(
  p2 + geom_point()+ 
      theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
  #   legend.title = element_blank(),
      legend.key=element_blank(), legend.position = "top",  legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
      scale_color_manual("Self-rated Health",values=c("black", "grey70"), labels=c("Very Good","Good"))+
      xlab("Excellent Self-rated Health Prevalence")+ xlim(0,0.42)+
      scale_y_continuous(name="Self-rated Health Prevalence", breaks=seq(0,0.6,0.2), labels=waiver(), limits=c(0,0.6))+
      stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.7, label.y.npc=1),
  
  p3 + geom_point()+ #cor1=0.47**, cor2=0.44**
      theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(),legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
      scale_color_manual("Self-rated Health", values=c("black", "grey70"), labels=c("Fair", "Poor"))+
      xlab("Chronic Diseases Prevalence")+ xlim(0.1,0.6)+
      scale_y_continuous(name="Self-rated Health Prevalence", breaks=seq(0,0.6,0.2), labels=waiver(), limits=c(0,0.6))+
      stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.5, label.y.npc=0.95),
  
  p1 + geom_point()+ #cor1=0.27, cor2=0.5***
      theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(), legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
    scale_color_manual("GALI", values=c("black", "grey70"), labels=c("Moderately Limited", "Severely Limited"))+
    xlab("Chronic Diseases Prevalence")+ xlim(0.1,0.6)+
    scale_y_continuous(name="GALI Prevalence", breaks=seq(0,0.6,0.2), labels=waiver(), limits=c(0,0.4))  +
    stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.1, label.y.npc=1),
  
  p4 + geom_point()+ # cor=0.15, cor=0.48**
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(), legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
    scale_color_manual("Dimensions' Levels:", values=c("black", "grey70"), labels=c("No,Not,Very Good", "Yes,Moderately,Good"))+
    xlab("No chronic, No limitations, Excellent Health Prevalence")+ 
    ylab("Multiple Dimensions Prevalence")+
    stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.7, label.y.npc=1),
  ncol = 2, nrow = 2)



dev.off()

###female
p1 <- ggplot(subset(bycountryf, (health=="GALI1"| health=="GALI2")) , aes(chr, prev,col=health, label=country))  
p2 <- ggplot(subset(bycountryf, (health=="level2" | health=="level3")) , aes(level1, prev,col=health, label=country))  
p3 <- ggplot(subset(bycountryf, (health=="level3" | health=="level4")) , aes(chr, prev,col=health, label=country))  
p4 <- ggplot(bycountrymf, aes(`11`, prev,col=health))  



dev.off()
pdf(file = "C:/Users/Magdalena/demography/happyindex/figures/prevalencefem.pdf", width = 10, height = 9)

ggarrange(
  p2 + geom_point()+ 
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
      #   legend.title = element_blank(),
      legend.key=element_blank(), legend.position = "top",  legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
    scale_color_manual("Self-rated Health",values=c("black", "grey70"), labels=c("Very Good","Good"))+
    xlab("Excellent Self-rated Health Prevalence")+ xlim(0,0.42)+
    scale_y_continuous(name="Self-rated Health Prevalence", breaks=seq(0,0.6,0.2), labels=waiver(), limits=c(0,0.6))+
    stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.7, label.y.npc=1),
  
  p3 + geom_point()+ #cor1=0.47**, cor2=0.44**
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(),legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
    scale_color_manual("Self-rated Health", values=c("black", "grey70"), labels=c("Fair", "Poor"))+
    xlab("Chronic Diseases Prevalence")+ xlim(0.1,0.6)+
    scale_y_continuous(name="Self-rated Health Prevalence", breaks=seq(0,0.6,0.2), labels=waiver(), limits=c(0,0.6))+
    stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.5, label.y.npc=0.95),
  
  p1 + geom_point()+ #cor1=0.27, cor2=0.5***
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(), legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
    scale_color_manual("GALI", values=c("black", "grey70"), labels=c("Moderately Limited", "Severely Limited"))+
    xlab("Chronic Diseases Prevalence")+ xlim(0.1,0.6)+
    scale_y_continuous(name="GALI Prevalence", breaks=seq(0,0.6,0.2), labels=waiver(), limits=c(0,0.4))  +
    stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.1, label.y.npc=1),
  
  p4 + geom_point()+ # cor=0.15, cor=0.48**
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(), legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
    scale_color_manual("Dimensions' Levels:", values=c("black", "grey70"), labels=c("No,Not,Very Good", "Yes,Moderately,Good"))+
    xlab("No chronic, No limitations, Excellent Health Prevalence")+ 
    ylab("Multiple Dimensions Prevalence")+
    stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.7, label.y.npc=1),
  ncol = 2, nrow = 2)



dev.off()


