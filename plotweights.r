#####################################################################################################
#################### Aim: make a scatterplot of weights
#####################################################################################################
rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Rcpp)
library(ggpubr)


data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\SILC"
setwd(data.dir)


########################################################################################################
#####average weights by country
########################################################################################################
chronicallw <- read.table(file="chronicweights.csv", sep=",", header=TRUE)
selfweights <- read.table(file="selfweights.csv", sep=",", header=TRUE)
GALIallw <- read.table(file="GALIweights.csv", sep=",", header=TRUE)
multiweights <- read.table(file="multiweights.csv", sep=",", header=TRUE)

#average by country
bycountrysr <-  selfweights %>%
  pivot_longer(cols =  SR2:SR5,
               names_to = "sr", values_to = "weight") %>%
  group_by(country,sr) %>%
  summarise(avg = mean(weight)) %>%
  pivot_wider(names_from = sr,
              values_from = avg)

bycountrychr <-  chronicallw %>%
  group_by(country) %>%
  summarise(chr = mean(chronic)) 

bycountryGALI <-  GALIallw %>%
  pivot_longer(cols = GALI1:GALI2,
               names_to = "GALI", values_to = "weight") %>%
  group_by(country,GALI) %>%
  summarise(avg = mean(weight)) %>%
  pivot_wider(names_from = GALI,
              values_from = avg)

bycountryall <- cbind(bycountrychr,bycountryGALI[,-1],bycountrysr[,-1])
write.table(bycountryall, file="allweightsaver.csv", sep=",")

bycountrymul <-  multiweights %>%
  group_by(country,level) %>%
  summarise(avg = mean(weight)) %>%
  pivot_wider(names_from = level,
              values_from = avg)
write.table(bycountrymul, file="bycountrymul.csv", sep=",")



##############################################################################################
######## plots for average by country
##############################################################################################
allweights <- read.table(file="allweightsaver.csv", sep=",", header=TRUE)


multiw <- read.table(file="bycountrymul.csv", sep=",", header=TRUE)
multiw[multiw <0] <- NA
sellevels <- c("country","X12","X123")
multiw2 <- multiw %>%
  select(sellevels)


allweights <- cbind(allweights,multiw2[,-1])

bycountry <-  allweights %>%
  pivot_longer(cols = c(3:ncol(allweights)),
               names_to = "health", values_to = "weight") 

p1 <- ggplot(subset(bycountry, (health=="GALI1"| health=="GALI2")) , aes(chr, weight,col=health, label=country))  
p2 <- ggplot(subset(bycountry, (health=="SR2" | health=="SR3")) , aes(chr, weight,col=health, label=country))  
p3 <- ggplot(subset(bycountry, (health=="SR4" | health=="SR5")) , aes(chr, weight,col=health, label=country))  
p4 <- ggplot(subset(bycountry, (health %in% sellevels)) , aes(chr, weight,col=health, label=country))  



dev.off()
pdf(file = "C:/Users/Magdalena/demography/happyindex/figures/weightsc.pdf", width = 10, height = 9)

ggarrange(
  p2 + geom_point()+ # cor1=0.10, cor2=0.37**
      theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
  #   legend.title = element_blank(),
      legend.key=element_blank(), legend.position = "top",  legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
      scale_color_manual("Self-rated Health",values=c("black", "grey70"), labels=c("Very Good","Good"))+
      xlab("Chronic Diseases Weights")+ xlim(0.78,0.93)+
      scale_y_continuous(name="Self-rated Health Weights", breaks=seq(0.35,0.95,0.2), labels=waiver(), limits=c(0.33,0.95))+
      stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.5, label.y.npc=0.2),
  
  p3 + geom_point()+ #cor1=0.47**, cor2=0.44**
      theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(),legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
      scale_color_manual("Self-rated Health", values=c("black", "grey70"), labels=c("Fair", "Poor"))+
      xlab("Chronic Diseases Weights")+ xlim(0.78,0.93)+
      scale_y_continuous(name="Self-rated Health Weights", breaks=seq(0.35,0.95,0.2), labels=waiver(), limits=c(0.33,0.95))+
      stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.5, label.y.npc=0.95),
  
  p1 + geom_point()+ #cor1=0.57***, cor2=0.68***
      theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(), legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
    scale_color_manual("GALI", values=c("black", "grey70"), labels=c("Moderately Limited", "Severely Limited"))+
    xlab("Chronic Diseases Weights")+ xlim(0.78,0.93)+
    scale_y_continuous(name="GALI Weights", breaks=seq(0.35,0.95,0.2), labels=waiver(), limits=c(0.33,0.95))+
    stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.5, label.y.npc=0.2),
  
  p4 + geom_point()+ # cor=0.15, cor=0.48**
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"), 
      legend.key=element_blank(), legend.position = "top", legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))+
    scale_color_manual("Dimensions' Levels:", values=c("black", "grey70"), labels=c("No,Not,Very Good", "Yes,Moderately,Good"))+
    xlab("Chronic Diseases Weights")+ xlim(0.78,0.93)+
    scale_y_continuous(name="Multiple Dimensions Weights", breaks=seq(0.35,0.95,0.2), labels=waiver(), limits=c(0.33,0.95))+
    stat_cor(method="pearson", cor.coef.name="r", r.accuracy=0.01, label.x.npc=0.5, label.y.npc=0.2),
  ncol = 2, nrow = 2)



dev.off()


########################################################################################################
#average by sex - not in tables/graphs
bysexsr <-  selfweights %>%
  pivot_longer(cols = SR2:SR5,
               names_to = "sr", values_to = "weight") %>%
  group_by(sex,sr) %>%
  summarise(avg = mean(weight)) %>%
  pivot_wider(names_from = sr,
              values_from = avg)

bysexchr <-  chronicallw %>%
  group_by(sex) %>%
  summarise(avg = mean(chronic)) 


bysexGALI <-  GALIallw %>%
  pivot_longer(cols = GALI1:GALI2,
               names_to = "GALI", values_to = "weight") %>%
  group_by(sex,GALI) %>%
  summarise(avg = mean(weight)) %>%
  pivot_wider(names_from = GALI,
              values_from = avg)


bysexmul <-  multiweights %>%
  group_by(sex,level) %>%
  summarise(avg = mean(weight)) %>%
  pivot_wider(names_from = level,
              values_from = avg)
