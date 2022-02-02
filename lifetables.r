#####################################################################################################
#################### Aim: download and prepare Lx and lx from life-tables
#####################################################################################################
rm(list=ls())
dev.off()

library(eurostat)
library(sjmisc)
library(tidyr)
library(dplyr)

data.dir <- "C:\\Users\\Magdalena\\demography\\happyindex\\data\\EUROSTAT"
setwd(data.dir)

ltdat <- get_eurostat("demo_mlifetable", time_format = "num") 

###the same country numbers as in SILC
csymb <- c("AT", "BE", "BG", "CH", "CY", "CZ", "DK", "EE","EL","ES","FI", "FR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "NL", "NO", "PL","PT",
            "RO", "RS", "SE", "SI", "SK","UK")
cno <- 1:29
ltdat$cno <- 0

for (i in 1:length(cno)){
  ltdat$cno[ltdat$geo==csymb[i]] <- cno[i]
}

Lx <- ltdat %>%
  dplyr::filter(indic_de == "PYLIVED", time==2018, geo %in% csymb, sex!="T")%>%
  mutate(age=replace(age,age=="Y_GE85", "Y85")) %>%
  mutate(age2=as.numeric(substr(age,2,3))) %>%
  mutate(agegr=0)%>%
  mutate(sex=replace(sex,sex=="F","W")) %>%  #the same order of sex as in SILC
  dplyr::filter(age2>=15) %>%
  select(geo,cno,sex,age2,agegr,values)
 
#replace age groups as in prevalence
ages <- c(15,30,40,seq(50,80,5))

for (i in 1:length(ages)){
  Lx$agegr[Lx$age2>=ages[i]] <- ages[i] 
}

l15 <- ltdat %>%
  dplyr::filter(indic_de == "SURVIVORS", time==2018, geo %in% csymb, age=="Y15",sex!="T") %>%
  mutate(sex=replace(sex,sex=="F","W")) %>%  #the same order of sex as in SILC
  select(geo,cno,sex,values) 
  

write.table(Lx, "Lx.csv", sep=",")
write.table(l15, "l15.csv", sep=",")

