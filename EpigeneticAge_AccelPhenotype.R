library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(kableExtra)
library(reshape)
library(corrplot)
library(Metrics)
library(factoextra)
library(NbClust)

dir_data <- "P:/UMICH/2022.01.12 - P01-MS-114/Data/"

epiage <- read.csv(paste0(dir_data, "P01_Epigenetic Age_202111.csv"))
epiage_new <- read.csv(paste0(dir_data, "epiage_newclocks.csv"))


sex <- read.csv(paste0(dir_data, "P01_Sexual_Maturation_201910.csv"))
sleep <- read.csv(paste0(dir_data, "P01_Sleep Calculations_Wide_202109.csv"))
pa <- read.csv(paste0(dir_data, "P01_Physical_Activity_Actigraph_202105.csv"))
gen <- read.csv(paste0(dir_data, "P01_General_Questionnaire_T1_202010.csv"))
pa_survey <- read.csv(paste0(dir_data, "P01_Physical_Activity_202105.csv"))
mat_edu <- read.csv(paste0(dir_data, "Historical_Mom_Educ_012022(T1 Only).csv"))
ses <- read.csv(paste0(dir_data, "P01_SES_201910.csv"))
age_decimal <- read.csv(paste0(dir_data,"P01_Sex_Age_202203.csv"))
anthro <- read.csv(paste0(dir_data, "P01_Anthropometry_withBMI_201911.csv"))

accel_wide <- read.csv(paste0(dir_data, "P01_Sleep Calculations_Wide_202109.csv"))


epiage_all <- merge(epiage, (epiage_new %>% filter(Sample_Group=="P01")), by="foliocc", all.x=TRUE)                  
age_vars <- epiage_all%>% select(foliocc, AgeAccelerationResidual, 
                                 DNAmAgeSkinBloodClockAdjAge,
                                 DNAmTLAdjAge,
                                 IEAA,
                                 EEAA,
                                 AgeAccelerationResidualHannum,
                                 AgeAccelGrim,
                                 AgeAccelPheno, 
                                 ageAcc2.PedBE,ageAcc2.Wu,
                                 PC1,PC2,PC3,PC4,PC5,
                                 CD8.naive,CD8pCD28nCD45RAn,PlasmaBlast,CD4T,NK,Mono,Gran) %>% 
  arrange(foliocc)

sex_vars=sex %>% filter(etapa=="T1") %>% select(etapa,foliocc, sexo, 
                                                p2_4a,p2_4b,pna5,v21) %>% 
  arrange(foliocc) %>% 
  mutate(left_test =cut(p2_4a, breaks=c(0,14.9,25),labels=c(0,1)),
         right_test=cut(p2_4b, breaks=c(0,14.9,25),labels=c(0,1)), 
         rl_test=ifelse((left_test==1 | right_test == 1),1,0))

sex_vars2=sex %>% filter(etapa=="T1") %>% select(etapa,foliocc, sexo, 
                                                 p2_2,p2_3,v21) %>% 
  mutate(p2_puberty=ifelse(sexo==1, ifelse(p2_3>3, 1, 0), ifelse(p2_2>3, 1,0 ))) %>% 
  arrange(foliocc)

sex_vars_all=sex %>% filter(etapa=="T1") %>% select(etapa,foliocc, sexo, 
                                                    p2_4a,p2_4b,pna5,v21,
                                                    p2_2,p2_3,v21) %>% 
  arrange(foliocc) %>% 
  mutate(left_test =cut(p2_4a, breaks=c(0,14.9,25),labels=c(0,1)),
         right_test=cut(p2_4b, breaks=c(0,14.9,25),labels=c(0,1)), 
         rl_test=ifelse((left_test==1 | right_test == 1),1,0)) %>% 
  mutate(p2_puberty=ifelse(sexo==1, ifelse(p2_3>3, 1, 0), ifelse(p2_2>3, 1,0 ))) %>% 
  arrange(foliocc)



pa_vars <- pa %>% filter(etapa=="T1") #%>% select(-c(TotalTime))
sleep_vars <- sleep %>% filter(etapa=="T1") %>% select(foliocc, sleep_dur_tot,midpoint_tot,fi_tot)
gen_vars <- gen %>% select(foliocc, edad, s2)%>% 
  mutate(age_cat=cut(edad, breaks=c(0,11.9, 13.9,15.9,18), labels=c("cat1","cat2","cat3","cat4")))

ses_vars <- ses %>% filter(etapa=="T1") %>% dplyr::select(foliocc, amai_13x6,amai_8x7) %>% 
  mutate(amai_13x6=as.factor(amai_13x6), amai_8x7=as.factor(amai_8x7))

anthro_vars <- anthro %>% filter(etapa=="T1") %>% 
  dplyr::mutate(foliocc=FOLIOCC) %>% 
  dplyr::select(foliocc, zbfa)

getScreenTime <- function(oldname){
  newvar=ifelse(oldname==1,0, ifelse(oldname==2,0.5,ifelse(oldname==3,1.5,
                                                           ifelse(oldname==4, 2.5, ifelse(oldname==5, 4.5, 
                                                                                          ifelse(oldname==6,6.5, ifelse(oldname==7,8,NA) ))))))
  return(newvar)
}

screen_vars <- pa_survey %>% filter(etapa=="T1") %>% 
  select(foliocc,acv4_1,acv4_2,acv4_3,acv4_4,acv4_5,acv4_6,acv4_7,
         acv5_1,acv5_2,acv7_1,acv7_2) %>% 
  mutate(scr1=getScreenTime(acv4_1),
         scr2=getScreenTime(acv4_2),
         scr3=getScreenTime(acv4_3),
         scr4=getScreenTime(acv4_4),
         scr5=getScreenTime(acv4_5),
         scr6=getScreenTime(acv4_6),
         scr7=getScreenTime(acv4_7),
         scr8=getScreenTime(acv5_1),
         scr9=getScreenTime(acv5_2),
         scr10=getScreenTime(acv7_1),
         scr11=getScreenTime(acv7_2),
         screentot= scr1 + scr2 +scr3 +scr4 +scr5 +  scr6+ scr7 + scr8*5 + scr9*5 + scr10*2 + scr11*2,
         screentot_cat=cut(screentot, breaks=c(0,22.5,32.5,48,116),labels=c("Q1","Q2","Q3","Q4")),
         screentot_TV=scr1 + scr2 +scr3 +scr4 +scr5 +  scr6+ scr7,
         screentot_movies= scr8*5+scr10*2,
         screentot_videogmes=scr9*5+scr11*2)

getCommuteTime <- function(oldname){
  newvar=ifelse(oldname==1,0.5, 
                ifelse(oldname==2,1.5,
                       ifelse(oldname==3,2.5,
                              ifelse(oldname==4, 3.5, 
                                     ifelse(oldname==5, 4,NA)))))
}

commute_vars=pa_survey %>% filter(etapa=="T1") %>% 
  select(foliocc,acv9_1,acv9_2) %>% 
  mutate(commute_wkday=getCommuteTime(acv9_1), 
         commute_wkend=getCommuteTime(acv9_2), 
         commute_tot=commute_wkday*5 + commute_wkend*2)

matedu_var <- mat_edu %>% 
  mutate(foliocc=FOLIOCC) %>% 
  mutate(school_cat=cut(SCHOOL_T1, breaks=c(0,8,11,12,max(SCHOOL_T1, na.rm=TRUE)))) %>% 
  select(foliocc,school_cat) 

age_dec_vars <- age_decimal %>% 
  dplyr::mutate(foliocc=FOLIOCC) %>% 
  dplyr::select(foliocc, Age_P01)



sleep_variance=accel_wide %>% filter(etapa=="T1") %>% select(foliocc,sleep_dur_tot_sd)
screen_vars_tot <- screen_vars %>% 
  dplyr::select(foliocc, screentot)

###Merge data sets

sex_test <- merge(merge(merge(merge(merge(merge(merge(age_vars, sex_vars_all, by="foliocc",all=T),
                                                screen_vars[c("foliocc","screentot_cat")], by="foliocc"),
                                          matedu_var, by="foliocc", all=T), 
                                    gen_vars[c(1:2,4)], by="foliocc",all=T),
                              ses_vars, by="foliocc", all=T),
                        age_dec_vars, by="foliocc", all=T), 
                  anthro_vars, by="foliocc", all=T)


epi_accel <- merge(merge(merge(merge(sex_test, sleep_vars, by="foliocc", all.x=T),
                               pa_vars, by="foliocc",all.x=T), 
                         sex_vars2, by="foliocc", all.x=T), 
                   anthro_vars, by="foliocc", all.x=T)



## Clustering
accel_vars=merge(merge(merge(pa_vars, sleep_vars, by="foliocc"),
                       sleep_variance, by="foliocc" ),
                 screen_vars_tot, by="foliocc")%>% 
  dplyr::select(-c(TotalTime))


### Select Number of Clusters
NbClust(data = scale(accel_vars[c(3:11)]),  distance = "euclidean",
        min.nc = 2, max.nc = 15, method = "kmeans") ## best =3 gps

NbClust(data = scale(accel_vars_boys[c(3:11)]),  distance = "euclidean",
        min.nc = 2, max.nc = 15, method = "kmeans")  ## best = 2 groups, 2nd best = 3gps

NbClust(data = scale(accel_vars_girls[c(3:11)]),  distance = "euclidean",
        min.nc = 2, max.nc = 15, method = "kmeans") ##best = 3 groups

## each chooses 2 or 3, majority chooses 3


## Create Clusters

set.seed(124);k_final=kmeans(accel_vars[-c(1:2)], centers = 3,
                             nstart = 25)
set.seed(124);k_2= kmeans(accel_vars[-c(1:2)], centers = 2,
                          nstart = 25)
accel_vars$clus2=k_2$cluster
accel_vars$cluster=k_final$cluster



#####################################
######### CLuster as Covariate #####
#####################################

epi_accel_clus <- merge(sex_test, accel_vars, by="foliocc", all.y=T)


x_age=c("AgeAccelerationResidual", "DNAmAgeSkinBloodClockAdjAge",
        "IEAA", "EEAA", "AgeAccelerationResidualHannum", "AgeAccelGrim", 
        "AgeAccelPheno", "ageAcc2.PedBE", "ageAcc2.Wu")


subs <- reshape2::melt(epi_accel_clus[c("foliocc",x_age,"cluster","clus2","sexo","edad","school_cat", "CD8.naive","CD4T","Gran", "amai_13x6", "amai_8x7", "v21","zbfa", "p2_puberty")], 
                       id.var=c("foliocc","cluster","clus2","sexo","edad","school_cat", "CD8.naive","CD4T","Gran","amai_13x6", "amai_8x7", "v21","zbfa","p2_puberty"), value.name="Epigenetic_Age") 


results1 <- subs%>% group_by(variable) %>% filter(!variable %in% c("IEAA","EEAA")) %>% 
  do(fitEpiAge = (tidy(lm(Epigenetic_Age ~ relevel(as.factor(cluster),3) +as.factor(sexo) +  CD8.naive + CD4T + Gran + school_cat +amai_8x7+p2_puberty+ zbfa , data = .),conf.int=TRUE))) %>%
  unnest(fitEpiAge) %>% filter(term %in% c("relevel(as.factor(cluster), 3)1","relevel(as.factor(cluster), 3)2"))  %>% 
  mutate(EpiAge_Outcome=variable) %>% 
  select(EpiAge_Outcome, term, estimate,p.value, conf.low, conf.high)

results_nocell <- subs%>% filter(variable %in% c("IEAA","EEAA")) %>% group_by(variable) %>%
  do(fitEpiAge = (tidy(lm(Epigenetic_Age ~ relevel(as.factor(cluster),3) +as.factor(sexo)+ school_cat + amai_8x7+p2_puberty+ zbfa , data = .),conf.int=TRUE))) %>%
  unnest(fitEpiAge) %>% filter(term %in% c("relevel(as.factor(cluster), 3)1","relevel(as.factor(cluster), 3)2"))  %>% 
  mutate(EpiAge_Outcome=variable) %>% 
  select(EpiAge_Outcome, term, estimate,p.value, conf.low, conf.high)

############################
## Sex-Specific Analysis ##
############################
accel_vars_sexstrat=merge(merge(merge(merge(pa_vars, sleep_vars, by="foliocc"),
                                      sleep_variance, by="foliocc" ),
                                screen_vars_tot, by="foliocc"),
                          sex_vars %>% dplyr::select(foliocc,sexo),by="foliocc")%>% 
  dplyr::select(-c(TotalTime))

accel_vars_girls <- accel_vars_sexstrat %>% filter(sexo==2) %>% select(-sexo)
accel_vars_boys <- accel_vars_sexstrat %>% filter(sexo==1) %>% select(-sexo)

epi_accel_clus_b <- merge(sex_test, accel_vars_boys, by="foliocc", all.y=T)
epi_accel_clus_g <- merge(sex_test, accel_vars_girls, by="foliocc", all.y=T)


set.seed(1234);k_final_b=kmeans(accel_vars_boys[c(3:11)], centers = 3,
                                nstart = 25)

set.seed(1234);k_final_g=kmeans(accel_vars_girls[c(3:11)], centers = 3,
                                nstart = 25)


epi_accel_clus_b$cluster3=k_final_b$cluster
epi_accel_clus_g $cluster3=k_final_g$cluster


## Cluster as Covariate

subs_b <- reshape2::melt(epi_accel_clus_b[c("foliocc",x_age,"cluster3","sexo","edad","school_cat", "CD8.naive","CD4T","Gran", "amai_13x6", "amai_8x7", "v21",
                                            "zbfa","p2_puberty")], 
                         id.var=c("foliocc","cluster3","sexo","edad","school_cat", "CD8.naive","CD4T","Gran","amai_13x6", "amai_8x7", "v21",
                                  "zbfa","p2_puberty"), value.name="Epigenetic_Age") 

subs_g <- reshape2::melt(epi_accel_clus_g[c("foliocc",x_age,"cluster3","sexo","edad","school_cat", "CD8.naive","CD4T","Gran", "amai_13x6", "amai_8x7", "v21",
                                            "zbfa","p2_puberty")], 
                         id.var=c("foliocc","cluster3","sexo","edad","school_cat", "CD8.naive","CD4T","Gran","amai_13x6", "amai_8x7", "v21",
                                  "zbfa","p2_puberty"), value.name="Epigenetic_Age") 



epi_accel_clus_b3  <- subs_b%>% filter(! variable %in% c("IEAA","EEAA"))%>% group_by(variable) %>%
  do(fitEpiAge = (tidy(lm(Epigenetic_Age ~ relevel(as.factor(cluster3),3)  + CD8.naive + CD4T + Gran + school_cat+ amai_8x7 + zbfa+ p2_puberty, data = .),conf.int=TRUE))) %>%
  unnest(fitEpiAge) %>% filter(term %in% c("relevel(as.factor(cluster3), 3)1","relevel(as.factor(cluster3), 3)2"))  %>% 
  mutate(EpiAge_Outcome=variable) %>% 
  select(EpiAge_Outcome, term, estimate,p.value, conf.low, conf.high)

epi_accel_clus_b3_nocell  <- subs_b%>% filter(variable %in% c("IEAA","EEAA"))%>% group_by(variable) %>%
  do(fitEpiAge = (tidy(lm(Epigenetic_Age ~ relevel(as.factor(cluster3),3)  +school_cat+ amai_8x7 + zbfa+ p2_puberty, data = .),conf.int=TRUE))) %>%
  unnest(fitEpiAge) %>% filter(term %in% c("relevel(as.factor(cluster3), 3)1","relevel(as.factor(cluster3), 3)2"))  %>% 
  mutate(EpiAge_Outcome=variable) %>% 
  select(EpiAge_Outcome, term, estimate,p.value, conf.low, conf.high)


epi_accel_clus_g3  <- subs_g%>% filter(! variable %in% c("IEAA","EEAA"))%>% group_by(variable) %>%
  do(fitEpiAge = (tidy(lm(Epigenetic_Age ~ relevel(as.factor(cluster3),3)  + CD8.naive + CD4T + Gran + school_cat+ amai_8x7 + zbfa+ p2_puberty, data = .),conf.int=TRUE))) %>%
  unnest(fitEpiAge) %>% filter(term %in% c("relevel(as.factor(cluster3), 3)1","relevel(as.factor(cluster3), 3)2"))  %>% 
  mutate(EpiAge_Outcome=variable) %>% 
  select(EpiAge_Outcome, term, estimate,p.value, conf.low, conf.high)

epi_accel_clus_g3_nocell  <- subs_g%>% filter(variable %in% c("IEAA","EEAA"))%>% group_by(variable) %>%
  do(fitEpiAge = (tidy(lm(Epigenetic_Age ~ relevel(as.factor(cluster3),3)  +school_cat+ amai_8x7 + zbfa+ p2_puberty, data = .),conf.int=TRUE))) %>%
  unnest(fitEpiAge) %>% filter(term %in% c("relevel(as.factor(cluster3), 3)1","relevel(as.factor(cluster3), 3)2"))  %>% 
  mutate(EpiAge_Outcome=variable) %>% 
  select(EpiAge_Outcome, term, estimate,p.value, conf.low, conf.high)
