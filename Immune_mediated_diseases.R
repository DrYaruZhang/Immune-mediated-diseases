####data process####
dementia<-ACD_AD_VD_covariants[,c(2,4,5,9,11,12,16,18,19,23,24:36,39)]
dementia_AID_0<-merge(FinalData[,c(1:56,60:65)],dementia,by="eid",all = TRUE)
dementia_AID_1<-subset(dementia_AID_0,dementia_AID_0$eid %in% HESIN_eid)
dementia_AID<-dementia_AID_1[,c(1:2,63:71,3:62,72:85)]

str(dementia_AID)
date_position<-grep("^Date",colnames(dementia_AID))
date_position
for (i in 1:length(date_position)) {
  dementia_AID[,date_position[i]]<-as.Date(dementia_AID[,date_position[i]], origin = "1970-01-01")
}
dementia_AID$X53.0.0<-as.Date(dementia_AID$X53.0.0,"%Y-%m-%d")
dementia_AID[, cols <- grep("_date$", names(dementia_AID))] <- 
  lapply(dementia_AID[, cols <- grep("_date$", names(dementia_AID))], as.Date, format = "%Y/%m/%d")
str(dementia_AID)

a<-subset(dementia_AID,dementia_AID$dementia_date=="2037-07-07" | dementia_AID$AD_date=="2037-07-07"
          | dementia_AID$VD_date=="2037-07-07")
table(a$dementia_status)
table(a$AD_status)
table(a$VD_status)
dementia_AID$dementia_date[dementia_AID$dementia_date=="2037-07-07"]<-"2021-04-07"
dementia_AID$dementia_years<-difftime(dementia_AID$dementia_date,dementia_AID$X53.0.0,units ="days")/365
dementia_AID$dementia_years<-as.numeric(dementia_AID$dementia_years)
dementia_AID$AD_years<-difftime(dementia_AID$AD_date,dementia_AID$X53.0.0,units ="days")/365
dementia_AID$AD_years<-as.numeric(dementia_AID$AD_years)
dementia_AID$VD_years<-difftime(dementia_AID$VD_date,dementia_AID$X53.0.0,units ="days")/365
dementia_AID$VD_years<-as.numeric(dementia_AID$VD_years)

#exclusion of dementia patients at baseline
ACD_losefl<-subset(dementia_AID,is.na(dementia_AID$dementia_years)|(dementia_AID$dementia_years<=0 & dementia_AID$dementia_status==0))
ACD_bl<-subset(dementia_AID,dementia_AID$dementia_years<=0 & dementia_AID$dementia_status==1)
dementia_AID$dementia_status[dementia_AID$dementia_years<=0]<-NA
x1<-subset(dementia_AID,!is.na(dementia_AID$dementia_status))

AD_losefl<-subset(x1,is.na(x1$AD_years)|(x1$AD_years<=0 & x1$AD_status==0))
AD_bl<-subset(x1,x1$AD_years<=0 & x1$AD_status==1)
x1$AD_status[x1$AD_years<=0]<-NA

VD_losefl<-subset(x1,is.na(x1$VD_years)|(x1$VD_years<=0 & x1$VD_status==0))
VD_bl<-subset(x1,x1$VD_years<=0 & x1$VD_status==1)
x1$VD_status[x1$VD_years<=0]<-NA

dementia_bl<-merge(ACD_bl,AD_bl,by="eid",all=TRUE)
dementia_losefl<-ACD_losefl

x2<-subset(x1,!is.na(x1$AD_status) & !is.na(x1$VD_status))

x2$AD_status[x2$AD_status==0 & x2$dementia_status==1]<-NA
x2$VD_status[x2$VD_status==0 & x2$dementia_status==1]<-NA
table(x2$dementia_status)
table(x2$AD_status)
table(x2$VD_status)
#write.csv(x2,"x2.csv")

for (i in c(1:375894))
{x2$min[i] <- min(t(x2[i,c(13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                         52, 55, 58, 61, 64, 67, 70)]),na.rm=TRUE)}
range(x2$min,na.rm = TRUE)

#write.csv(x2,"/share/inspurStorage/home1/Royce/zyr_data/Immune_Mediated_Diseases/AID_Date.csv")

###Any_AID before dementia##
x2<-AID_Date[,2:87]

str(x2)
x2$X53.0.0<-as.Date(x2$X53.0.0,"%Y-%m-%d")
x2[, cols <- grep("_date$", names(x2))] <- 
  lapply(x2[, cols <- grep("_date$", names(x2))], as.Date, format = "%Y-%m-%d")
x2[, cols <- grep("^Date_", names(x2))] <- 
  lapply(x2[, cols <- grep("^Date_", names(x2))], as.Date, format = "%Y-%m-%d")
str(x2)

x2$first_AID<-as.numeric((difftime(x2$dementia_date,x2$min,units = c("days")))/365)
range(x2$first_AID,na.rm = TRUE)
x2[,12:71][x2[,12:71]==""]<-NA
x2$sum_all_AID<-rowSums(x2[,c(14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 
                              53, 56, 59, 62, 65, 68, 71)],na.rm = TRUE)
table(x2$sum_all_AID)
sum(x2$sum_all_AID>0,na.rm = TRUE)
x2$Any_AID<-ifelse(x2$sum_all_AID>0 & x2$first_AID>1,1,0)#1years
table(x2$Any_AID)

x2$Any_AID2<-x2$Any_AID#1-2years
x2$Any_AID2[x2$first_AID>2]<-NA
table(x2$Any_AID2)
x2$Any_AID3<-x2$Any_AID#1-3years
x2$Any_AID3[x2$first_AID>3]<-NA
table(x2$Any_AID3)
x2$Any_AID4<-x2$Any_AID#1-4years
x2$Any_AID4[x2$first_AID>4]<-NA
table(x2$Any_AID4)
x2$Any_AID5<-x2$Any_AID#1-5years
x2$Any_AID5[x2$first_AID>5]<-NA
table(x2$Any_AID5)
x2$Any_AID10<-x2$Any_AID#1-10years
x2$Any_AID10[x2$first_AID>10]<-NA
table(x2$Any_AID10)
x2$Any_AID15<-x2$Any_AID#1-15years
x2$Any_AID15[x2$first_AID>15]<-NA
table(x2$Any_AID15)
x2$Any_AID20<-x2$Any_AID#1-20years
x2$Any_AID20[x2$first_AID>20]<-NA
table(x2$Any_AID20)


date_position<-grep("^Date",colnames(x2))
date_position
diag_position<-date_position+1
diag_position

x2[,c(14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65, 68, 
      71)][is.na(x2[,c(14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 
                       59, 62, 65, 68, 71)])]<-0
for (i in 1:length(date_position)) {
  diag_name<-substring(names(x2)[date_position[i]],6:nchar(names(x2)[date_position[i]]))[1]
  diff_name<-paste0("dementia","_",diag_name)
  x2[diff_name]<-as.numeric(difftime(x2$dementia_date,x2[,date_position[i]],units = c("days"))/365)
  x2[diag_name][is.na(x2[diag_name])==TRUE|x2[diff_name]<1]<-0
  x2[diag_name][x2[diag_name]==0 & x2$Any_AID==1]<-NA
}

x2$sum_individual_AID<-rowSums(x2[,c(14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 
                                     53, 56, 59, 62, 65, 68, 71)],na.rm = TRUE)
table(x2$sum_individual_AID)
sum(x2$sum_individual_AID>0,na.rm = TRUE)
x3<-x2[,c(1:11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 
          53, 56, 59, 62, 65, 68, 71, 72:96,117)]

#write.csv(x3,"Dementia_AID_all.csv")

####demographic####
a<-Dementia_AID

library(table1)
library(lubridate)

table(a$APOE4)
a$APOE4[a$APOE4==2]<-1
table(a$APOE4)

table(a$Qualification)
a$Qualification[a$Qualification==6]<-1
a$Qualification[a$Qualification<6 & a$Qualification>1]<-0
table(a$Qualification)

table(a$race_group)
a$race_group[a$race_group>1]<-0
table(a$race_group)

a$dementia_status<-factor(a$dementia_status,levels = c(0,1),labels = c("No Incident Dementia","Incident ACD"))
a$AD_status<-factor(a$AD_status,levels = c(0,1),labels = c("No Incident Dementia","Incident AD"))
a$VD_status<-factor(a$VD_status,levels = c(0,1),labels = c("No Incident Dementia","Incident VD"))
a$Sex<-factor(a$Sex,levels = c(0,1),labels = c("Female, n (%)","Males"))
a$APOE4<-factor(a$APOE4,levels = c(0,1),labels = c("none APOE4","ApoE ε4 carriers, n (%)"))
a$Qualification<-factor(a$Qualification,levels = c(0,1),labels = c("Lower education","Higher education, n (%)"))
a$race_group<-factor(a$race_group,levels = c(0,1),labels = c("non_White","Ethnicity_White, n (%)"))
a$Smoking<-factor(a$Smoking,levels = c(0,1,2),labels = c("Never, n (%)","Previous, n (%)","Current, n (%)"))
a$Alcohol<-factor(a$Alcohol,levels = c(0,1,2),labels = c("Never, n (%)","Previous, n (%)","Current, n (%)"))
a$Any_AID<-factor(a$Any_AID,levels = c(0,1),labels = c("None, n (%)","Any, n (%)"))

label(a$Age)<- "Age, y"
label(a$Sex)<- "Female, n (%)"
label(a$APOE4)<- "ApoE ε4 carriers, n (%)"
label(a$Qualification)<- "Higher education, n (%)"
label(a$race_group)<- "Ethnicity_White, n (%)"
label(a$BMI)<- "BMI, Kg/m2"
label(a$Townsend)<- "Townsend deprivation score"
label(a$Smoking)<- "Smoking"
label(a$Alcohol)<- "Alcohol consumption"
label(a$Any_AID)<- "Immune-mediated diseases"

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=4), c("",
                                                           "Mean (SD)"=sprintf("%s (%s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.2f%%)", FREQ, PCT))))
}

table1(~ Age + Sex + APOE4 + Qualification + race_group + BMI
       + Townsend + Smoking + Alcohol + Any_AID| dementia_status, data=a, overall=F,
       render.continuous=my.render.cont, render.categorical=my.render.cat)

table1(~ Age + Sex + APOE4 + Qualification + race_group + BMI
       + Townsend + Smoking + Alcohol+ Any_AID| AD_status, data=a, overall=F,
       render.continuous=my.render.cont, render.categorical=my.render.cat)

table1(~ Age + Sex + APOE4 + Qualification + race_group + BMI
       + Townsend + Smoking + Alcohol+ Any_AID| VD_status, data=a, overall=F,
       render.continuous=my.render.cont, render.categorical=my.render.cat)

####COX####
a<-Dementia_AID

####______cox_any####
name_a<-unique(colnames(Dementia_AID[,c(13:32,51:57,50)]))
name_a
library(dplyr)
library(survival)
library(ggfortify)
library(survminer)
library(forcats)

res.cox<-coxph(Surv(dementia_years, dementia_status==1)~Any_AID,data=a)
summary(res.cox)
res.cox<-coxph(Surv(dementia_years, dementia_status==1)~Any_AID+Age+Sex+APOE4+Qualification,data=a)
summary(res.cox)
res.cox<-coxph(Surv(dementia_years, dementia_status==1)~Any_AID+Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,data=a)
summary(res.cox)

res.cox<-coxph(Surv(AD_years, AD_status==1)~Any_AID,data=a)
summary(res.cox)
res.cox<-coxph(Surv(AD_years, AD_status==1)~Any_AID+Age+Sex+APOE4+Qualification,data=a)
summary(res.cox)
res.cox<-coxph(Surv(AD_years, AD_status==1)~Any_AID+Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,data=a)
summary(res.cox)

res.cox<-coxph(Surv(VD_years, VD_status==1)~Any_AID,data=a)
summary(res.cox)
res.cox<-coxph(Surv(VD_years, VD_status==1)~Any_AID+Age+Sex+APOE4+Qualification,data=a)
summary(res.cox)
res.cox<-coxph(Surv(VD_years, VD_status==1)~Any_AID+Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,data=a)
summary(res.cox)

####survival curve###
library(ggfortify)
fit<- survfit(Surv(dementia_years,dementia_status==1)~Any_AID,data = a)
ggsurvplot(fit,conf.int = FALSE,size = 1,palette = c("#4169E1", "#FF6600"),
           risk.table =TRUE,risk.table.fontsize=3,tables.y.text=F,
           xlim=c(0,12),break.x.by=2,ylim=c(0.95,1),break.y.by=0.01,censor = F,
           xlab ="Follow up (Years)", ylab='Non-progression Proportion',legend.labs=c('None','Any'), 
           legend.title="Immune-mediated diseases",ncensor.plot = FALSE,test.for.trend= TRUE,
           ggtheme = theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
             theme(panel.background=element_blank()))

fit<- survfit(Surv(AD_years,AD_status==1)~Any_AID,data = a)
ggsurvplot(fit,conf.int = FALSE,size = 1,palette = c("#4169E1", "#FF6600"),
           risk.table =TRUE,risk.table.fontsize=3,tables.y.text=F,
           xlim=c(0,12),break.x.by=2,ylim=c(0.95,1),break.y.by=0.01,censor = F,
           xlab ="Follow up (Years)", ylab='Non-progression Proportion',legend.labs=c('None','Any'), 
           legend.title="Immune-mediated diseases",ncensor.plot = FALSE,test.for.trend= TRUE,
           ggtheme = theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
             theme(panel.background=element_blank()))

fit<- survfit(Surv(VD_years,VD_status==1)~Any_AID,data = a)
ggsurvplot(fit,conf.int = FALSE,size = 1,palette = c("#4169E1", "#FF6600"),
           risk.table =TRUE,risk.table.fontsize=3,tables.y.text=F,
           xlim=c(0,12),break.x.by=2,ylim=c(0.95,1),break.y.by=0.01,censor = F,
           xlab ="Follow up (Years)", ylab='Non-progression Proportion',legend.labs=c('None','Any'), 
           legend.title="Immune-mediated diseases",ncensor.plot = FALSE,test.for.trend= TRUE,
           ggtheme = theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
             theme(panel.background=element_blank()))

####______cox_individual####
#unadjusted
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a,!is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i)))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a,!is.na(AD_status),!is.na(AD_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i)))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a,!is.na(VD_status),!is.na(VD_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i)))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2


#adjusted for Age+Sex+ApoE4+education+race+smoking+alcohol+BMI+Townsend
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

####Mediation####
data1<-AID_dementia_blood[,2:65]
str(data1)

data1_all<-subset(data1,data1$T_dementia_blood>0 & data1$T_blood_AID>0)
summary(data1_all$first_AID)
table(data1_all$Any_AID)
data<-data1_all


library(mediation)
data$dementia_status<-factor(data$dementia_status,levels = c(0,1))
data$AD_status<-factor(data$AD_status,levels = c(0,1))
data$VD_status<-factor(data$VD_status,levels = c(0,1))

####ACD_AnyAID_Neutrophill##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ Any_AID 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Neutrophill) ~ Any_AID 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ scale(Neutrophill) 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Neutrophill) + Any_AID
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='Any_AID', mediator='scale(Neutrophill)')
summary(results2)
####ACD_AnyAID_Lymphocyte##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ Any_AID 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Lymphocyte) ~ Any_AID 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ scale(Lymphocyte) 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Lymphocyte) + Any_AID
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='Any_AID', mediator='scale(Lymphocyte)')
summary(results2)
####ACD_DiabetesMellitus..type.one._Neutrophill##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ DiabetesMellitus..type.one. 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Neutrophill) ~ DiabetesMellitus..type.one. 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Neutrophill) + DiabetesMellitus..type.one.
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='DiabetesMellitus..type.one.', mediator='scale(Neutrophill)')
summary(results2)
####ACD_DiabetesMellitus..type.one._Lymphocyte##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ DiabetesMellitus..type.one. 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Lymphocyte) ~ DiabetesMellitus..type.one. 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Lymphocyte) + DiabetesMellitus..type.one.
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='DiabetesMellitus..type.one.', mediator='scale(Lymphocyte)')
summary(results2)
####ACD_Rheumatic.fever...rheumatic.heart.diseases_Neutrophill##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ Rheumatic.fever...rheumatic.heart.diseases 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Neutrophill) ~ Rheumatic.fever...rheumatic.heart.diseases 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Neutrophill) + Rheumatic.fever...rheumatic.heart.diseases
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='Rheumatic.fever...rheumatic.heart.diseases', mediator='scale(Neutrophill)')
summary(results2)
####ACD_Rheumatic.fever...rheumatic.heart.diseases_Lymphocyte##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ Rheumatic.fever...rheumatic.heart.diseases 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Lymphocyte) ~ Rheumatic.fever...rheumatic.heart.diseases 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Lymphocyte) + Rheumatic.fever...rheumatic.heart.diseases
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='Rheumatic.fever...rheumatic.heart.diseases', mediator='scale(Lymphocyte)')
summary(results2)
####ACD_Multiple.sclerosis_Neutrophill##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ Multiple.sclerosis 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Neutrophill) ~ Multiple.sclerosis 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Neutrophill) + Multiple.sclerosis
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='Multiple.sclerosis', mediator='scale(Neutrophill)')
summary(results2)
####ACD_Multiple.sclerosis_Lymphocyte##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ Multiple.sclerosis 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Lymphocyte) ~ Multiple.sclerosis 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Lymphocyte) + Multiple.sclerosis
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='Multiple.sclerosis', mediator='scale(Lymphocyte)')
summary(results2)
####ACD_Necrotizing.vasculopathies_Neutrophill##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ Necrotizing.vasculopathies 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Neutrophill) ~ Necrotizing.vasculopathies 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Neutrophill) + Necrotizing.vasculopathies
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='Necrotizing.vasculopathies', mediator='scale(Neutrophill)')
summary(results2)
####ACD_Necrotizing.vasculopathies_Lymphocyte##
fit.totaleffect<-glm(relevel(dementia_status,ref="0") ~ Necrotizing.vasculopathies 
                     +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                     family = binomial(link = "logit"),data = data)
summary(fit.totaleffect)
fit.mediator<-lm(scale(Lymphocyte) ~ Necrotizing.vasculopathies 
                 +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
                 data = data)
summary(fit.mediator)
fit.dv=glm(relevel(dementia_status,ref="0") ~ scale(Lymphocyte) + Necrotizing.vasculopathies
           +Age+Sex+APOE4+Qualification+race_group+Smoking+Alcohol+BMI+Townsend,
           family = binomial(link = "logit"),data = data)
summary(fit.dv)
results2 = mediate(fit.mediator, fit.dv, sims=1000, treat='Necrotizing.vasculopathies', mediator='scale(Lymphocyte)')
summary(results2)
