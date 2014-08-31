#############################################################
#script to calculate scaling coeffients between different####
#measures of logging intensity and damage####################
#############################################################

#first load packages
library(nlme)
library(ggplot2)
library(MuMIn)
library(lme4)
library(nlme)

#load data
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Prop_dam<-read.csv("prop_damage.csv")
head(Prop_dam)
colnames(Prop_dam)<-c("Study","Site_ID","Age","Method","BA_log","Prop_BA_log","Vol_log","Tree_ex_ha","Dam_tree","Dam_ha","Sev_per_tree","Severe_ha","BA_dam","Prop_ba_dam","Prop_seve_dam","Prop_dam","ID","All","Region","N_logged","Plot","Notes")
Rich<-read.csv("Logged_Richness.csv")
AGB<-read.csv("Logged_AGB_rev.csv")

###############################################################
#organise all data so that there is data on number of trees####
#logged per hectare as well as the volume logged and region####
###############################################################

#first for stem damage
keeps1<-c("Study","Site_ID","Vol_log","Tree_ex_ha","Region")
Prop_dam2<-Prop_dam[keeps1]
colnames(Prop_dam2)<-c("Study","Site","Vol","Trees","Region")
#now for richness
keeps2<-c("Study.ID","Site.ID","Volume.extracted","Trees_extracted","Region")
Rich2<-Rich[keeps2]
colnames(Rich2)<-c("Study","Site","Vol","Trees","Region")
#now for AGB
str(AGB)
keeps3<-c("Study","Site","Vol","Trees_logged","Region")
AGB2<-AGB[keeps3]
colnames(AGB2)<-c("Study","Site","Vol","Trees","Region")
Int<-rbind(Prop_dam2,Rich2,AGB2)

#loop to change field values so they all conform
Int[Int==levels(Int$Region)[4]]<-"Asia"
Int$Region<-factor(Int$Region)
Int[Int==-9999]<-NA
#create variable to group both Africa and Asia together as these
#regions tend to have larger tree sizes and thus greater volume
#per tree
levels(Int$Region)[which(levels(Int$Region)=="Asia")] <- "Asia/Africa"
levels(Int$Region)[which(levels(Int$Region)=="Africa")] <- "Asia/Africa"

ggplot(Int,aes(x=Trees,y=Vol,colour=Region,group=Region))+geom_point(size=3)+geom_smooth(se=F,method="lm")


################################################################
#model of volume as a function of number of trees logged per ha#
################################################################

#create dataset with only complete cases in explanatory variables
Prop_dam_CC<-Int[complete.cases(Int[,c(3,4)]),]
Prop_dam_CC$Tree_corr<-Prop_dam_CC$Trees-mean(Prop_dam_CC$Trees)
str(Prop_dam_CC)
M1<-lme(log(Vol)~Tree_corr*Region,random=~1|Study,Prop_dam_CC)
summary(M1)

#look at diagnostic plots
plot(M1)
plot(log(Prop_dam_CC$Vol),(predict(M1)))
abline(a=0,b=1)
plot((Prop_dam_CC$Vol),exp(predict(M1)))
abline(a=0,b=1)

#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7)

#create predictions from model
Tree_pred1<-data.frame(Site=Int$Site,Tree_corr=Int$Trees-mean(Prop_dam_CC$Trees),Region=Int$Region)
Tree_pred1<-Tree_pred1[complete.cases(Tree_pred1),]
Tree_pred1$Pred<-exp(predict(avgm,newdata=Tree_pred1,level=0))

#merge onto dataset of other data
M_Damage<-merge(Int,Tree_pred1,by="Site",all=T)
#put in predicted volume values where there are none where number of treesextracted is <20 per ha
for (i in 1:nrow(M_Damage)){
  M_Damage$Vol2[i]<-ifelse(is.na(M_Damage$Vol[i])&M_Damage$Trees[i]<20,M_Damage$Pred[i],M_Damage$Vol[i])
}

summary(M_Damage)

ggplot(M_Damage,aes(x=Trees,y=Vol2,colour=Region.x))+geom_point(size=4)


######################################################################
#model to predict the number of trees removed from the volume logged##
######################################################################

M1<-lme(log(Trees)~Vol*Region,random=~1|Study,Prop_dam_CC)
summary(M1)
#diagnostic plots
plot(M1)
qqnorm(M1)

#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7)


#create dataframe for use with predictions
Tree_pred1<-data.frame(Site=Int$Site,Vol=Int$Vol,Region=Int$Region)
Tree_pred1<-Tree_pred1[complete.cases(Tree_pred1),]
Tree_pred1$Pred<-exp(predict(avgm,newdata=Tree_pred1,level=0))
M_Damage2<-merge(M_Damage,Tree_pred1,by="Site",all=T)

str(M_Damage2)

#put in the number of trees extracted per ha where there are none
#and volume is within the range of the model
#so that there is no crazy extrapolation
for (i in 1:nrow(M_Damage2)){
  M_Damage2$Trees2[i]<-ifelse(is.na(M_Damage2$Trees[i])&M_Damage2$Vol.x[i]<110,M_Damage2$Pred.y[i],M_Damage2$Trees[i])
}

str(M_Damage2)

#now remove columns which are no use
keeps <- c("Study","Site","Vol.x","Vol2","Trees","Trees2","Region.x")
M_Damage3<-M_Damage2[keeps]
#and remove duplicate rows
M_Damage4<-M_Damage3[!duplicated(M_Damage3), ]

#save this data for use in other analyses
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
write.csv(M_Damage4,"Log_intens.csv",row.names=F)


###################################################################
#merge this data with that from AGB, Richness and damage analyses##
###################################################################
#first AGB
head(M_Damage4)
AGB_final<-merge(AGB,M_Damage4,by="Site")[-c(19,20,22)]
write.csv(AGB_final,"AGB_intens.csv",row.names=F)

head(Rich)

#richness
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Rich_final<-merge(Rich,M_Damage4,by.x="Site.ID",by.y="Site")[,-c(5:7,15,18,19,20,21,23)]
str(Rich_final)
colnames(Rich_final)<-c("Site","Study","Age","Method","M_UL","V_UL","SS_UL","M_L","V_L","SS_L","VarT","Rare","Region","Vol","Trees","Region2")
Rich_final2<-Rich_final[!duplicated(Rich_final), ]
write.csv(Rich_final2,"Rich_intens.csv",row.names=F)

#now damage
Damage_final<-merge(Prop_dam,M_Damage4,by.x="Site_ID",by.y="Site")[,-c(3,5:15,18,23,24,26)]
str(Damage_final)
colnames(Damage_final)<-c("Site","Study","Method","Prop_dam","ID","Region","N_logged","Plot","Notes","Vol2","Trees2","Region2")
write.csv(Damage_final,"Dam_intens.csv",row.names=F)
