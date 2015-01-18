#############################################################
#script to calculate scaling coeffients between different####
#measures of logging intensity and damage####################
#############################################################

rm(list=ls())

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

################################################################
#model of volume as a function of number of trees logged per ha#
################################################################

#create dataset with only complete cases in explanatory variables
Prop_dam_CC<-Prop_dam[complete.cases(Prop_dam[,c(7,8)]),]
Prop_dam_CC$Tree_corr<-Prop_dam_CC$Tree_ex_ha-mean(Prop_dam_CC$Tree_ex_ha)
str(Prop_dam_CC)
Prop_dam_CC<-subset(Prop_dam_CC,Tree_ex_ha<20)
M1<-lm(Vol_log~Tree_corr*Region-1,Prop_dam_CC,na.action="na.fail")
summary(M1)

#look at diagnostic plots
plot(M1)
plot(Prop_dam_CC$Vol_log,(predict(M1)))
abline(a=0,b=1)

#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7)

#create predictions from model
tapply(Prop_dam_CC$Tree_ex_ha,Prop_dam_CC$Region,summary)

Tree_pred1<-data.frame(Trees=c(seq(0.4,7.4,0.1),seq(0.2,16,0.1),seq(4.7,10.7,0.1)),Region=c(rep("Africa",71),rep("Americas",159),rep("Asia",61)))
Tree_pred1$Tree_corr<-Tree_pred1$Trees-mean(Prop_dam_CC$Tree_ex_ha)
Tree_pred1<-Tree_pred1[complete.cases(Tree_pred1),]
Tree_pred1$Pred<-(predict(M1,newdata=Tree_pred1,level=0))

#plot this prediction
theme_set(theme_bw(base_size=12))
Scaling_plot1<-ggplot(Prop_dam_CC,aes(x=Tree_ex_ha,y=Vol_log,colour=Region))+geom_point(size=3,shape=1)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Scaling_plot2<-Scaling_plot1+xlab(expression(paste("Number of trees extracted ",ha^-1)))+ylab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
Scaling_plot2+geom_line(data=Tree_pred1,aes(x=Trees,y=Pred,group=Region,colour=Region),size=2)
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Volume_trees.jpeg",height=6,width=8,dpi=1200)


#merge onto dataset of other data
Tree_pred2<-data.frame(Site_ID=Prop_dam$Site_ID,Tree_corr=Prop_dam$Tree_ex_ha-mean(Prop_dam_CC$Tree_ex_ha),Region=Prop_dam$Region)
Tree_pred2<-Tree_pred2[complete.cases(Tree_pred2),]
Tree_pred2$Pred<-(predict(M1,newdata=Tree_pred2,level=0))
M_Damage<-merge(Prop_dam,Tree_pred2,by="Site_ID",all=T)
head(M_Damage)
#put in predicted volume values where there are none where number of treesextracted is <25 per ha
M_Damage$Vol_log2<-NULL
for (i in 1:nrow(M_Damage)){
  M_Damage$Vol_log2[i]<-ifelse(is.na(M_Damage$Vol_log[i])&M_Damage$Tree_ex_ha[i]<25,M_Damage$Pred[i],M_Damage$Vol_log[i])
}

######################################################################
#model to predict the number of trees removed from the volume logged##
######################################################################

M1<-lme(Tree_ex_ha~Vol_log*Region,random=~1|Study,Prop_dam_CC)

r.squaredGLMM(M1)

#diagnostic plots
plot(M1)
qqnorm(M1)

#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7)


#create dataframe for use with predictions
Tree_pred1<-data.frame(Site_ID=Prop_dam$Site_ID,Vol_log=Prop_dam$Vol_log,Region=Prop_dam$Region)
Tree_pred1<-Tree_pred1[complete.cases(Tree_pred1),]
Tree_pred1$Pred<-exp(predict(M1,newdata=Tree_pred1,level=0))
M_Damage2<-merge(M_Damage,Tree_pred1,by="Site_ID",all=T)


#put in the number of trees extracted per ha where there are none
for (i in 1:nrow(M_Damage2)){
  M_Damage2$Tree_ex_ha2[i]<-ifelse(is.na(M_Damage2$Tree_ex_ha[i])&M_Damage2$Vol_log.x[i]<150,M_Damage2$Pred.y[i],M_Damage2$Tree_ex_ha[i])
}


################################################################################################
#model to predict the proportion of stems damaged from the number of trees damaged per hectare##
################################################################################################
Prop_dam_CC2<-subset(Prop_dam_CC,!is.na(Dam_ha))
Prop_dam_CC2$log_ha<-log(Prop_dam_CC2$Dam_ha)
M1<-lm(qlogis(Prop_dam)~log_ha,data=Prop_dam_CC2)
summary(M1)
log_ha<-data.frame(log_ha=log(Prop_dam_CC2$Dam_ha))
Prop_dam_CC2$Pred<-plogis(predict(M1,newdata=log_ha))
M_Damage3<-merge(M_Damage2,Prop_dam_CC2,all=T)

M_Damage3$Prop_dam2<-NULL
for (i in 1:nrow(M_Damage3)){
  M_Damage3$Prop_dam2[i]<-ifelse(is.na(M_Damage3$Prop_dam[i]),M_Damage3$Pred[i],M_Damage3$Prop_dam[i])
}


Preds<-data.frame(Dam_ha=exp(log_ha),Pred=plogis(predict(M1,newdata=log_ha)))
colnames(Preds)<-c("Dam_ha","Pred")
head(Preds)

#make a plot of the relationship
Scaling_plot1<-ggplot(M_Damage3,aes(x=Dam_ha,y=Prop_dam))+geom_point(size=3,shape=1)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Scaling_plot2<-Scaling_plot1+xlab(expression(paste("Number of trees damaged ",ha^-1)))+ylab("Proportion of residual trees damaged")
Scaling_plot2+geom_line(data=Preds,aes(x=Dam_ha,y=Pred))
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Volume_trees.jpeg",height=6,width=8,dpi=1200)


#now remove columns which are no use
keeps <- c("Study","Site_ID","Vol_log2","Tree_ex_ha2","Method","Prop_dam","Prop_dam2","Prop_seve_dam","Region")
M_Damage4<-M_Damage3[keeps]
M_Damage4<-subset(M_Damage4,Vol_log2<200)
M_Damage4<-subset(M_Damage4,Study!="Guariguata et al 2009")
M_Damage4$Prop_dam2<-ifelse(M_Damage4$Study=="Whitman et al 1997",0.048,M_Damage4$Prop_dam2)

setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
write.csv(M_Damage4,"Dam_intens.csv",row.names=F)

