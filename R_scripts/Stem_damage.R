#############################################################
#script to calculate scaling coeffients between different####
#measures of logging intensity and damage####################
#############################################################

rm(list=ls())

#first load packages
library(ggplot2)
library(MuMIn)
library(lme4)
library(plyr)
library(reshape2)

#load data
Prop_dam<-read.csv("Data/prop_damage_edit2.csv")
head(Prop_dam)

length(which(!is.na(Prop_dam$Damage_prop)))

colnames(Prop_dam)<-c("Study","Site_ID","ID","Method","Vol","Tree_ex","Dam_ha","Prop_dam","Region","Plot","Replicates")

################################################################
#model of volume as a function of number of trees logged per ha#
################################################################

#create dataset with only complete cases in explanatory variables
Prop_dam_CC<-Prop_dam[complete.cases(Prop_dam[,c(5,6)]),]
Prop_dam_CC$Tree_corr<-(Prop_dam_CC$Tree_ex-mean(Prop_dam_CC$Tree_ex,na.rm = T))/sd(Prop_dam_CC$Tree_ex,na.rm = T)

ggplot(Prop_dam_CC,aes(x=Tree_ex,y=Vol,colour=Region))+geom_point()+facet_wrap(~Study)

M1<-lmer(Vol~Tree_corr*Region-1+(Tree_corr|Study),Prop_dam_CC,na.action="na.fail")
summary(M1)

#look at diagnostic plots
plot(M1)
plot(Prop_dam_CC$Tree_corr,(predict(M1,re.form=NA)))

#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7)
summary(avgm)

#create predictions from model
tapply(Prop_dam_CC$Tree_ex,Prop_dam_CC$Region,summary)

Tree_pred1<-data.frame(Trees=c(seq(0.4,7.4,length.out = 500),seq(1.85,16,length.out = 500),seq(3,10.7,length.out = 500)),Region=c(rep("Africa",500),rep("Americas",500),rep("Asia",500)))
Tree_pred1$Tree_corr<-(Tree_pred1$Trees-mean(Prop_dam_CC$Tree_ex,na.rm = T))/sd(Prop_dam_CC$Tree_ex,na.rm = T)
Tree_pred1$Pred<-(predict(avgm,newdata=Tree_pred1,re.form=NA))

#plot this prediction
theme_set(theme_bw(base_size=12))
Scaling_plot1<-ggplot(Prop_dam_CC,aes(x=Tree_ex,y=Vol,colour=Region))+geom_point(size=3,shape=1,alpha=0.5)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Scaling_plot2<-Scaling_plot1+xlab(expression(paste("Number of trees extracted ",ha^-1)))+ylab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
Scaling_plot2+geom_line(data=Tree_pred1,aes(x=Trees,y=Pred,group=Region,colour=Region),size=2,alpha=0.5)
ggsave("Figures/Volume_trees.png",height=6,width=8,dpi=400)


#merge onto dataset of other data
Tree_pred2<-data.frame(ID=Prop_dam$ID,Tree_corr=(Prop_dam$Tree_ex-mean(Prop_dam_CC$Tree_ex,na.rm = T))/sd(Prop_dam_CC$Tree_ex,na.rm = T),Region=Prop_dam$Region)
Tree_pred2<-Tree_pred2[complete.cases(Tree_pred2),]
Tree_pred2$Pred<-(predict(avgm,newdata=Tree_pred2,re.form=NA))
M_Damage<-merge(Prop_dam,Tree_pred2,by="ID",all=T)
#put in predicted volume values where there are none where number of treesextracted is <25 per ha
M_Damage$Vol2<-NULL
for (i in 1:nrow(M_Damage)){
  M_Damage$Vol2[i]<-ifelse(is.na(M_Damage$Vol[i]),M_Damage$Pred[i],M_Damage$Vol[i])
}


################################################################################################
#model to predict the proportion of stems damaged from the number of trees damaged per hectare##
################################################################################################
Prop_dam_CC2<-subset(Prop_dam,!is.na(Dam_ha)&!is.na(Prop_dam))

M1<-lmer(qlogis(Prop_dam)~(Dam_ha)*Region-1+(Dam_ha|Study),data=Prop_dam_CC2,na.action="na.fail")
#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7)
summary(avgm)

#add predictions to data frame for use in further analyses
Prop_dam_CC2$Pred<-plogis(predict(avgm,Prop_dam_CC2,re.form=NA))
M_Damage2<-merge(M_Damage,Prop_dam_CC2,by="ID",all=T)
head(M_Damage2)
M_Damage2$Prop_dam2<-NULL
for (i in 1:nrow(M_Damage2)){
  M_Damage2$Prop_dam2[i]<-ifelse(is.na(M_Damage2$Prop_dam.x[i]),M_Damage2$Pred.y[i],M_Damage2$Prop_dam.x[i])
}

#create predictions for figures
ddply(Prop_dam_CC2,.(Region),summarise,Dam_min=min(Dam_ha),Dam_max=max(Dam_ha))
Dam_ha<-rbind(data.frame(Dam_ha=seq(3.5,31,length.out = 50),Region="Africa"),
              data.frame(Dam_ha=seq(10.075,191.4,length.out = 50),Region="Americas"),
              data.frame(Dam_ha=seq(18.7,189.9,length.out = 50),Region="Asia"))
Dam_ha$Pred<-plogis(predict(avgm,newdata=Dam_ha,re.form=NA))

#make predictions for plot

qplot(Prop_dam_CC2$Dam_ha,qlogis(Prop_dam_CC2$Prop_dam),colour=Prop_dam_CC2$Region)

#make a plot of the relationship
Scaling_plot1<-ggplot(Prop_dam_CC2,aes(x=Dam_ha,y=Prop_dam,colour=Region))+geom_point(size=3,shape=1)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Scaling_plot2<-Scaling_plot1+xlab(expression(paste("Number of trees damaged ",ha^-1)))+ylab("Proportion of residual trees damaged")
Scaling_plot2+geom_line(data=Dam_ha,aes(x=Dam_ha,y=Pred))
ggsave("Figures/Damage_coeff.png",height=6,width=8,dpi=400)

str(M_Damage2)
#now remove columns which are no use
keeps <- c("Study.x","Site_ID.x","ID","Vol2","Method.x","Prop_dam2","Region")
M_Damage3<-M_Damage2[keeps]
write.csv(M_Damage3,"Data/Dam_intens.csv",row.names=F)

