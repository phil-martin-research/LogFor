#############################################################
#script to predict proportion of trees damaged using volume##
#############################################################

#first load packages
library(nlme)
library(ggplot2)
library(MuMIn)
library(lme4)
library(nlme)

#clear objects
rm(list=ls())

#import data
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Dam<-read.csv("Dam_intens.csv")
head(Dam)
colnames(Dam)[3]<-"Vol2"
Dam<-Dam[complete.cases(Dam[,7]),]

#######################################################
#now predict proportion of trees damaged using volume##
#######################################################

#create column for Volume squared and log volume
Dam$Vol_sq<-Dam$Vol2^2
Dam$Vol_log<-log(Dam$Vol2)
Dam$Vol_log2<-(log(Dam$Vol2))^2

nrow(Dam)

#test for the effects of volume logged, non-linear volume logged (squared and log), differences in method
#and a null model
M1<-lmer(qlogis(Prop_dam2)~Vol2*Method+Vol_sq*Method+(1|Study),Dam,na.action="na.fail")
M2<-lmer(qlogis(Prop_dam2)~Vol2+Vol_sq+(1|Study),Dam,na.action="na.fail")
M3<-lmer(qlogis(Prop_dam2)~Vol2+(1|Study),Dam,na.action="na.fail")
M4<-lmer(qlogis(Prop_dam2)~Vol_log*Method+Vol_log2*Method+(1|Study),Dam,na.action="na.fail")
M5<-lmer(qlogis(Prop_dam2)~Vol_log*Method+(1|Study),Dam,na.action="na.fail")
M6<-lmer(qlogis(Prop_dam2)~Method+(1|Study),Dam,na.action="na.fail")
M0<-lmer(qlogis(Prop_dam2)~1+(1|Study),Dam,na.action="na.fail")
AICc(M1,M2,M3,M4,M5,M6,M0)
summary(M1)
All_mods<-list(M1,M2,M3,M4,M5,M6,M0)

#diagnostic plots

#model averaging
ms1<-mod.sel(object=All_mods,rank="AICc",fit=T,trace=T,subset=dc(Vol2,Vol_sq))
ms1$r2<-c(r.squaredGLMM(M5)[1],r.squaredGLMM(M6)[1],r.squaredGLMM(M3)[1],r.squaredGLMM(M0)[1],r.squaredGLMM(M4)[1],r.squaredGLMM(M2)[1],r.squaredGLMM(M1)[1])
#write this table to csv
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.csv(ms1,"Tree_damage_models.csv")


delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7,se.fit=T)


#predict from model averaged stuff
tapply(Dam$Vol2,Dam$Method,summary)
length(seq(5,104,1))


Damage_pred<-data.frame(Vol2=c(seq(11,153,1),seq(5,104,1)),Method=c(rep("Conventional",143),rep("RIL",100)))
Damage_pred$Vol_sq<-Damage_pred$Vol2^2
Damage_pred$Vol_log<-log(Damage_pred$Vol2)
Damage_pred$Vol_log2<-Damage_pred$Vol_log^2


Damage_pred2<-predict(avgm,newdata=Damage_pred,level=0,se.fit=T)
Damage_pred2<-cbind(Damage_pred,Damage_pred2)
head(Damage_pred2)

#plot these results
theme_set(theme_bw(base_size=12))
Dam_1<-ggplot(Dam,aes(x=Vol2,y=Prop_dam,colour=Method))+geom_point(shape=16,size=3)+geom_line(data=Damage_pred2,aes(x=Vol2,y=plogis(fit),colour=Method),size=2)
Dam_2<-Dam_1+geom_line(data=Damage_pred2,aes(x=Vol2,y=plogis(fit+(1.96*se.fit)),colour=Method),lty=2)+geom_line(data=Damage_pred2,aes(x=Vol2,y=plogis(fit-(1.96*se.fit)),colour=Method),lty=2)
Dam_3<-Dam_2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Dam_3+xlab("Volume logged per hectare")+ylab("Proportion of residual tree stems damaged")+scale_colour_brewer(palette = "Set1")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_damaged_vol.jpeg",height=6,width=8,dpi=1200)


