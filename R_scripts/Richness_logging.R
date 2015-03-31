#######################################################################################
#Script for meta-analysis of changes in species richness with logging##################
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:13/03/2013

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(metafor)
library(MuMIn)
library(boot)

#import data
Richness<-read.csv("Data/Rich_intens.csv")
head(Richness)


#calculate SDs
#unlogged
Richness$SDU<-ifelse(Richness$VarT=="SE",Richness$V_UL*sqrt(Richness$SS_UL),Richness$V_UL)
Richness$SDU<-ifelse(Richness$VarT=="CI",(Richness$V_UL/1.96)*sqrt(Richness$SS_UL),Richness$SDU)
#logged
Richness$SDL<-ifelse(Richness$VarT=="SE",Richness$V_L*sqrt(Richness$SS_L),Richness$V_L)
Richness$SDL<-ifelse(Richness$VarT=="CI",(Richness$V_L/1.96)*sqrt(Richness$SS_L),Richness$SDL)

#impute missing standard deviation values
#based on coefficient of variation
#following ideas in Koricheva et al 2013
Richness2<-subset(Richness,Richness$SDU>0)
(Richness2)
Imp_U<-(sum(Richness2$SDU))/(sum(Richness2$M_UL))
Imp_L<-(sum(Richness2$SDL))/(sum(Richness2$M_L))
Richness$SDU<-ifelse(Richness$SDU<0,Richness$M_UL*Imp_U,Richness$SDU)
Richness$SDL<-ifelse(Richness$SDL<0,Richness$M_L*Imp_L,Richness$SDL)


#change rarefied column
head(Richness)
Richness$Rare2<-ifelse(Richness$Rare=="Not rarefied","NR",NA)
Richness$Rare2<-ifelse(Richness$Rare=="Rarefied - Area","R",Richness$Rare2)
Richness$Rare2<-ifelse(Richness$Rare=="Rarefied - Individuals","R",Richness$Rare2)
Richness$Rare2<-as.factor(Richness$Rare2)

#calculate the log ratio
ROM<-escalc(data=Richness,measure="ROM",m2i=M_UL,sd2i=SDU,n2i=SS_UL,m1i=M_L,sd1i=SDL,n1i=SS_L,append=T)
ROM$Age<-ifelse(is.na(ROM$Age),mean(ROM$Age,na.rm=T),ROM$Age)

write.csv(ROM,"Data/Richness_studies.csv")


############################################################
#Analysis for studies including volume######################
#accounting for age and survey method diffs#################
############################################################
ROM$Vol<-ifelse(is.na(ROM$Vol),mean(ROM$Vol),ROM$Vol)

ROM_vol<-subset(ROM,!is.na(Vol))
ROM_vol<-subset(ROM,Vol!=-9999)
Rich_vol<-subset(ROM_vol,!is.na(vi))
sum(Rich_vol$SS_UL)
sum(Rich_vol$SS_L)


#models of richness change including volume

#standardise volume and age using Zuurs methods
Rich_vol$Vol_std<-(Rich_vol$Vol-mean(Rich_vol$Vol))/sd(Rich_vol$Vol)
Rich_vol$Age_std<-(Rich_vol$Age-mean(Rich_vol$Age))/sd(Rich_vol$Age)


head(Rich_vol)
Model0_Vol<-rma.mv(yi,vi,mods=~1,random=list(~1|Study, ~ 1 | Rare,~1|as.factor(M_UL)),method="ML",data=Rich_vol)
Model1_Vol<-rma.mv(yi,vi,mods=~Vol_std,random=list(~1|Study, ~ 1 | Rare,~1|as.factor(M_UL)),method="ML",data=Rich_vol)
Model2_Vol<-rma.mv(yi,vi,mods=~Vol_std+I(Vol_std^2),random=list(~1|Study, ~ 1 | Rare,~1|as.factor(M_UL)),method="ML",data=Rich_vol)
Model3_Vol<-rma.mv(yi,vi,mods=~Method,random=list(~1|Study, ~ 1 | Rare,~1|as.factor(M_UL)),method="ML",data=Rich_vol)
Model4_Vol<-rma.mv(yi,vi,mods=~Method*Vol,random=list(~1|Study, ~ 1 | Rare,~1|as.factor(M_UL)),method="ML",data=Rich_vol)
Model5_Vol<-rma.mv(yi,vi,mods=~Vol_std*Age_std,random=list(~1|Study, ~ 1 | Rare,~1|as.factor(M_UL)),method="ML",data=Rich_vol)
Model6_Vol<-rma.mv(yi,vi,mods=~Age_std,random=list(~1|Study, ~ 1 | Rare,~1|as.factor(M_UL)),method="ML",data=Rich_vol)
Model7_Vol<-rma.mv(yi,vi,mods=~Method+Vol,random=list(~1|Study, ~ 1 | Rare,~1|as.factor(M_UL)),method="ML",data=Rich_vol)


ggplot(Rich_vol,aes(x=Vol_std,y=yi,size=1/vi))+geom_point()

#work out model AICc
Model_AIC<-data.frame(AICc=c(Model0_Vol$fit.stats$ML[5],Model1_Vol$fit.stats$ML[5],Model2_Vol$fit.stats$ML[5],Model3_Vol$fit.stats$ML[5],Model4_Vol$fit.stats$ML[5],Model5_Vol$fit.stats$ML[5],Model6_Vol$fit.stats$ML[5]))

Model_AIC$Vars<-c("Null","Volume","Volume+Volume^2",
                   "Method","Volume*Method","Volume*Age","Age")

(sum(Model0_Vol$sigma2) - sum(Model4_Vol$sigma2)) / sum(Model0_Vol$sigma2)

#calculate r squared
str(Model0_Vol)

Null_sigma<-sum(Model0_Vol$sigma2)
Model_AIC$sigma<-c(sum(Model0_Vol$sigma2),sum(Model1_Vol$sigma2),
                   sum(Model2_Vol$sigma2),sum(Model3_Vol$sigma2),
                   sum(Model4_Vol$sigma2),sum(Model5_Vol$sigma2))

c(Null_sigma-Model_AIC$sigma)/Null_sigma


#reorder from lowest to highest
Model_AIC<-Model_AIC[order(Model_AIC$AICc),]
#calculate AICc delta
Model_AIC$delta<-Model_AIC$AICc-Model_AIC$AICc[1]

#calculate the relative likelihood of model
Model_AIC$rel_lik<-exp((Model_AIC$AICc[1]-Model_AIC$AICc)/2)
#calculate the AICc weight
Model_AIC$weight<-Model_AIC$rel_lik/(sum(Model_AIC$rel_lik))



setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.table(Model_AIC,file="Rich_with_vol.csv",sep=",")

#re-do model with REML
Model1_reml<-rma.mv(yi,vi,mods=~Vol+Age,random=list(~1|Study, ~ 1 | Rare2),method="REML",data=Rich_vol)
summary(Model1_reml)

#create dataframe for predictions

all<-data.frame(yi=ROM_vol$yi,vi=ROM_vol$vi,Vol=ROM_vol$Vol,Method=ROM_vol$Method,Age=ROM_vol$Age)
summary(ROM_vol$Vol)

summary(ROM_vol)

newdat<-expand.grid(Vol=seq(5,118,length.out=500),Age=c(0,10,20))
preds<-predict.rma(Model1_reml,newmods=cbind(c(5,50,100),c(0,0,0)),addx=T)
?predict.rma

new_preds<-data.frame(preds=preds$pred,ci.lb=preds$ci.lb,ci.ub=preds$ci.ub,Vol=newdat$Vol,Age=newdat$Age)

#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
theme_set(theme_bw(base_size=25))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2
vol_plot3
vol_plot4<-vol_plot3+ylab("Proportional change in tree species richness following logging")+geom_point(shape=16,aes(x=Vol,y=exp(yi)-1,colour=Age,size=(1/vi)*5))

vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=new_preds,aes(x=Vol,y=exp(preds)-1,colour=Age,group=as.factor(Age)),size=2)+geom_hline(y=0,lty=2,size=1)
vol_plot5
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot7<-vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
rich_vol_plot<-vol_plot7+geom_line(data=new_preds,aes(y=exp(ci.lb)-1,x=Vol),lty=3,size=1)+geom_line(data=new_preds,aes(y=exp(ci.ub)-1,x=Vol),lty=3,size=1)+theme(legend.position="none")+scale_colour_brewer(palette="Set1")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("SR_volume.png",height=12,width=12,dpi=400)


#create funnel plot with residuals

plot(resid(Model1_reml),1/sqrt(Rich_vol$vi))
