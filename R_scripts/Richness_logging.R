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

#import data
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Richness<-read.csv("Rich_intens.csv")
head(Richness)


#calculate SDs
#unlogged
Richness$SDU<-ifelse(Richness$VarT=="SE",Richness$V_UL*sqrt(Richness$SS_UL),Richness$V_UL)
Richness$SDU<-ifelse(Richness$VarT=="CI",(Richness$V_UL/1.96)*sqrt(Richness$SS_UL),Richness$SD_UL)
#logged
Richness$SDL<-ifelse(Richness$VarT=="SE",Richness$V_L*sqrt(Richness$SS_L),Richness$V_L)
Richness$SDL<-ifelse(Richness$VarT=="CI",(Richness$V_L/1.96)*sqrt(Richness$SS_L),Richness$SD_L)

#impute missing standard deviation values
#based on coefficient of variation
#following ideas in Koricheva et al 2013
Richness2<-subset(Richness,Richness$SDU>0)
(Richness2)
Imp_U<-(sum(Richness2$SDU))/(sum(Richness2$MU))
Imp_L<-(sum(Richness2$SDL))/(sum(Richness2$ML))
Richness$SDU<-ifelse(Richness$SDU<0,Richness$MU*Imp_U,Richness$SDU)
Richness$SDL<-ifelse(Richness$SDL<0,Richness$ML*Imp_L,Richness$SDL)

#change rarefied column
head(Richness)
Richness$Rare2<-ifelse(Richness$Rare=="Not rarefied","NR",NA)
Richness$Rare2<-ifelse(Richness$Rare=="Rarefied - Area","R",Richness$Rare2)
Richness$Rare2<-ifelse(Richness$Rare=="Rarefied - Individuals","R",Richness$Rare2)
Richness$Rare2<-as.factor(Richness$Rare2)

#calculate the log ratio
ROM<-escalc(data=Richness,measure="ROM",m2i=M_UL,sd2i=SDU,n2i=SS_UL,m1i=M_L,sd1i=SDL,n1i=SS_L,append=T)
ROM$Age<-ifelse(is.na(ROM$Age),mean(ROM$Age,na.rm=T),ROM$Age)

setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
write.csv(ROM,"Richness_studies.csv")


############################################################
#Analysis for studies including volume######################
#accounting for age and survey method diffs#################
############################################################
ROM_vol<-subset(ROM,!is.na(Vol))
ROM_vol<-subset(ROM,Vol!=-9999)
Rich_vol<-subset(ROM_vol,!is.na(vi))


#models of richness change including volume
Model0_Vol<-rma.mv(yi,vi,mods=~1,random=~list(~1|Rare),method="ML",data=Rich_vol)
Model1_Vol<-rma.mv(yi,vi,mods=~Vol,random=~list(~1|Rare),method="ML",data=Rich_vol)
Model2_Vol<-rma.mv(yi,vi,mods=~Vol+I(Vol^2),random=~list(~1|Rare),method="ML",data=Rich_vol)
Model3_Vol<-rma.mv(yi,vi,mods=~Method,random=~list(~1|Rare),method="ML",data=Rich_vol)
Model4_Vol<-rma.mv(yi,vi,mods=~Method+Vol,random=~list(~1|Rare),method="ML",data=Rich_vol)


#work out model AICc
Model_AIC<-data.frame(AICc=c(Model0_Vol$fit.stats$ML[5],Model1_Vol$fit.stats$ML[5],Model2_Vol$fit.stats$ML[5],Model3_Vol$fit.stats$ML[5],Model4_Vol$fit.stats$ML[5]))

Model_AIC$Vars<-c("Null","Volume","Volume+Volume^2",
                   "Method","Volume+Method")

#calculate r squared
str(Model0_Vol)

Null_sigma<-sum(Model0_Vol$sigma2)
Model_AIC$sigma<-c(sum(Model0_Vol$sigma2),sum(Model1_Vol$sigma2),
                   sum(Model2_Vol$sigma2),sum(Model3_Vol$sigma2),
                   sum(Model4_Vol$sigma2))

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
Model1_reml<-rma.mv(yi,vi,mods=~Vol,random=~(1|as.factor(Age)),method="REML",data=Rich_vol)
summary(Model1_reml)

#create dataframe for predictions

all<-data.frame(yi=ROM_vol$yi,vi=ROM_vol$vi,Vol=ROM_vol$Vol,Method=ROM_vol$Method)
summary(ROM_vol$Vol)

summary(ROM_vol)

Vol<-seq(5,118,length.out=500)
preds<-predict.rma(Model1_reml,newmods=Vol,addx=T)
head(preds)

new_preds<-data.frame(preds=preds$pred,ci.lb=preds$ci.lb,ci.ub=preds$ci.ub,Vol=Vol)

#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
theme_set(theme_bw(base_size=25))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2
vol_plot3
vol_plot4<-vol_plot3+ylab("Proportional change in tree species richness following logging")+geom_point(shape=16,aes(x=Vol,y=exp(yi)-1,colour=Method,size=(1/vi)*5))

vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=new_preds,aes(x=Vol,y=exp(preds)-1),size=2)+geom_hline(y=0,lty=2,size=1)
vol_plot5
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot7<-vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
rich_vol_plot<-vol_plot7+geom_line(data=new_preds,aes(y=exp(ci.lb)-1,x=Vol),lty=3,size=1)+geom_line(data=new_preds,aes(y=exp(ci.ub)-1,x=Vol),lty=3,size=1)+theme(legend.position="none")+scale_colour_brewer(palette="Set1")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("SR_volume.png",height=12,width=12,dpi=400)


#create funnel plot with residuals

plot(resid(Model1_reml),1/sqrt(Rich_vol$vi))
