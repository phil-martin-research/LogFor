#######################################################################################
#Script for producing figures of how biomass & richness vary with logging intensity#####
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:22/11/2013

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(metafor)
library(GGally)
library(multcomp)


###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################

setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")

AGB<-read.csv("Logged_AGB.csv")


#recalculate SDs
#unlogged
AGB$SDU<-ifelse(AGB$VarT=="SE",AGB$VU*sqrt(AGB$SSU),AGB$VU)
AGB$SDU<-ifelse(AGB$VarT=="CI",(AGB$VU/1.96)*sqrt(AGB$SSU),AGB$SDU)
#logged
AGB$SDL<-ifelse(AGB$VarT=="SE",AGB$VL*sqrt(AGB$SSL),AGB$VL)
AGB$SDL<-ifelse(AGB$VarT=="CI",(AGB$VL/1.96)*sqrt(AGB$SSL),AGB$SDL)
head(AGB)
AGB$N_Logged2<-as.factor(AGB$N_Logged)
AGB<-subset(AGB,Age!=18)

#subset data to remove those without vol
AGB_vol<-subset(AGB,Vol>0)

#import data
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Richness<-read.csv("Logged_Richness.csv")

#import data
head(Richness)
colnames(Richness)<-c("Study","Site","Age","Method","Vol","Type","MU","VU","SSU","ML","VL","SSL","Vtype","Tax","Rar","Region")

#recalculate SDs

#unlogged
Richness$SDU<-ifelse(Richness$Vtype=="SE",Richness$VU*sqrt(Richness$SSU),Richness$VU)
Richness$SDU<-ifelse(Richness$Vtype=="CI",(Richness$VU/1.96)*sqrt(Richness$SSU),Richness$SDU)

#logged
Richness$SDL<-ifelse(Richness$Vtype=="SE",Richness$VL*sqrt(Richness$SSL),Richness$VL)
Richness$SDL<-ifelse(Richness$Vtype=="CI",(Richness$VL/1.96)*sqrt(Richness$SSL),Richness$SDL)
head(Richness)

#impute missing standard deviation values
Richness2<-subset(Richness,Richness$SDU>0)
head(Richness2)
Imp_U<-(sum(Richness2$SDU))/(sum(Richness2$MU))
Imp_L<-(sum(Richness2$SDL))/(sum(Richness2$ML))

Richness$SDU<-ifelse(Richness$SDU<0,Richness$MU*Imp_U,Richness$SDU)
Richness$SDL<-ifelse(Richness$SDL<0,Richness$ML*Imp_L,Richness$SDL)

#change rarefied column
Richness$Rare<-ifelse(Richness$Rar=="Not rarefied","NR",NA)
Richness$Rare<-ifelse(Richness$Rar=="Rarefied - Area","R",Richness$Rare)
Richness$Rare<-ifelse(Richness$Rar=="Rarefied - Individuals","R",Richness$Rare)
Richness$Rare<-as.factor(Richness$Rare)


#################################################################################
#analysis of the effect of volume on post logging biomass########################
#################################################################################

#log ratio effect size calculation for results with volume
ROM2<-escalc(data=AGB_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
head(ROM2)

#run best model
Model8_reml<-rma.mv(yi,vi,mods=~Vol+I(Vol^2),random=~(1|ID),method="REML",data=ROM2)

#create dataframe for predictions

all<-data.frame(yi=ROM2$yi,vi=ROM2$vi,Vol=ROM2$Vol,Method=ROM2$Method,Metric="Aboveground biomass")

Vol<-seq(8.11,179,length.out=500)
preds<-predict.rma(Model8_reml,newmods=cbind(Vol,Vol^2),addx=T)
head(preds)

new_preds_AGB<-data.frame(preds=preds$pred,ci.lb=preds$ci.lb,ci.ub=preds$ci.ub,Vol=Vol,Metric="Aboveground biomass")

#####do the same for richness#######
Rich_vol<-subset(Richness,Vol!=-9999)
Rich_vol<-subset(Rich_vol,Rare!="NR")

#calculate the log ratio
ROM_vol<-escalc(data=Rich_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)

#best model
Model1_reml<-rma.uni(yi,vi,mods=~Vol,method="REML",data=ROM_vol)

#create dataframe for predictions

all2<-data.frame(yi=ROM_vol$yi,vi=ROM_vol$vi,Vol=ROM_vol$Vol,Method=ROM_vol$Method,Metric="Species richness")

Vol2<-seq(18,92.4,length.out=500)
preds<-predict.rma(Model1_reml,newmods=Vol2,addx=T)

new_preds_Rich<-data.frame(preds=preds$pred,ci.lb=preds$ci.lb,ci.ub=preds$ci.ub,Vol=Vol2,Metric="Species richness")
plot(new_preds_Rich$Vol,new_preds_Rich$preds)

#combine data and predicitons
head(new_preds_Rich)
new_preds<-rbind(new_preds_AGB,new_preds_Rich)
head(new_preds)
all3<-rbind(all,all2)


#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
theme_set(theme_bw(base_size=10))
vol_plot<-ggplot(data=all3)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2+theme(legend.position="none")
vol_plot3
vol_plot4<-vol_plot3+ylab("Proportional change following logging")+geom_point(shape=1,aes(x=Vol,y=exp(yi)-1,colour=Method,size=4))


vol_plot5<-vol_plot4+geom_line(data=new_preds,aes(x=Vol,y=exp(preds)-1))+geom_hline(y=0,lty=2,size=1)
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot7<-vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
vol_plot8<-vol_plot7+geom_line(data=new_preds,aes(y=exp(ci.lb)-1,x=Vol),lty=3)+geom_line(data=new_preds,aes(y=exp(ci.ub)-1,x=Vol),lty=3)
vol_plot8+facet_wrap(~Metric,scales="free")
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Biomas_richness_volume.pdf",height=5,width=10,dpi=1200)

