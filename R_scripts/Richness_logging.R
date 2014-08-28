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
setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Richness<-read.csv("Logged_Richness.csv")
colnames(Richness)<-c("Study","Site","Age","Method","Vol","Type","MU","VU","SSU","ML","VL","SSL","Vtype","Tax","Rar","Region")

#calculate SDs
#unlogged
Richness$SDU<-ifelse(Richness$Vtype=="SE",Richness$VU*sqrt(Richness$SSU),Richness$VU)
Richness$SDU<-ifelse(Richness$Vtype=="CI",(Richness$VU/1.96)*sqrt(Richness$SSU),Richness$SDU)
#logged
Richness$SDL<-ifelse(Richness$Vtype=="SE",Richness$VL*sqrt(Richness$SSL),Richness$VL)
Richness$SDL<-ifelse(Richness$Vtype=="CI",(Richness$VL/1.96)*sqrt(Richness$SSL),Richness$SDL)

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
Richness$Rare<-ifelse(Richness$Rar=="Not rarefied","NR",NA)
Richness$Rare<-ifelse(Richness$Rar=="Rarefied - Area","R",Richness$Rare)
Richness$Rare<-ifelse(Richness$Rar=="Rarefied - Individuals","R",Richness$Rare)
Richness$Rare<-as.factor(Richness$Rare)

#calculate the log ratio
ROM<-escalc(data=Richness,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
ROM$Age<-ifelse(is.na(ROM$Age),mean(ROM$Age,na.rm=T),ROM$Age)

setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.csv(ROM,"Richness_studies.csv")

###############################################
#summary analysis##############################
###############################################

#random effects meta-analysis
ROM.ma<-rma.mv(yi,vi,random=list(~1|Rare),method="REML",data=ROM)
summary(ROM.ma)
exp(-.0914)-1
exp(-0.1778)-1
exp(-0.0049)-1

#forest plot of this
theme_set(theme_bw(base_size=25))
forrest_data2<-rbind(data.frame(ES=ROM.ma$yi,SE=sqrt(ROM.ma$vi),Type="Site",Study=Richness$Site,Method=Richness$Method),data.frame(ES=ROM.ma$b,SE=ROM.ma$se,Type="Summary",Study="Summary",Method="Summary"))
forrest_data2$Study2<-factor(forrest_data2$Study, levels=rev(levels(forrest_data2$Study)) )
levels(forrest_data2$Study2)
plot1<-ggplot(data=forrest_data2,aes(x=Study2,y=exp(ES)-1,ymax=exp(ES+(1.96*SE))-1,ymin=exp(ES-(1.96*SE))-1,size=(1/SE)/20,colour=factor(Type)))+geom_pointrange(shape=15)
plot2<-plot1+coord_flip()+geom_hline(aes(x=0), lty=2,size=1)
plot3<-plot2+xlab("Study")+ylab("Proportional change in \n species richness following logging")+scale_colour_manual(values=c("grey","black"))
plot4_rich<-plot3+theme(legend.position="none")+scale_size_continuous(range=c(1,2.5))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+theme(legend.position="none")+ labs(title="Species richness")
setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Forest_Richness_poster.jpeg",height=12,width=12,dpi=1200)

##################################################################
#model to test for region and method differences##################
##################################################################

Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|Rare,~1|Tax,~1|Age),method="ML",data=ROM)
Model1<-rma.mv(yi,vi,mods=~Method-1+Region,random=list(~1|Rare,~1|Tax,~1|Age),method="ML",data=ROM)
Model2<-rma.mv(yi,vi,mods=~Method-1,random=list(~1|Rare,~1|Tax,~1|Age),method="ML",data=ROM)
Model3<-rma.mv(yi,vi,mods=~Region-1,random=list(~1|Rare,~1|Tax,~1|Age),method="ML",data=ROM)
Model4<-rma.mv(yi,vi,mods=~Tax-1,random=list(~1|Rare,~1|Age),method="ML",data=ROM)

ggplot(ROM,aes(x=Tax,y=yi,colour=Tax))+geom_jitter()

#work out model AICc
Model_AIC<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],
                             Model2$fit.stats$ML[5],Model3$fit.stats$ML[5]))


Model_AIC$Vars<-c("Null","Method + Region","Method","Region")

#calculate pseudo r squared
Null_sigma1<-Model0$sigma2[1]   
Null_sigma2<-Model0$sigma2[2]

(sum(Model0$sigma2) - sum(Model1$sigma2)) / sum(Model0$sigma2)

Model_AIC$sigma1<-c(Model0$sigma2[1],Model1$sigma2[1],
              Model2$sigma2[1],Model3$sigma2[1])
Model_AIC$sigma2<-c(Model0$sigma2[2],Model1$sigma2[2],
              Model2$sigma2[2],Model3$sigma2[2])
Model_AIC$sumsigma1<-ifelse(Null_sigma1-Model_AIC$sigma1<0,0,Null_sigma1-Model_AIC$sigma1)
Model_AIC$sumsigma2<-ifelse(Null_sigma2-Model_AIC$sigma2<0,0,Null_sigma2-Model_AIC$sigma2)

Model_AIC$R_squared<-1-((Null_sigma1+Null_sigma2)-(Model_AIC$sumsigma1+Model_AIC$sumsigma2))/(Null_sigma1+Null_sigma2)

#reorder from lowest to highest
Model_AIC<-Model_AIC[order(Model_AIC$AIC),]
#calculate AICc delta
Model_AIC$delta<-Model_AIC$AIC-Model_AIC$AIC[1]
#drop last models with delta >7
AIC_sel<-subset(Model_AIC,delta<=7)
#calculate the relative likelihood of model
AIC_sel$rel_lik<-exp((AIC_sel$AIC[1]-AIC_sel$AIC)/2)
#calculate the AICc weight
AIC_sel$weight<-AIC_sel$rel_lik/(sum(AIC_sel$rel_lik))
#write this table to directory
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.table(AIC_sel,file="Rich_without_vol.csv",sep=",")

#dummy coding for variable importance
AIC_sel$Region<-c(0,0,1,1)
AIC_sel$Method<-c(0,1,0,1)

Importance<-c(sum(AIC_sel$Region*AIC_sel$weight),sum(AIC_sel$Method*AIC_sel$weight))
Importance_vals<-data.frame(Variable=c("Region","Logging method"),Importance=Importance)
#export importance values
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.table(Importance_vals,file="Imp_Rich_no_vol.csv",sep=",")




############################################################
#Analysis for studies including volume######################
#accounting for age and survey method diffs#################
############################################################
ROM_vol<-subset(ROM,!is.na(Vol))
ROM_vol<-subset(ROM,Vol!=-9999)
Rich_vol<-subset(Rich_vol,Rare!="NR")

#calculate the log ratio
ROM_vol<-escalc(data=ROM_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)

#models of richness change including volume
Model0_Vol<-rma.mv(yi,vi,mods=~1,random=~list(~1|Rare),method="ML",data=ROM_vol)
Model1_Vol<-rma.mv(yi,vi,mods=~Vol,random=~list(~1|Rare),method="ML",data=ROM_vol)
Model2_Vol<-rma.mv(yi,vi,mods=~Method,random=~list(~1|Rare),method="ML",data=ROM_vol)

plot(ROM_vol$Vol,predict(Model1_Vol)$pred)

#work out model AICc
Model_AIC<-data.frame(AICc=c(Model0_Vol$fit.stats$ML[5],Model1_Vol$fit.stats$ML[5],Model2_Vol$fit.stats$ML[5]))

Model_AIC$Vars<-c("Null","Volume",
                   "Method")

ggplot(ROM_vol,aes(x=Vol,y=yi))+geom_point()


pred.null<-predict(Model0_Vol)$pred
pred <- predict(Model2_Vol)$pred

1 - sum((ROM_vol$yi-pred)^2) / sum((ROM_vol$yi-pred.null)^2)

#calculate pseudo r squared
Null_sigma1<-Model0_Vol$sigma2[1]   
Null_sigma2<-Model0_Vol$sigma2[2]

(sum(Model0_Vol$sigma2) - sum(Model1_Vol$sigma2)) / sum(Model0_Vol$sigma2)

Model_AIC$sigma1<-c(Model0$sigma2[1],Model1$sigma2[1],
              Model2$sigma2[1],Model3$sigma2[1])
Model_AIC$sigma2<-c(Model0$sigma2[2],Model1$sigma2[2],
              Model2$sigma2[2],Model3$sigma2[2])
Model_AIC$sumsigma1<-ifelse(Null_sigma1-Model_AIC$sigma1<0,0,Null_sigma1-Model_AIC$sigma1)
Model_AIC$sumsigma2<-ifelse(Null_sigma2-Model_AIC$sigma2<0,0,Null_sigma2-Model_AIC$sigma2)

Model_AIC$R_squared<-1-((Null_sigma1+Null_sigma2)-(Model_AIC$sumsigma1+Model_AIC$sumsigma2))/(Null_sigma1+Null_sigma2)

#reorder from lowest to highest
Model_AIC<-Model_AIC[order(Model_AIC$AIC),]
#calculate AICc delta
Model_AIC$delta<-Model_AIC$AIC-Model_AIC$AIC[1]

#drop last models with delta >7
AIC_sel<-subset(Model_AIC,delta<=7)
#calculate the relative likelihood of model
AIC_sel$rel_lik<-exp((AIC_sel$AIC[1]-AIC_sel$AIC)/2)
#calculate the AICc weight
AIC_sel$weight<-AIC_sel$rel_lik/(sum(AIC_sel$rel_lik))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.table(AIC_sel,file="Rich_with_vol.csv",sep=",")

#dummy coding for variable importance
AIC_sel$Vol<-c(1,0,0)
AIC_sel$Method<-c(0,0,1)

sum(AIC_sel$Vol*AIC_sel$weight)
sum(AIC_sel$Method*AIC_sel$weight)

#write this table to directory
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.table(AIC_sel,file="Rich_with_vol.csv",sep=",")

#re-do model with REML
Model1_reml<-rma.mv(yi,vi,mods=~Vol,random=~(1|as.factor(Age)),method="REML",data=ROM_vol)
summary(Model1_reml)

#create dataframe for predictions

all<-data.frame(yi=ROM_vol$yi,vi=ROM_vol$vi,Vol=ROM_vol$Vol,Method=ROM_vol$Method)
summary(ROM_vol$Vol)

Vol<-seq(18,92.4,length.out=500)
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
vol_plot4<-vol_plot3+ylab("Proportional change in /nrichness following logging")+geom_point(shape=16,aes(x=Vol,y=exp(yi)-1,colour=Method,size=(1/vi)*5),alpha=0.5)

vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=new_preds,aes(x=Vol,y=exp(preds)-1),size=2)+geom_hline(y=0,lty=2,size=1)
vol_plot5
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot7<-vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
rich_vol_plot<-vol_plot7+geom_line(data=new_preds,aes(y=exp(ci.lb)-1,x=Vol),lty=3,size=1)+geom_line(data=new_preds,aes(y=exp(ci.ub)-1,x=Vol),lty=3,size=1)+xlim(18,100)+ylim(-.30,.20)+ guides(colour = guide_legend(override.aes = list(size=18)))

setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_volume.jpeg",height=12,width=12,dpi=1200)
