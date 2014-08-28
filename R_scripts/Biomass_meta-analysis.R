#######################################################################################
#Script for meta-analysis of changes in aboveground biomass with logging###############
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:12/03/2013

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(metafor)
library(GGally)
library(multcomp)
library(MuMIn)
library(plyr)
library(regress.r)


###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################

#try new mixed effects meta-analysis

setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
AGB<-read.csv("Logged_AGB_rev.csv")

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
AGB_NoSafe<-subset(AGB,N_Logged<2)

#subset data to remove those without vol
AGB_vol<-subset(AGB,Vol>0)


###########################################################################################
####Carry out meta-analysis for reviewers to look at plot size and variability#############
###########################################################################################

#calculate effect sizes
#log ratio
ROM<-escalc(data=AGB,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
ROM

setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
write.csv(ROM,"AGB_studies.csv")

#grain size
theme_set(theme_bw(base_size=12))
a<-ggplot(ROM,aes(x=Plot_size,y=1/vi))+geom_point(size=4,shape=1)+scale_x_log10()+scale_y_log10()
b<-a+xlab("grain size (ha)")+ylab("relative weight contributed\nto meta-analysis")
b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("AGB_grain.pdf",height=4,width=6,dpi=1200)

#extent size
theme_set(theme_bw(base_size=12))
a<-ggplot(ROM,aes(x=Plot_size*(SSU+SSL),y=1/vi))+geom_point()+scale_x_log10()
a
b<-a+xlab("study extent (ha)")+ylab("relative weight contributed\nto meta-analysis")
b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("AGB_grain.pdf",height=4,width=6,dpi=1200)

#now calculate cv based on plot size
Rom_plot<-subset(ROM,Plot_size>0&SDU>0)

#first for unlogged plots
Plot_MU<-lm(log(I(SDU/MU))~(Plot_size),data=Rom_plot)
summary(Plot_MU)
plot(Plot_MU)
plot((Rom_plot$Plot_size),-1/sqrt(Rom_plot$SDU/Rom_plot$MU),)
plot((Rom_plot$Plot_size),(Rom_plot$SDL/Rom_plot$MU),)

predict.plot((I(SDU/MU)) ~ Plot_size, Rom_plot)


#now for logged plots
Plot_ML<-lm(I(SDL/ML)~(Plot_size)+I(Plot_size^2),data=Rom_plot)
summary(Plot_ML)
plot(Plot_ML)

Preds1<-(data.frame(PS=(Rom_plot$Plot_size),Pred=predict(Plot_ML)))
head(Preds1)
Preds1<-Preds1[with(Preds1, order(PS)), ]
plot((Rom_plot$Plot_size),(Rom_plot$SDL/Rom_plot$ML),)
lines(Preds1)

.13904+(-.1685*log(6.25))
AGB$SDU

#now create loop to fill in gaps
head(AGB)

for (i in 1:nrow(AGB)){
  AGB$SDU<-ifelse(is.na(AGB$SDU[i])&AGB$Plot_size[i]>1,AGB$MU*0.15,AGB$SDU[i])
  AGB$SDU<-ifelse(is.na(AGB$SDU[i])&AGB$Plot_size[i]<1,AGB$MU*0.25,AGB$SDU[i])
  AGB$SDL<-ifelse(is.na(AGB$SDL[i])&AGB$Plot_size[i]>1,AGB$ML*0.2,AGB$SDU[i])
  AGB$SDL<-ifelse(is.na(AGB$SDL[i])&AGB$Plot_size[i]<1,AGB$ML*0.5,AGB$SDU[i])
}


#now re-do the meta-analysis

ROM_add<-escalc(data=AGB,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
ROM_add_nosafe<-subset(ROM_add,N_Logged==1)
#runs a random effects meta-analysis for log ratio data
ME_summary<-rma.mv(yi,vi,random=list(~1|ID,~1|Age),method="REML",data=ROM_add)
ME_summary2<-rma.mv(yi,vi,random=list(~1|ID,~1|Age),method="REML",data=ROM_add_nosafe)
ME_summary2<-rma.mv(yi,vi,random=list(~1|ID),method="REML",data=ROM_add)

summary(ME_summary)
summary(ME_summary2)
exp(-0.5691)-1
exp(-0.3036)-1
exp(-0.8357)-1

#forrest plot of this
library(ggplot2)
theme_set(theme_bw(base_size=25))
forrest_data<-rbind(data.frame(ES=ME_summary$yi,SE=sqrt(ME_summary$vi),Type="Site",Study=AGB$Site,Site=AGB$Site),data.frame(ES=ME_summary$b,SE=ME_summary$se,Type="Summary",Study="Summary",Site="Summary"))
forrest_data$Study2<-factor(forrest_data$Study, levels=rev(levels(forrest_data$Study)) )
forrest_data$Study3<-factor(forrest_data$Site, levels=rev(levels(forrest_data$Site)) )
levels(forrest_data$Study2)
forrest_data$ord<-seq(1:nrow(forrest_data))
forrest_data$Method<-c(AGB$Method,"Summary")
forrest_data$Method<-c(AGB$Method,"Summary")
forrest_data$caption<-c(paste(AGB$Site," (",AGB$Plot_size,", ",AGB$SSL,")",sep=""),"Summary")
forrest_data$caption<-as.factor(forrest_data$caption)
forrest_data$caption<-(factor(forrest_data$caption,levels(forrest_data$caption)[c(41,40:1)]))
plot1<-ggplot(data=forrest_data,aes(x=reorder(Study,-ord),y=exp(ES)-1,ymax=exp(ES+(1.96*SE))-1,ymin=exp(ES-(1.96*SE))-1,colour=factor(Type)))+geom_pointrange(shape=15,size=2)
plot2<-plot1+coord_flip()+geom_hline(aes(x=0), lty=2,size=1)
plot3<-plot2+xlab("Site")+ylab("Proportional change in biomass \nfollowing logging")+scale_colour_manual(values=c("grey","black"))
plot4<-plot3+theme(legend.position="none")+scale_size_continuous(range=c(0.25,1))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ labs(title="Aboveground biomass")

setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Forest_BM.jpeg",height=12,width=12,dpi=1200)

##################################################################################
##analysis of how effects ########################################################
#on biomass vary by region, method, age and number of times logged################
##################################################################################

#compare models with region and method as covariates
model0<-rma.mv(yi,vi,mods=~1,random=list(~1|ID,~1|Age),method="ML",data=ROM_add)
model1<-rma.mv(yi,vi,mods=~Region-1,random=list(~1|ID,~1|Age),method="ML",data=ROM_add)
model2<-rma.mv(yi,vi,mods=~Method-1,random=list(~1|ID,~1|Age),method="ML",data=ROM_add)
model3<-rma.mv(yi,vi,mods=~Method+Region-1,random=list(~1|ID,~1|Age),method="ML",data=ROM_add)

#look at which model is the best fit
model0$fit.stats$ML[5]
Model_AICc<-data.frame(AICc=c(model0$fit.stats$ML[5],model1$fit.stats$ML[5],model2$fit.stats$ML[5],model3$fit.stats$ML[5]))
Model_AICc$model<-c("Null","Model1","Model2","Model3")

#calculate AICc delta
Model_AICc$delta<-Model_AICc$AIC-min(Model_AICc$AICc)

#calculate pseudo r squared
Null_sigma<-(sum(model0$sigma2))
 
Model_AICc$sigma<-c(sum(model0$sigma2),sum(model1$sigma2),
              sum(model2$sigma2),sum(model3$sigma2))

Model_AICc$R_squared<-c(Null_sigma-Model_AICc$sigma)/Null_sigma

#drop models with delta >7
AICc_sel<-subset(Model_AICc,delta<=7)
#calculate the relative likelihood of model
AICc_sel$rel_lik<-exp((AICc_sel$AIC[1]-AICc_sel$AIC)/2)
#calculate the AICc weight
AICc_sel$weight<-AICc_sel$rel_lik/(sum(AICc_sel$rel_lik))
#reorder sorting
AICc_sel<-AICc_sel[order(AICc_sel$AIC),]
#put in model parameters
AICc_sel$Vars<-c("Region + logging method","Region","Logging method","Null_model")
rownames(AICc_sel)<-NULL
#output this as a table
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.table(AICc_sel,file="AGB_without_vol.csv",sep=",")

#now we calculate the relative importance of each 
#first dummy coding for variables
AICc_sel$Region<-c(1,1,0,0)
AICc_sel$Method<-c(1,0,1,0)

Importance_val<-c(sum(AICc_sel$Region*AICc_sel$weight),sum(AICc_sel$Method*AICc_sel$weight))

Importance_AGB_no_vol<-data.frame(Variable=c("Region","Method"),Importance=Importance_val)
write.csv(Importance_AGB_no_vol,file="Imp_AGB_no_vol.csv")

#it looks like the one including region and method is best
#so we recalculate the parameter estimates using REML
REML_model5<-rma.mv(yi,vi,mods=~Region+Method-1,random=list(~1|ID,~1|Age),method="REML",data=ROM_add)
summary(REML_model5)

Preds_AGB<-predict(REML_model5,addx=T)

Preds_AGB2<-data.frame(pred=Preds_AGB$pred,SE=Preds_AGB$se,ci.lb=Preds_AGB$ci.lb,ci.ub=Preds_AGB$ci.ub)

Preds_unique<-unique(Preds_AGB2)
Preds_unique$Region<-c("SE Asia","Americas","Africa")
Preds_unique$Method<-rep(c("Conventional","RIL"),each=3)


#now we can plot this - the effect of region and number of times logged
theme_set(theme_bw(base_size=10))
a<-ggplot(data=Preds_unique,aes(y=exp(pred)-1,ymin=exp(ci.lb)-1,ymax=exp(ci.ub)-1,x=Region,group=Method,colour=Method))+geom_pointrange(size=1,position=position_dodge(width = 0.2))
b<-a+geom_hline(x=0,lty=2,size=1)+ylab("Proportional change after logging")+xlab("Regions")
c<-b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1,colour="black",fill=NA))
c+scale_color_brewer(name="Logging method",palette="Set1")
setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("AGB_Region.pdf",height=4,width=6,dpi=1200)

#output of how logging intensity varies by region
a<-ggplot(data=AGB_vol,aes(x=Region,y=Vol,fill=N_Logged2))+geom_boxplot(size=0.2)
b<-a+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1,colour="black",fill=NA))
c<-b+ylab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_fill_brewer(name="Number of \ntimes logged",palette="Set1")
c
ggsave("Region_vol.pdf",height=4,width=6,dpi=1200)


#################################################################################
#analysis of the effect of volume on post logging biomass########################
#################################################################################

#log ratio effect size calculation for results with volume
ROM2<-escalc(data=AGB_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
head(ROM2)

#different models relating volume and method to post logging change
Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|ID,~1|Age),method="ML",data=ROM2)
Model1<-rma.mv(yi,vi,mods=~~Vol+Method,random=list(~1|ID,~1|Age),method="ML",data=ROM2)
Model2<-rma.mv(yi,vi,mods=~Vol,random=list(~1|ID,~1|Age),method="ML",data=ROM2)
Model3<-rma.mv(yi,vi,mods=~I(Vol/MU),random=list(~1|ID,~1|Age),method="ML",data=ROM2)
Model4<-rma.mv(yi,vi,mods=~Method,random=list(~1|ID,~1|Age),method="ML",data=ROM2)
Model5<-rma.mv(yi,vi,mods=~Vol+I(Vol^2),random=list(~1|ID,~1|Age),method="ML",data=ROM2)
Model6<-rma.mv(yi,vi,mods=~Vol*Method,random=list(~1|ID,~1|Age),method="ML",data=ROM2)

Model_AICc<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],Model3$fit.stats$ML[5],Model4$fit.stats$ML[5],Model5$fit.stats$ML[5],Model6$fit.stats$ML[5]))
Model_AICc$model<-c("Null","Model1","Model2","Model3","Model4","Model5","Model6")
#calculate AICc delta
Model_AICc$delta<-Model_AICc$AICc-min(Model_AICc$AICc)
#calculate pseudo r squared for each model
Null_sigma<-(sum(Model0$sigma2))
 
Model_AICc$sigma<-c(sum(Model0$sigma2),sum(Model1$sigma2),
              sum(Model2$sigma2),sum(Model3$sigma2),
              sum(Model4$sigma2),sum(Model5$sigma2),
              sum(Model6$sigma2))

Model_AICc$R_squared<-c(Null_sigma-Model_AICc$sigma)/Null_sigma

#drop last models with delta >7
AICc_sel<-Model_AICc
#calculate the relative likelihood of model
AICc_sel$rel_lik<-exp((AICc_sel$AICc[1]-AICc_sel$AICc)/2)
AICc_sel$rel_lik<-round(AICc_sel$rel_lik,2)
#calculate the AICc weight
AICc_sel$weight<-AICc_sel$rel_lik/(sum(AICc_sel$rel_lik))
AICc_sel$weight<-round(AICc_sel$weight,3)
#reorder sorting
AICc_sel<-AICc_sel[order(AICc_sel$AICc),]
#put in variables
AICc_sel$vars<-c("Vol+Vol^2","Vol","Vol*Method","Vol+Method","Vol/MU","Null","Method")

setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.csv(AICc_sel,file="AGB_vol.csv")

#re-do model with REML
Model8_reml<-rma.mv(yi,vi,mods=~Vol+I(Vol^2),random=list(~1|ID,~1|Age),method="REML",data=ROM2)

#create dataframe for predictions
all<-data.frame(yi=ROM2$yi,vi=ROM2$vi,Vol=ROM2$Vol,Method=ROM2$Method,Logged=ROM2$N_Logged,MU=ROM2$MU)
all$preds<-(predict(Model8_reml,level=0))$pred
all$ci.lb<-(predict(Model8_reml,level=0))$ci.lb
all$ci.ub<-(predict(Model8_reml,level=0))$ci.ub
Vol<-seq(8.11,179,length.out=500)
preds<-predict.rma(Model8_reml,newmods=cbind(Vol,Vol^2),addx=T)
head(preds)
preds
new_preds<-data.frame(preds=preds$pred,ci.lb=preds$ci.lb,ci.ub=preds$ci.ub,Vol=Vol)

#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
AGB_ax<-expression(paste("Unlogged biomass (Mg",ha^-1,")"))
theme_set(theme_bw(base_size=25))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2
vol_plot4<-vol_plot3+ylab("Proportional change in \nbiomass following logging")
vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=new_preds,aes(x=Vol,y=exp(preds)-1),size=3)+geom_hline(y=0,lty=2,size=2)
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot7<-vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
biommass_vol_plot<-vol_plot7+geom_line(data=new_preds,aes(y=exp(ci.lb)-1,x=Vol),lty=3,size=2)+geom_line(data=new_preds,aes(y=exp(ci.ub)-1,x=Vol),lty=3,size=2)
biommass_vol_plot+geom_point(shape=16,aes(x=Vol,y=exp(yi)-1,colour=Method,size=1/vi),alpha=0.5)+scale_size_continuous(range=c(8,16))+ guides(colour = guide_legend(override.aes = list(size=18)))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_volume2.jpeg",height=12,width=12,dpi=1200)

#plot of volume logged vs unlogged biomass
theme_set(theme_bw(base_size=10))
AGB_Vol<-ggplot(AGB_vol,aes(x=MU,y=Vol,colour=Region))+geom_point(size=4)
AGB_Vol2<-AGB_Vol+ylab((expression(paste("Volume of wood logged (",m^3,ha^-1,")"))))
AGB_vol3<-AGB_Vol2+xlab("Mean unlogged biomass in reference plots")
AGB_vol3+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Volume_AGB.pdf",height=4,width=6,dpi=1200)