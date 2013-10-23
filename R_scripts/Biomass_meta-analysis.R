#######################################################################################
#Script for meta-analysis of changes in aboveground biomass with logging###############
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:23/10/2013

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

#try new mixed effects meta-analysis

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

###########################################################################################
####Carry out meta-analysis################################################################
###########################################################################################

#calculate effect sizes
#log ratio
ROM<-escalc(data=AGB,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)

#runs a random effects meta-analysis for log ratio data
ME_summary<-rma.mv(yi,vi,random=~(1|ID),method="REML",data=ROM)
summary(ME_summary)

#forrest plot of this
theme_set(theme_bw(base_size=10))
forrest_data<-rbind(data.frame(ES=ME_summary$yi,SE=sqrt(ME_summary$vi),Type="Site",Study=AGB$Site),data.frame(ES=ME_summary$b,SE=ME_summary$se,Type="Summary",Study="Summary"))
forrest_data$Study2<-factor(forrest_data$Study, levels=rev(levels(forrest_data$Study)) )
levels(forrest_data$Study2)
plot1<-ggplot(data=forrest_data,aes(x=Study2,y=exp(ES)-1,ymax=exp(ES+(1.96*SE))-1,ymin=exp(ES-(1.96*SE))-1,size=(1/SE)/4,colour=factor(Type)))+geom_pointrange(shape=15)
plot2<-plot1+coord_flip()+geom_hline(aes(x=0), lty=2,size=1)
plot3<-plot2+xlab("Site")+ylab("Proportional change")+scale_colour_manual(values=c("grey","black"))
plot3+theme(legend.position="none")+scale_size_continuous(range=c(0.5,1.5))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Forest_BM.pdf",height=6,width=6,dpi=1200)

##################################################################################
##analysis of how effects ########################################################
#on biomass vary by region, method, age and number of times logged################
##################################################################################

#compare models with region and method as covariates
model0<-rma.mv(yi,vi,mods=~1,random=~(1|ID),method="ML",data=ROM)
model1<-rma.mv(yi,vi,mods=~Region-1,random=~(1|ID),method="ML",data=ROM)
model2<-rma.mv(yi,vi,mods=~Method-1,random=~(1|ID),method="ML",data=ROM)
model3<-rma.mv(yi,vi,mods=~Age,random=~(1|ID),method="ML",data=ROM)
model4<-rma.mv(yi,vi,mods=~N_Logged2,random=~(1|ID),method="ML",data=ROM)
model5<-rma.mv(yi,vi,mods=~N_Logged2-1+Region-1,random=~(1|ID),method="ML",data=ROM)
model6<-rma.mv(yi,vi,mods=~N_Logged2-1+Method-1,random=~(1|ID),method="ML",data=ROM)
model7<-rma.mv(yi,vi,mods=~N_Logged2-1+Age,random=~(1|ID),method="ML",data=ROM)

#look at which model is the best fit
model0$fit.stats$ML[5]
Model_AICc<-data.frame(AICc=c(model0$fit.stats$ML[5],model1$fit.stats$ML[5],model2$fit.stats$ML[5],model3$fit.stats$ML[5],model4$fit.stats$ML[5],model5$fit.stats$ML[5],model6$fit.stats$ML[5],model7$fit.stats$ML[5]))
Model_AICc$model<-c("Null","Model1","Model2","Model3","Model4","Model5","Model6","Model7")

#calculate AICc delta
Model_AICc$delta<-Model_AICc$AIC-min(Model_AICc$AICc)

#calculate pseudo r squared
Null_dev<-logLik(model0)[1]
Model_AICc$logLik<-c(logLik(model0)[1],logLik(model1)[1],logLik(model2)[1],logLik(model3)[1],logLik(model4)[1],logLik(model5)[1],logLik(model6)[1],logLik(model7)[1])
Model_AICc$pseudo_R2<-(Null_dev-Model_AICc$logLik)/Null_dev

#drop models with delta >7
AICc_sel<-subset(Model_AICc,delta<=7)
#calculate the relative likelihood of model
AICc_sel$rel_lik<-exp((AICc_sel$AIC[1]-AICc_sel$AIC)/2)
#calculate the AICc weight
AICc_sel$weight<-AICc_sel$rel_lik/(sum(AICc_sel$rel_lik))
#reorder sorting
AICc_sel<-AICc_sel[order(AICc_sel$AIC),]
#put in model parameters
AICc_sel$Vars<-c("No. times logged & Region","No. times logged & Logging method","No. times logged","No. times logged & Age")
rownames(AICc_sel)<-NULL
#output this as a table
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.table(AICc_sel,file="AGB_without_vol.csv",sep=",")

#now we calculate the relative importance of each 
#first dummy coding for variables
AICc_sel$N_logged<-c(1,1,1,1)
AICc_sel$Region<-c(1,0,0,0)
AICc_sel$Method<-c(0,1,0,0)
AICc_sel$Age<-c(0,0,0,1)

Importance_val<-c(sum(AICc_sel$N_logged*AICc_sel$weight),sum(AICc_sel$Region*AICc_sel$weight),
sum(AICc_sel$Method*AICc_sel$weight),sum(AICc_sel$Age*AICc_sel$weight))

Importance_AGB_no_vol<-data.frame(Variable=c("No. of times logged","Region","Method","Age"),Importance=Importance_val)
write.csv(Importance_AGB_no_vol,file="Imp_AGB_no_vol.csv")


#it looks like the one including region and number of times logged is best
#so we recalculate the parameter estimates using REML
REML_model5<-rma(yi,vi,mods=~Region-1+N_Logged2,method="REML",data=ROM)
summary(REML_model5)

coef(summary(REML_model5))

Preds_AGB<-predict(REML_model5,addx=T)

Preds_AGB2<-data.frame(Preds=Preds_AGB[1],SE=Preds_AGB[2],ci.lb=Preds_AGB[3],ci.ub=Preds_AGB[4])

Preds_unique<-unique(Preds_AGB2)
Preds_unique$Region<-c("SE Asia & Australasia","Americas","Africa","SE Asia & Australasia")
Preds_unique$N_Logged<-c(1,1,1,2)

#do model averaging
# once logged
Mod_av1_1<-coef(summary(model5))[1,1]*AICc_sel$weight[1]
Mod_av2_1<-coef(summary(model6))[1,1]*AICc_sel$weight[2]
Mod_av3_1<-coef(summary(model4))[1,1]*AICc_sel$weight[3]
Mod_av4_1<-coef(summary(model7))[1,1]*AICc_sel$weight[4]

Once_logged<-(Mod_av1_1+Mod_av2_1+Mod_av3_1+Mod_av4_1)

Mod_av1_2<-coef(summary(model5))[2,1]*AICc_sel$weight[1]
Mod_av2_2<-coef(summary(model6))[2,1]*AICc_sel$weight[2]
Mod_av3_2<-coef(summary(model4))[2,1]*AICc_sel$weight[3]
Mod_av4_2<-coef(summary(model7))[2,1]*AICc_sel$weight[4]

Twice_logged<-(Mod_av1_2+Mod_av2_2+Mod_av3_2+Mod_av4_2)



#now we can plot this - the effect of region and number of times logged
theme_set(theme_bw(base_size=10))
a<-ggplot(data=Preds_unique,aes(y=exp(pred)-1,ymin=exp(ci.lb)-1,ymax=exp(ci.ub)-1,x=Region,colour=factor(N_Logged)))+geom_pointrange(size=1,position=position_dodge(width = 0.2))
b<-a+geom_hline(x=0,lty=2,size=1)+ylab("Proportional change after logging")+xlab("Regions")
c<-b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1,colour="black",fill=NA))
c+scale_color_brewer(name="Number of \ntimes logged",palette="Set1")
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("AGB_Region_N_logged.pdf",height=4,width=6,dpi=1200)

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
Model0<-rma.mv(yi,vi,mods=~1,random=~(1|ID),method="ML",data=ROM2)
Model1<-rma.mv(yi,vi,mods=~I(Vol/MU)+Method,random=~(1|ID),method="ML",data=ROM2)
Model2<-rma.mv(yi,vi,mods=~~Vol+Method,random=~(1|ID),method="ML",data=ROM2)
Model3<-rma.mv(yi,vi,mods=~Vol,random=~(1|ID),method="ML",data=ROM2)
Model4<-rma.mv(yi,vi,mods=~I(Vol/MU),random=~(1|ID),method="ML",data=ROM2)
Model5<-rma.mv(yi,vi,mods=~Vol+MU,random=~(1|ID),method="ML",data=ROM2)
Model6<-rma.mv(yi,vi,mods=~MU,random=~(1|ID),method="ML",data=ROM2)
Model7<-rma.mv(yi,vi,mods=~MU+Method,random=~(1|ID),method="ML",data=ROM2)
Model8<-rma.mv(yi,vi,mods=~Vol+I(Vol^2),random=~(1|ID),method="ML",data=ROM2)
Model9<-rma.mv(yi,vi,mods=~Vol+I(Vol^2)+MU,random=~(1|ID),method="ML",data=ROM2)
Model10<-rma.mv(yi,vi,mods=~Vol*Method,random=~(1|ID),method="ML",data=ROM2)
plot(fitted(Model10),resid(Model10))

Model_AICc<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],Model3$fit.stats$ML[5],Model4$fit.stats$ML[5],Model5$fit.stats$ML[5],Model6$fit.stats$ML[5],Model7$fit.stats$ML[5],Model8$fit.stats$ML[5],Model9$fit.stats$ML[5]))
Model_AICc$model<-c("Null","Model1","Model2","Model3","Model4","Model5","Model6","Model7","Model8","Model9")
#calculate AICc delta
Model_AICc$delta<-Model_AICc$AIC-min(Model_AICc$AICc)
#calculate pseudo r squared for each model
Null_dev<-logLik(Model0)[1]
Model_AICc$logLik<-c(logLik(Model0)[1],logLik(Model1)[1],logLik(Model2)[1],logLik(Model3)[1],logLik(Model4)[1],logLik(Model5)[1],logLik(Model6)[1],logLik(Model7)[1],logLik(Model8)[1],logLik(Model9)[1])
Model_AICc$pseudo_R2<-(Null_dev-Model_AICc$logLik)/Null_dev

#drop last models with delta >7
AICc_sel<-subset(Model_AICc,delta<=7)
#calculate the relative likelihood of model
AICc_sel$rel_lik<-exp((AICc_sel$AICc[1]-AICc_sel$AICc)/2)
#calculate the AICc weight
AICc_sel$weight<-AICc_sel$rel_lik/(sum(AICc_sel$rel_lik))
#reorder sorting
AICc_sel<-AICc_sel[order(AICc_sel$AICc),]

setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.csv(AICc_sel,file="AGB_vol.csv")

#re-do model with REML
Model8_reml<-rma.mv(yi,vi,mods=~Vol+I(Vol^2),random=~(1|ID),method="REML",data=ROM2)

#create dataframe for predictions

all<-data.frame(yi=ROM2$yi,vi=ROM2$vi,Vol=ROM2$Vol,Method=ROM2$Method,Logged=ROM2$N_Logged,MU=ROM2$MU)
all$preds<-(predict(Model8_reml))$pred
all$ci.lb<-(predict(Model8_reml))$ci.lb
all$ci.ub<-(predict(Model8_reml))$ci.ub
plot(all$yi,all$preds)
plot(fitted(Model8_reml),resid(Model8_reml))

#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
AGB_ax<-expression(paste("Unlogged biomass (Mg",ha^-1,")"))
theme_set(theme_bw(base_size=10))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2+theme(legend.position="none")
vol_plot3
vol_plot4<-vol_plot3+ylab("Proportional change in biomass following logging")+geom_point(shape=1,aes(x=Vol,y=exp(yi)-1,colour=Method,size=1/vi))

vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=all,aes(x=Vol,y=exp(preds)-1))+geom_hline(y=0,lty=2,size=1)
vol_plot5
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_volume.pdf",height=4,width=6,dpi=1200)

#plot of volume logged vs unlogged biomass
theme_set(theme_bw(base_size=10))
AGB_Vol<-ggplot(AGB_vol,aes(x=MU,y=Vol,colour=Region))+geom_point(size=4)
AGB_Vol2<-AGB_Vol+ylab((expression(paste("Volume of wood logged (",m^3,ha^-1,")"))))
AGB_vol3<-AGB_Vol2+xlab("Mean unlogged biomass in reference plots")
AGB_vol3+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Volume_AGB.pdf",height=4,width=6,dpi=1200)