#######################################################################################
#Script for meta-analysis of changes in aboveground biomass with logging###############
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:16/09/2013

#open packages
library(RODBC)
library(ggplot2)
library(metafor)
library(GGally)
library(multcomp)


###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################

#connect to database
log <- odbcConnect("Logging")
sqlTables(log)

#import data
AGB<- sqlFetch(log, "AGB query")
head(AGB)
colnames(AGB)<-c("Study","Site","Age","Method","Vol","Type","MU","VU","SSU","ML","VL","SSL","Vtype","Scale","Allom","Region")

#recalculate SDs
#unlogged
AGB$SDU<-ifelse(AGB$Vtype=="SE",AGB$VU*sqrt(AGB$SSU),AGB$VU)
AGB$SDU<-ifelse(AGB$Vtype=="CI",(AGB$VU/1.96)*sqrt(AGB$SSU),AGB$SDU)
#logged
AGB$SDL<-ifelse(AGB$Vtype=="SE",AGB$VL*sqrt(AGB$SSL),AGB$VL)
AGB$SDL<-ifelse(AGB$Vtype=="CI",(AGB$VL/1.96)*sqrt(AGB$SSL),AGB$SDL)
head(AGB)

#remove pseudoreplication of Berry study
AGB<-subset(AGB,Age!=18)

#impute missing standard deviation
AGB2<-subset(AGB,AGB$SDU>0)
head(AGB2)
Imp_U<-(sum(AGB2$SDU))/(sum(AGB2$MU))
Imp_L<-(sum(AGB2$SDL))/(sum(AGB2$ML))


#subset data to remove those without vol
AGB_vol<-subset(AGB,Vol>0)

###########################################################################################
####Carry out meta-analysis################################################################
###########################################################################################

#calculate effect sizes
#log ratio
ROM<-escalc(data=AGB,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)

#runs a random effects meta-analysis for log ratio data
ROM.ma<-rma.uni(yi,vi,method="REML",data=ROM)
summary(ROM.ma)

#forrest plot of this
theme_set(theme_bw(base_size=10))
forrest_data<-rbind(data.frame(ES=ROM.ma$yi,SE=sqrt(ROM.ma$vi),Type="Site",Study=AGB$Site),data.frame(ES=ROM.ma$b,SE=ROM.ma$se,Type="Summary",Study="Summary"))
forrest_data$Study2<-factor(forrest_data$Study, levels=rev(levels(forrest_data$Study)) )
levels(forrest_data$Study2)
plot1<-ggplot(data=forrest_data,aes(x=Study2,y=exp(ES)-1,ymax=exp(ES+(1.96*SE))-1,ymin=exp(ES-(1.96*SE))-1,size=(1/SE)/4,colour=factor(Type)))+geom_pointrange(shape=15)
plot2<-plot1+coord_flip()+geom_hline(aes(x=0), lty=2,size=1)
plot3<-plot2+xlab("Study")+ylab("Proportional change")+scale_colour_manual(values=c("grey","black"))
plot3+theme(legend.position="none")+scale_size_continuous(range=c(0.5,2))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Forrest_BM.pdf",height=4,width=6,dpi=1200)

##################################################################################
##analysis of how effects on biomass vary by region, method and age###############
##################################################################################

#compare models with region and method as covariates
model0<-rma(yi,vi,mods=~1,method="ML",data=ROM)
model1<-rma(yi,vi,mods=~Region-1,method="ML",data=ROM)
model2<-rma.uni(yi,vi,mods=~Method-1,method="ML",measure="ROM",data=ROM)
model3<-rma.uni(yi,vi,mods=~Age,method="ML",measure="ROM",data=ROM)

Model_AIC<-data.frame(AIC=c(AIC(model0),AIC(model1),AIC(model2),AIC(model3)))
Model_AIC$model<-c("Null","Model1","Model2","Model3")
#calculate AIC delta
Model_AIC$delta<-Model_AIC$AIC-min(Model_AIC$AIC)
#calculate pseudo r squared for each model
Model_AIC$R_squared<-c(1-(model0$tau2/model0$tau2),1-(model1$tau2/model0$tau2),1-(model2$tau2/model0$tau2),1-(model3$tau2/model0$tau2))
#drop last models with delta >7
AIC_sel<-subset(Model_AIC,delta<=7)
#calculate the relative likelihood of model
AIC_sel$rel_lik<-exp((AIC_sel$AIC[1]-AIC_sel$AIC)/2)
#calculate the AICc weight
AIC_sel$weight<-AIC_sel$rel_lik/(sum(AIC_sel$rel_lik))
#reorder sorting
AIC_sel<-AIC_sel[order(AIC_sel$AIC),]
#put in model parameters
AIC_sel$Vars<-c("Region","Logging method","Null","Age")

#it looks like the one including just region is best
#so we recalculate the parameter estimates using REML
ROM.ma1<-rma(yi,vi,mods=~Region-1,method="REML",data=ROM)
summary(ROM.ma1)

#now we can plot this
Region<-data.frame(coef(summary(ROM.ma1)))
Region$Region<-as.factor(c("Africa","C & S America", "SE Asia & Australasia"))
theme_set(theme_bw(base_size=10))
a<-ggplot(data=Region,aes(y=estimate,ymin=ci.lb,ymax=ci.ub,x=Region))+geom_pointrange(size=1)
b<-a+geom_hline(x=0,lty=2,size=1)+ylab("Proportional change after logging")+xlab("Regions")
b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Region.pdf",height=4,width=6,dpi=1200)

#output of how logging intensity varies by region
a<-ggplot(data=AGB_vol,aes(x=Region,y=Vol))+geom_boxplot()
b<-a+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1,colour="black",fill=NA))
b+ylab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
ggsave("Region_vol.pdf",height=4,width=6,dpi=1200)



#################################################################################
#analysis of the effect of volume on post logging biomass########################
#################################################################################

#log ratio effect size calculation for results with volume
ROM2<-escalc(data=AGB_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
head(ROM2)

#different models relating volume and method to post logging change
Model0<-rma.uni(yi,vi,mods=~1,method="ML",data=ROM2)
Model1<-rma.uni(yi,vi,mods=~I(Vol/MU)+Method,method="ML",data=ROM2)
Model2<-rma.uni(yi,vi,mods=~~Vol+Method,method="ML",data=ROM2)
Model3<-rma.uni(yi,vi,mods=~Vol,method="ML",data=ROM2)
Model4<-rma.uni(yi,vi,mods=~I(Vol/MU),method="ML",data=ROM2)
Model5<-rma.uni(yi,vi,mods=~Vol+MU,method="ML",data=ROM2)
Model6<-rma.uni(yi,vi,mods=~MU,method="ML",data=ROM2)
Model7<-rma.uni(yi,vi,mods=~MU+Method,method="ML",data=ROM2)

AIC(Model0)
AIC(Model1)
AIC(Model2)
AIC(Model3)
AIC(Model4)
AIC(Model5)
AIC(Model6)
AIC(Model7)

Model_AIC<-data.frame(AIC=c(AIC(Model0),AIC(Model1),AIC(Model2),AIC(Model3),AIC(Model4),AIC(Model5),AIC(Model6),AIC(Model7)))
Model_AIC$model<-c("Null","Model1","Model2","Model3","Model4","Model5","Model6","Model7")
#calculate AICc delta
Model_AIC$delta<-Model_AIC$AIC-min(Model_AIC$AIC)
#calculate pseudo r squared for each model
Model_AIC$R_squared<-c(1-(Model0$tau2/Model0$tau2),1-(Model1$tau2/Model0$tau2),1-(Model2$tau2/Model0$tau2),1-(Model3$tau2/Model0$tau2),1-(Model4$tau2/Model0$tau2),1-(Model5$tau2/Model0$tau2),1-(Model6$tau2/Model0$tau2),1-(Model7$tau2/Model0$tau2))
#drop last models with delta >7
AIC_sel<-subset(Model_AIC,delta<=7)
#calculate the realtive likelihood of model
AIC_sel$rel_lik<-exp((AIC_sel$AIC[1]-AIC_sel$AIC)/2)
#calculate the AICc weight
AIC_sel$weight<-AIC_sel$rel_lik/(sum(AIC_sel$rel_lik))
#reorder sorting
AIC_sel<-AIC_sel[order(AIC_sel$AIC),]

#dummy coding for variable importance
AIC_sel$Vol<-c(1,1,1)
AIC_sel$UM<-c(0,1,0)
AIC_sel$Method<-c(0,0,1)

sum(AIC_sel$Vol*AIC_sel$weight)
sum(AIC_sel$UM*AIC_sel$weight)
sum(AIC_sel$Method*AIC_sel$weight)

          
#model 3 is the most parsimonious      
summary(Model3)

#calculate pseudo-rsquared
ROM.ma2<-rma.uni(yi,vi,method="ML",data=ROM2)
summary(ROM.ma2)
1-(deviance(Model3)/deviance(ROM.ma2))

#reset model as REML to get unbiased parameter estimates
Model5<-rma.uni(yi,vi,mods=~Vol,method="REML",data=ROM2)
summary(Model5)

#create dataframe for predictions
head(ROM2)
range(ROM2$Vol)
preds<-predict.rma(Model5,newmods=seq(8.11,179,.1))
Vol<-seq(8.11,179,.1)
preds2<-data.frame(preds=preds$pred,se=preds$se,Vol=Vol,lower=preds$ci.lb,upper=preds$ci.ub)

all<-data.frame(yi=ROM2$yi,vi=ROM2$vi,Vol=ROM2$Vol,Method=ROM2$Method)

#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
AGB_ax<-expression(paste("Unlogged biomass (Mg",ha^-1,")"))
theme_set(theme_bw(base_size=10))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2+theme(legend.position="none")
vol_plot3
vol_plot4<-vol_plot3+ylab("Proportional change in biomass following logging")+geom_point(shape=16,aes(x=Vol,y=exp(yi)-1,colour=Method,size=1/vi))
vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=preds2,aes(x=Vol,y=exp(preds)-1),size=1.5)+geom_hline(y=0,lty=2,size=1)
vol_plot5
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_volume.pdf",height=4,width=6,dpi=1200)

#plot of volume logged vs unlogged biomass

exp(-0.1)-1

theme_set(theme_bw(base_size=10))
AGB_Vol<-ggplot(AGB_vol,aes(x=MU,y=Vol,colour=Region))+geom_point(size=4)
AGB_Vol2<-AGB_Vol+ylab((expression(paste("Volume of wood logged (",m^3,ha^-1,")"))))
AGB_vol3<-AGB_Vol2+xlab("Mean unlogged biomass in reference plots")
AGB_vol3+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Volume_AGB.pdf",height=4,width=6,dpi=1200)