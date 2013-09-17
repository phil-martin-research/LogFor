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

exp(coef(ROM.ma))-1

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
ggsave("Forrest_BM.jpeg",height=4,width=6,dpi=1200)


#look at region as a covariate
ROM_reg1<-rma(yi,vi,mods=~Region-1,method="ML",data=ROM)
ROM_reg2<-rma(yi,vi,mods=~1,method="ML",data=ROM)

AIC(ROM_reg1)
AIC(ROM_reg2)


#it looks like the one including just region is best
#so we recalculate the parameter estimates using REML
ROM.ma1<-rma(yi,vi,mods=~Region-1,method="REML",data=ROM)
summary(ROM.ma1)

#now we can plot this
Region<-data.frame(Mean=c(-0.2,-0.2884,-0.6622),SE=c(0.2209,0.1117,0.0921))
Region$Region<-as.factor(c("Africa","C & S America", "SE Asia & Australasia"))
ggplot(data=Region,aes(y=Mean,ymin=Mean-(1.96*SE),ymax=Mean+(1.96*SE),x=Region))+geom_pointrange(size=1)

################################################################################
#analysis of logging methods - Conv vs RIL######################################
################################################################################

#mixed effects meta-analysis of Method
ROM_meth1<-rma.uni(yi,vi,mods=~Method-1,method="ML",measure="ROM",data=ROM)
ROM_meth2<-rma.uni(yi,vi,mods=~1,method="ML",measure="ROM",data=ROM)
summary(ROM_meth1)

#put coefficient estimates into a dataframe
Methods<-data.frame(coef(summary(ROM_meth1)))
Methods$methods<-c("Conventional","RIL")

#plot results
theme_set(theme_bw(base_size=10))
a<-ggplot(Methods,aes(x=methods,y=exp(estimate)-1,ymin=exp(estimate-(se*1.96))-1,ymax=exp(estimate+(se*1.96))-1))+geom_pointrange(size=1.5)
b<-a+coord_flip()+geom_hline(x=0,lty=2,size=2)+ylab("Proportional change after logging")+xlab("Methods")
b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("RIL_Conv.jpeg",height=4,width=6,dpi=1200)

#################################################################################
#analysis of the effect of volume on post logging biomass########################
#################################################################################

#log ratio effect size calculation for results with volume
ROM2<-escalc(data=AGB_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
head(ROM2)

#different models relating volume and method to post logging change
Model1<-rma.uni(yi,vi,mods=~I(Vol/MU)+Method,method="ML",data=ROM2)
Model2<-rma.uni(yi,vi,mods=~~Vol+Method,method="ML",data=ROM2)
Model3<-rma.uni(yi,vi,mods=~Vol,method="ML",data=ROM2)

AIC(Model1)
AIC(Model2)
AIC(Model3)
          
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
preds<-predict.rma(Model5,newmods=seq(8.11,179,0.1))
preds2<-data.frame(preds=preds$pred,se=preds$se,Vol=seq(8.11,179,0.1),lower=preds$ci.lb,upper=preds$ci.ub)

all<-data.frame(yi=ROM2$yi,vi=ROM2$vi,Vol=ROM2$Vol,preds=(predict.rma(Model3)$pred),upper=(predict.rma(Model3)$ci.ub),lower=(predict.rma(Model3)$ci.lb),Method=ROM2$Method)

#plot results
#first crate x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
AGB_ax<-expression(paste("Unlogged biomass (Mg",ha^-1,")"))
theme_set(theme_bw(base_size=10))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2+theme(legend.position="none")
vol_plot3
vol_plot4<-vol_plot3+ylab("Proportional change in biomass following logging")+geom_point(shape=16,aes(x=Vol,y=exp(yi)-1,colour=Method,size=1/vi))
vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=preds2,aes(x=Vol,y=exp(preds)-1),size=1.5)+geom_hline(y=0,lty=2,size=2)
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_volume.jpeg",height=4,width=6,dpi=1200)


