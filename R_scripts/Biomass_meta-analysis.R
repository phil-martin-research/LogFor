#######################################################################################
#Script for meta-analysis of changes in aboveground biomass with logging###############
#and plots for paper###################################################################
#######################################################################################

#test

#name: Phil Martin
#date:22/07/2013

#open packages
library(RODBC)
library(ggplot2)
library(metafor)

###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################

#connect to database
log <- odbcConnect("Logging")
sqlTables(log)

#import data
AGB<- sqlFetch(log, "AGB query")
head(AGB)
colnames(AGB)<-c("Study","Site","Age","Method","Vol","Type","MU","VU","SSU","ML","VL","SSL","Vtype","Scale","Allom")

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

plot(ROM.ma)

exp(-.4998)-1

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

#cumulative meta-analysis-organise data
order1<-order(sort(ROM.ma$vi))
cum_meta<-cumul(ROM.ma,order=order(sort(ROM.ma$vi)))
cum_data<-data.frame(estimate=as.numeric(cum_meta$estimate),ci.ub=cum_meta$ci.ub,ci.lb=cum_meta$ci.lb,order=as.numeric(seq(1,23)))

#plot cumulative meta-analysis
theme_set(theme_bw(base_size=26))
cum_plot<-ggplot(cum_data,aes(x=order,y=exp(estimate)-1,ymax=exp(ci.ub)-1,ymin=exp(ci.lb)-1))+geom_ribbon(alpha=0.5)
cum_plot+ylab("Proportional change in biomass following logging")+xlab("Number of studies")+geom_line(aes(x=order,y=exp(estimate)-1))

#look at method and allometric method as covariates
ROM.ma1<-rma.uni(yi,vi,mods=~Method+Allom,method="ML",data=ROM)
ROM.ma2<-rma.uni(yi,vi,mods=~Method,method="ML",data=ROM)


#test which model is the most parsimonious
AIC(ROM.ma1)
AIC(ROM.ma2)

#looks like it's model 2, let's recalculate it using REML
#so there are non-biased estimates of the covariates
ROM.ma2<-rma.uni(yi,vi,mods=~Method,method="REML",measure="ROM",data=ROM)
summary(ROM.ma2)
exp(-0.5763)-1
exp(-0.5763+.3003)-1

#put into a dataframe
Methods<-data.frame(coef(summary(ROM.ma2)))
Methods$methods<-c("Conventional","RIL")
Methods[2,1]<-Methods[2,1]+Methods[1,1]
#plot results
theme_set(theme_bw(base_size=10))
a<-ggplot(Methods,aes(x=methods,y=exp(estimate)-1,ymin=exp(estimate-(se*1.96))-1,ymax=exp(estimate+(se*1.96))-1))+geom_pointrange(size=1.5)
b<-a+coord_flip()+geom_hline(x=0,lty=2,size=2)+ylab("Proportional change after logging")+xlab("Methods")
b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("RIL_Conv.jpeg",height=4,width=6,dpi=1200)

#log ratio of results with volume
ROM2<-escalc(data=AGB_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
head(ROM2)

#look at proportional volume change
Model1<-rma.uni(yi,vi,mods=~I(Vol/MU)*Method,method="ML",data=ROM2)
Model1
Model2<-rma.uni(yi,vi,mods=~I(Vol/MU)+Method,method="ML",data=ROM2)
summary(Model2)
Model3<-rma.uni(yi,vi,mods=~I(Vol/MU),method="ML",data=ROM2)
Model3

plot(ROM2$MU,ROM2$yi)
AIC(Model1)
AIC(Model2)
AIC(Model3)
          
#it looks like model 3 is the best          
summary(Model3)

#calculate pseudo-rsquared
ROM.ma2<-rma.uni(yi,vi,method="ML",data=ROM2)
summary(ROM.ma2)
1-(49.11761/56)

#reset model as REML
Model5<-rma.uni(yi,vi,mods=~I(Vol/MU),method="REML",data=ROM2)
summary(Model5)

#create dataframe for predictions
head(ROM2)
range(ROM2$Vol/ROM2$MU)
preds<-predict.rma(Model5,newmods=seq(0.019,0.47,0.001))
preds2<-data.frame(preds=preds$pred,se=preds$se,Vol=seq(0.019,0.47,0.001),lower=preds$ci.lb,upper=preds$ci.ub)

all<-data.frame(yi=ROM2$yi,vi=ROM2$vi,Vol=ROM2$Vol/ROM2$MU,preds=(predict.rma(Model3)$pred),upper=(predict.rma(Model3)$ci.ub),lower=(predict.rma(Model3)$ci.lb),Method=ROM2$Method)

#plot results
#first crate x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
AGB_ax<-expression(paste("Unlogged biomass (Mg",ha^-1,")"))
theme_set(theme_bw(base_size=10))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2+theme(legend.position="none")
vol_plot3
vol_plot4<-vol_plot3+ylab("Change in biomass following logging")+geom_point(shape=16,aes(x=Vol,y=exp(yi)-1,colour=Method,size=1/vi))
vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=preds2,aes(x=Vol,y=exp(preds)-1),size=1.5)+geom_hline(y=0,lty=2,size=2)
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot6+xlab(expression(frac(paste("Volume of wood logged (",m^3,ha^-1,")")
                               ,paste("Unlogged biomass (Mg ",ha^-1,")"))))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_volume.jpeg",height=4,width=6,dpi=1200)

?as.expression
expression(alpha <= beta)
?expression
?paste
