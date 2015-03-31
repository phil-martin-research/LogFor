#######################################################################################
#Script for meta-analysis of changes in aboveground biomass with logging###############
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:31/08/2014

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(metafor)
library(GGally)
library(multcomp)
library(MuMIn)
library(plyr)


###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################


setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
AGB<-read.csv("AGB_intens.csv")

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



###########################################################################################
####Carry out meta-analysis for reviewers to look at plot size and variability#############
###########################################################################################

#calculate effect sizes
#log ratio
ROM<-escalc(data=AGB,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
ROM

setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
write.csv(ROM,"AGB_studies.csv")

#plot how grain size (plot size) affects study between sample variability
theme_set(theme_bw(base_size=12))
a<-ggplot(ROM,aes(x=Plot_size,y=1/vi))+geom_point(size=4,shape=1)+scale_x_log10()+scale_y_log10()
b<-a+xlab("grain size (ha)")+ylab("relative weight contributed\nto meta-analysis")
b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("AGB_grain.jpeg",height=4,width=6,dpi=600)

#plot how extent of studies affects study between sample variability
theme_set(theme_bw(base_size=12))
a<-ggplot(ROM,aes(x=Plot_size*(SSU+SSL),y=1/vi))+geom_point(shape=1,size=4)+scale_x_log10()
b<-a+xlab("study extent (ha)")+ylab("relative weight contributed\nto meta-analysis")
b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
ggsave("AGB_extent.jpeg",height=4,width=6,dpi=600)


#########################################################################################
#imputation of standard deviations#######################################################
#########################################################################################

#now calculate cv based on plot size
Rom_plot<-subset(ROM,Plot_size>0&SDU>0)

#first for unlogged plots
Plot_MU<-lm(log(I(SDU/MU))~(Plot_size),data=Rom_plot)
summary(Plot_MU)
plot(Plot_MU)
plot((Rom_plot$Plot_size),-1/sqrt(Rom_plot$SDU/Rom_plot$MU),)
plot((Rom_plot$Plot_size),(Rom_plot$SDL/Rom_plot$MU),)

#now for logged plots
Plot_ML<-lm(I(SDL/ML)~(Plot_size)+I(Plot_size^2),data=Rom_plot)
summary(Plot_ML)
plot(Plot_ML)

Preds1<-(data.frame(PS=(Rom_plot$Plot_size),Pred=predict(Plot_ML)))
head(Preds1)
Preds1<-Preds1[with(Preds1, order(PS)), ]
plot((Rom_plot$Plot_size),(Rom_plot$SDL/Rom_plot$ML),)
lines(Preds1)

#now create loop to fill in gaps
head(AGB)

for (i in 1:nrow(AGB)){
  AGB$SDU[i]<-ifelse(is.na(AGB$SDU[i])&AGB$Plot_size[i]>1,AGB$MU*0.15,AGB$SDU[i])
  AGB$SDU[i]<-ifelse(is.na(AGB$SDU[i])&AGB$Plot_size[i]<1,AGB$MU*0.25,AGB$SDU[i])
  AGB$SDL[i]<-ifelse(is.na(AGB$SDL[i])&AGB$Plot_size[i]>1,AGB$ML*0.2,AGB$SDU[i])
  AGB$SDL[i]<-ifelse(is.na(AGB$SDL[i])&AGB$Plot_size[i]<1,AGB$ML*0.5,AGB$SDU[i])
}


#################################################################################
#analysis of the effect of volume on post logging biomass########################
#################################################################################

#subset data to remove those without vol
AGB_vol<-subset(AGB,Vol2>0)

#log ratio effect size calculation for results with volume
ROM2<-escalc(data=AGB_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
ROM2<-subset(ROM2,Vol>0)

setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.csv(ROM2,"AGB_studies_vol.csv")


#replace conventional with 1 and RIl with 0
ROM2$Method2<-as.factor(ifelse(ROM2$Method=="Conventional",1,0))

#different models relating volume and method to post logging change
Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|ID),method="ML",data=ROM2)
Model1<-rma.mv(yi,vi,mods=~Age,random=list(~1|ID),method="ML",data=ROM2)
Model2<-rma.mv(yi,vi,mods=~Vol2,random=list(~1|ID),method="ML",data=ROM2)
Model3<-rma.mv(yi,vi,mods=~Method,random=list(~1|ID),method="ML",data=ROM2)
Model4<-rma.mv(yi,vi,mods=~Vol2*Method,random=list(~1|ID),method="ML",data=ROM2)
Model5<-rma.mv(yi,vi,mods=~Vol2*Age+Vol2*Method,random=list(~1|ID),method="ML",data=ROM2)
Model6<-rma.mv(yi,vi,mods=~Vol2+I(Vol2^2),random=list(~1|ID),method="ML",data=ROM2)
Model7<-rma.mv(yi,vi,mods=~Vol2*Method+I(Vol2^2)*Method,random=list(~1|ID),method="ML",data=ROM2)
Model8<-rma.mv(yi,vi,mods=~Vol2*Age,random=list(~1|ID),method="ML",data=ROM2)
Model9<-rma.mv(yi,vi,mods=~Vol2*Method+I(Vol2^2),random=list(~1|ID),method="ML",data=ROM2)
Model10<-rma.mv(yi,vi,mods=~Vol2*Region,random=list(~1|ID),method="ML",data=ROM2)

(Model0$sigma-Model2$sigma)/Model0$sigma

1-(Model4$sigma/Model0$sigma)

nrow(ROM)
summary(ROM$Age)


Model_AICc<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],Model3$fit.stats$ML[5],Model4$fit.stats$ML[5],Model5$fit.stats$ML[5],Model6$fit.stats$ML[5],Model7$fit.stats$ML[5],Model8$fit.stats$ML[5],Model9$fit.stats$ML[5],Model10$fit.stats$ML[5]))
Model_AICc$model<-c("Null","Model1","Model2","Model3","Model4","Model5","Model6","Model7","Model8","Model9","Model10")

#calculate AICc delta
Model_AICc$delta<-Model_AICc$AICc-min(Model_AICc$AICc)
#calculate pseudo r squared for each model
Null_sigma<-(sum(Model0$sigma2))
 
Model_AICc$sigma<-c(sum(Model0$sigma2),sum(Model1$sigma2),
              sum(Model2$sigma2),sum(Model3$sigma2),
              sum(Model4$sigma2),sum(Model5$sigma2),
              sum(Model6$sigma2),sum(Model7$sigma2),
              sum(Model8$sigma2),sum(Model9$sigma2),
              sum(Model10$sigma2))

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
AICc_sel$vars<-c("Vol*Method","Vol","Vol+Vol^2","Vol+Vol^2+Vol^3","Vol*Age","Vol*Age+Vol*Method","Vol+Vol^2+Vol^3+Method","Age","Null model","Method")

setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.csv(AICc_sel,file="AGB_vol.csv")

#re-do model with REML
Model8_reml<-rma.mv(yi,vi,mods=~Vol2*Method,random=list(~1|ID),method="REML",data=ROM2)

plot(predict(Model8_reml)$pred,resid(Model8_reml))

#create dataframe for predictions
all<-data.frame(yi=ROM2$yi,vi=ROM2$vi,Vol=ROM2$Vol2,Method=ROM2$Method,Age=ROM2$Age,MU=ROM2$MU,Study=ROM2$Study.x,Trees=ROM2$Trees2)
all$preds<-(predict(Model8_reml,level=0))$pred
all$ci.lb<-(predict(Model8_reml,level=0))$ci.lb
all$ci.ub<-(predict(Model8_reml,level=0))$ci.ub
Vol<-seq(0,179,length.out=500)



Method<-rep(c("Conventional","RIL"),times = 250)
new_df<-data.frame(Vol2=c(0,50,100,150,0,50,100,150),
           Method2=as.factor(c(1,1,1,1,0,0,0,0)))


preds<-predict(Model8_reml,addx=T)
preds[,7]
preds2<-data.frame(print(preds))                                           

                                           

new_preds<-data.frame(preds=as.numeric(preds2$pred),ci.lb=as.numeric(preds2$ci.lb),
                      ci.ub=as.numeric(preds2$ci.ub),Vol=as.numeric(preds2$X.Vol2),Method=all$Method)


#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
theme_set(theme_bw(base_size=25))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2
vol_plot4<-vol_plot3+ylab("Proportional change in biomass following logging")
vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=new_preds,aes(x=Vol,y=exp(preds)-1,group=Method,colour=Method),size=3)+geom_hline(y=0,lty=2,size=2)
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot7<-vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
biommass_vol_plot<-vol_plot7+geom_line(data=new_preds,aes(y=exp(ci.lb)-1,x=Vol,group=Method,colour=Method),lty=3,size=2)+geom_line(data=new_preds,aes(y=exp(ci.ub)-1,x=Vol,group=Method,colour=Method),lty=3,size=2)
biommass_vol_plot+geom_point(data=all,shape=16,aes(x=Vol,y=exp(yi)-1,colour=Method,size=1/vi),alpha=0.5)+scale_size_continuous(range=c(8,16))+ guides(colour = guide_legend(override.aes = list(size=18)))+ theme(legend.position="none")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_volume2.png",height=12,width=12,dpi=400)



########################################################################
#now run model subsetting the data to removed pinard study##############
########################################################################

head(ROM2)

nrow(ROM2)

ROM3<-subset(ROM2,Study.x!="Pinard et al 1996")

#different models relating volume and method to post logging change
Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|ID),method="ML",data=ROM3)
Model1<-rma.mv(yi,vi,mods=~~Age,random=list(~1|ID),method="ML",data=ROM3)
Model2<-rma.mv(yi,vi,mods=~Vol2,random=list(~1|ID),method="ML",data=ROM3)
Model3<-rma.mv(yi,vi,mods=~Method,random=list(~1|ID),method="ML",data=ROM3)
Model4<-rma.mv(yi,vi,mods=~Vol2*Method,random=list(~1|ID),method="ML",data=ROM3)
Model5<-rma.mv(yi,vi,mods=~Vol2*Age+Vol2*Method,random=list(~1|ID),method="ML",data=ROM3)
Model6<-rma.mv(yi,vi,mods=~Vol2+I(Vol2^2),random=list(~1|ID),method="ML",data=ROM3)
Model7<-rma.mv(yi,vi,mods=~Vol2*Method+I(Vol2^2)*Method,random=list(~1|ID),method="ML",data=ROM3)
Model8<-rma.mv(yi,vi,mods=~Vol2*Age-1,random=list(~1|ID),method="ML",data=ROM3)
Model9<-rma.mv(yi,vi,mods=~Vol2*Method+I(Vol2^2),random=list(~1|ID),method="ML",data=ROM3)
Model10<-rma.mv(yi,vi,mods=~Vol2*Region-1,random=list(~1|ID),method="ML",data=ROM3)



Model_AICc<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],Model3$fit.stats$ML[5],Model4$fit.stats$ML[5],Model5$fit.stats$ML[5],Model6$fit.stats$ML[5],Model7$fit.stats$ML[5],Model8$fit.stats$ML[5],Model9$fit.stats$ML[5],Model10$fit.stats$ML[5]))
Model_AICc$model<-c("Null","Model1","Model2","Model3","Model4","Model5","Model6","Model7","Model8","Model9","Model10")

#calculate AICc delta
Model_AICc$delta<-Model_AICc$AICc-min(Model_AICc$AICc)
#calculate pseudo r squared for each model
Null_sigma<-(sum(Model0$sigma2))

Model_AICc$sigma<-c(sum(Model0$sigma2),sum(Model1$sigma2),
                    sum(Model2$sigma2),sum(Model3$sigma2),
                    sum(Model4$sigma2),sum(Model5$sigma2),
                    sum(Model6$sigma2),sum(Model7$sigma2),
                    sum(Model8$sigma2),sum(Model9$sigma2),
                    sum(Model10$sigma2))

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
AICc_sel$vars<-c("Vol*Method","Vol","Vol*Method+Vol^2","Vol+*Age+Vol*Method","Vol*Age","Vol*Vol^2","Vol*Region","Vol*Method=Vol^2*Method","Null model","Age","Method")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables")
write.csv(AICc_sel,file="AGB_vol_no_pinard.csv")

#re-do model with REML
Model8_reml<-rma.mv(yi,vi,mods=~Vol2*Method-1,random=list(~1|ID),method="REML",data=ROM3)

plot(predict(Model8_reml)$pred,resid(Model8_reml))

#create dataframe for predictions
all<-data.frame(yi=ROM3$yi,vi=ROM3$vi,Vol=ROM3$Vol2,Method=ROM3$Method,Age=ROM3$Age,MU=ROM3$MU,Study=ROM3$Study.x,Trees=ROM3$Trees2)
all$preds<-(predict(Model8_reml,level=0))$pred
all$ci.lb<-(predict(Model8_reml,level=0))$ci.lb
all$ci.ub<-(predict(Model8_reml,level=0))$ci.ub

preds<-predict(Model8_reml,addx=T)
preds[,7]
preds2<-data.frame(print(preds))                                           



new_preds<-data.frame(preds=as.numeric(preds2$pred),ci.lb=as.numeric(preds2$ci.lb),
                      ci.ub=as.numeric(preds2$ci.ub),Vol=as.numeric(preds2$X.Vol2),Method=all$Method)



#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
theme_set(theme_bw(base_size=25))
vol_plot<-ggplot(data=all)
vol_plot2<-vol_plot
vol_plot3<-vol_plot2
vol_plot4<-vol_plot3+ylab("Proportional change in biomass following logging")
vol_plot5<-vol_plot4+scale_size_continuous(range=c(5,10))+geom_line(data=new_preds,aes(x=Vol,y=exp(preds)-1,group=Method,colour=Method),size=3)+geom_hline(y=0,lty=2,size=2)
vol_plot6<-vol_plot5+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot7<-vol_plot6+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
biommass_vol_plot<-vol_plot7+geom_line(data=new_preds,aes(y=exp(ci.lb)-1,x=Vol,group=Method,colour=Method),lty=3,size=2)+geom_line(data=new_preds,aes(y=exp(ci.ub)-1,x=Vol,group=Method,colour=Method),lty=3,size=2)
biommass_vol_plot+geom_point(data=all,shape=16,aes(x=Vol,y=exp(yi)-1,colour=Method,size=1/vi),alpha=0.5)+scale_size_continuous(range=c(8,16))+ guides(colour = guide_legend(override.aes = list(size=18)))+ theme(legend.position="none")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_volume2_no_pinard.jpeg",height=12,width=12,dpi=1200)

#plot to whos the difference in logging intesnity for RIl and conventional
ggplot(data=ROM2,aes(x=Vol2,fill=Method))+geom_density()+facet_wrap(~Method,scales="free_x")

