#######################################################################################
#Script for meta-analysis of changes in species richness with logging##################
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:13/04/2015

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(metafor)
library(MuMIn)
library(boot)
library(plyr)
library(reshape2)
#and functions needed
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


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

Site_unique<-unique(Rich_vol$M_UL)
Model_AIC_summary<-NULL
for (i in 1:1000){
print(i)
Rich_samp<-NULL
for (j in 1:length(Site_unique)){
  Rich_sub<-subset(Rich_vol,M_UL==Site_unique[j])
  Rich_sub<-Rich_sub[sample(nrow(Rich_sub), 1), ]
  Rich_samp<-rbind(Rich_sub,Rich_samp)
}
Model0_Vol<-rma.mv(yi,vi,mods=~1,random=list(~ 1 | Rare),method="ML",data=Rich_samp)
Model1_Vol<-rma.mv(yi,vi,mods=~Vol_std,random=list( ~ 1 | Rare),method="ML",data=Rich_samp)
Model2_Vol<-rma.mv(yi,vi,mods=~Method,random=list( ~ 1 | Rare),method="ML",data=Rich_samp)
Model3_Vol<-rma.mv(yi,vi,mods=~Method*Vol_std,random=list( ~ 1 | Rare),method="ML",data=Rich_samp)
Model4_Vol<-rma.mv(yi,vi,mods=~Age_std,random=list( ~ 1 | Rare),method="ML",data=Rich_samp)
Model5_Vol<-rma.mv(yi,vi,mods=~Method+Vol_std,random=list( ~ 1 | Rare),method="ML",data=Rich_samp)
Model_AIC<-data.frame(AICc=c(Model0_Vol$fit.stats$ML[5],Model1_Vol$fit.stats$ML[5],Model2_Vol$fit.stats$ML[5],Model3_Vol$fit.stats$ML[5],Model4_Vol$fit.stats$ML[5],Model5_Vol$fit.stats$ML[5]))
Model_AIC$Vars<-c("Null","Volume",
                  "Method","Volume*Method","Age","Method+Volume")
Model_AIC$logLik<-c(Model0_Vol$fit.stats$ML[1],Model1_Vol$fit.stats$ML[1],Model2_Vol$fit.stats$ML[1],Model3_Vol$fit.stats$ML[1],Model4_Vol$fit.stats$ML[1],Model5_Vol$fit.stats$ML[1])
Null_log<-Model0_Vol$fit.stats$ML[1]
Model_AIC$R2<-1-(Model_AIC$logLik/Null_log)
Model_AIC<-Model_AIC[order(Model_AIC$AICc),] #reorder from lowest to highest
Model_AIC$delta<-Model_AIC$AICc-Model_AIC$AICc[1]#calculate AICc delta
Model_AIC$rel_lik<-exp((Model_AIC$AICc[1]-Model_AIC$AICc)/2)#calculate the relative likelihood of model
Model_AIC$weight<-Model_AIC$rel_lik/(sum(Model_AIC$rel_lik))
Model_AIC$Run<-i
Model_AIC$Rank<-seq(1,6,1)
Model_AIC_summary<-rbind(Model_AIC,Model_AIC_summary)
}

head(Model_AIC_summary)
Model_AIC_summary$Rank1<-ifelse(Model_AIC_summary$Rank==1,1,0)


Model_sel_boot<-ddply(Model_AIC_summary,.(Vars),summarise,Modal_rank=Mode(Rank),Prop_rank=sum(Rank1)/1000,log_liklihood=median(logLik),AICc_med=median(AICc),
      delta_med=median(delta),R2_med=median(R2))

setwd("Tables")
write.table(Model_sel_boot,file="Rich_vol_model_sel.csv",sep=",")


#re-run top model using REML
#bootstrapping 10,000 times to get estimates

Site_unique<-unique(Rich_vol$M_UL)
Param_boot<-NULL
for (i in 1:10000){
  print(i)
  Rich_samp<-NULL
  for (j in 1:length(Site_unique)){
    Rich_sub<-subset(Rich_vol,M_UL==Site_unique[j])
    Rich_sub<-Rich_sub[sample(nrow(Rich_sub), 1), ]
    Rich_samp<-rbind(Rich_sub,Rich_samp)
  }
  Model1_Vol<-rma.mv(yi,vi,mods=~Vol,random=list( ~ 1 | Rare),method="REML",data=Rich_samp)
  Param_vals<-data.frame(Parameter=c("Intercept","Vol_slope"),estimate=coef(summary(Model1_Vol))[1],se=coef(summary(Model1_Vol))[2],
             pval=coef(summary(Model1_Vol))[4],ci_lb=coef(summary(Model1_Vol))[5],ci_ub=coef(summary(Model1_Vol))[6])
  Param_boot<-rbind(Param_vals,Param_boot)
}

head(Param_boot)

Param_boot2<-subset(Param_boot,Parameter=="Vol_slope")
Param_boot3<-subset(Param_boot,Parameter=="Intercept")

hist(Param_boot2$estimate)
quantile(Param_boot2$estimate,probs=c(0.025,0.5,0.975))
hist(Param_boot3$estimate)
quantile(Param_boot3$estimate,probs=c(0.025,0.5,0.975))


ddply(Param_boot,.(Parameter),summarise,coef_estimate=median(estimate),lower=median(ci.lb),
      upper=median(ci.ub),med_pval=median(pval))

summary(Rich_vol$Vol_std)

#create dataframe for predictions
newdat<-data.frame(Vol2=seq(-1.4270,1.6150,length.out=500))
newdat$Vol<-(newdat$Vol2*sd(Rich_vol$Vol))+mean(Rich_vol$Vol)
newdat$yi<-(0.078061590)+(-0.001484976*newdat$Vol)
newdat$UCI<-(0.1205197213)+(-0.0009727783*newdat$Vol)
newdat$LCI<-(0.033975199)+(-0.001997173*newdat$Vol)

(-0.037209107)+(-0.05790845*-1.427000)
(-0.027583018)+(-0.09668837*-1.427000)
(0.002441329)+(-0.02890387*-1.427000)


head(newdat)
ggplot(newdat,aes(x=Vol,y=exp(yi)-1,ymax=UCI,ymin=LCI))+geom_ribbon()+geom_line()

#plot results
#first create x axis labels
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
theme_set(theme_bw(base_size=25))
vol_plot<-ggplot(newdat,aes(x=Vol,y=exp(yi)-1,ymax=UCI,ymin=LCI))+geom_ribbon(alpha=0.2)+geom_line()
vol_plot2<-vol_plot+geom_point(data=ROM_vol,aes(ymax=NULL,ymin=NULL,colour=Method,size=1/vi),shape=1)
vol_plot3<-vol_plot2+ylab("Proportional change in tree species richness following logging")
vol_plot4<-vol_plot3+scale_size_continuous(range=c(5,10))+geom_hline(y=0,lty=2,size=1)
vol_plot5<-vol_plot4+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot6<-vol_plot5+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
rich_vol_plot<-vol_plot6+theme(legend.position="none")+scale_colour_brewer(palette="Set1")
ggsave("Figures/SR_volume.png",height=12,width=12,dpi=400)


#create funnel plot with residuals

plot(resid(Model1_reml),1/sqrt(Rich_vol$vi))
