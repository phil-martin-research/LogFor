#######################################################################################
#Script for meta-analysis of changes in aboveground biomass with logging###############
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:29/06/2015

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(metafor)
library(GGally)
library(MuMIn)
library(plyr)
#and functions
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################
AGB<-read.csv("Data/AGB.csv")

#recalculate Standard deviations depending on the error type used in the paper
#unlogged sites
AGB$SDU<-ifelse(AGB$VarT=="SE",AGB$VU*sqrt(AGB$SSU),AGB$VU)
AGB$SDU<-ifelse(AGB$VarT=="CI",(AGB$VU/1.96)*sqrt(AGB$SSU),AGB$SDU)
#logged sites
AGB$SDL<-ifelse(AGB$VarT=="SE",AGB$VL*sqrt(AGB$SSL),AGB$VL)
AGB$SDL<-ifelse(AGB$VarT=="CI",(AGB$VL/1.96)*sqrt(AGB$SSL),AGB$SDL)
#calculate log ratio effect size
ROM<-escalc(data=AGB,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)


#########################################################################################
#imputation of standard deviations#######################################################
#########################################################################################

#now calculate coeffient of variation based on plot size
Rom_plot<-subset(ROM,Plot_size>0&SDU>0&!is.na(vi))
Rom_plot$U_CV<-Rom_plot$SDU/Rom_plot$MU
Rom_plot$L_CV<-Rom_plot$SDL/Rom_plot$ML

#now model coefficient of varion for unlogged plots so it can be estimated where SE is not reported
Plot_UL1<-lm(U_CV~Plot_size+I(Plot_size^2),data=Rom_plot)
Plot_UL2<-lm(U_CV~log(Plot_size),data=Rom_plot)
#now model coefficient of varion for logged plots
Plot_ML2<-lm(L_CV~log(Plot_size),data=Rom_plot)
#now use a loop to predict coeffient of variation where no SE was reported
for (i in 1:nrow(AGB)){
  AGB$SDL[i]<-ifelse(AGB$Plot_size[i]<=4&is.na(AGB$SDU[i]),(0.13333+(-0.17945*(log(AGB$Plot_size[i]))))*AGB$ML,AGB$SDU[i])
  AGB$SDU[i]<-ifelse(AGB$Plot_size[i]<=4&is.na(AGB$SDU[i]),(0.13870+(-0.12938*(log(AGB$Plot_size[i]))))*AGB$MUL,AGB$SDU[i]) 
  AGB$SDU[i]<-ifelse(AGB$Plot_size[i]>=4&is.na(AGB$SDU[i]),0.15*AGB$MUL,AGB$SDU[i])
  AGB$SDL[i]<-ifelse(AGB$Plot_size[i]>=4&is.na(AGB$SDU[i]),0.1*AGB$ML,AGB$SDL[i]) 
}

#################################################################################
#analysis of the effect of volume on post logging biomass########################
#################################################################################
#subset the dataset to remove studies without measure of logging intensity
AGB_vol<-subset(AGB,Vol2>0)
#log ratio effect size calculation for results with volume
ROM2<-escalc(data=AGB_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
ROM2<-subset(ROM2,!is.na(vi))
write.csv(ROM2,"Data/AGB_studies_vol.csv")
#replace conventional with 1 and RIL with 0
ROM2$Method2<-as.factor(ifelse(ROM2$Method=="Conventional",1,0))

#boostrap different models relating intensity, method and age to post logging change in biomass
Site_unique<-unique(ROM2$MU)
Model_AIC_summary<-NULL
for (i in 1:10000){ #boostrap 10,000 times
  print(i)
  AGB_samp<-NULL
  for (j in 1:length(Site_unique)){#sample one site for each study so that no reference site is used more than once
    AGB_sub<-subset(ROM2,MU==Site_unique[j])
    AGB_sub<-AGB_sub[sample(nrow(AGB_sub), 1), ] 
    AGB_samp<-rbind(AGB_sub,AGB_samp)
  } #run different models with the subsampled data
  Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model1<-rma.mv(yi,vi,mods=~Age,random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model2<-rma.mv(yi,vi,mods=~Vol2,random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model3<-rma.mv(yi,vi,mods=~Method,random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model4<-rma.mv(yi,vi,mods=~Vol2*Method,random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model5<-rma.mv(yi,vi,mods=~Vol2*Age+Vol2*Method,random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model6<-rma.mv(yi,vi,mods=~Vol2+I(Vol2^2),random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model7<-rma.mv(yi,vi,mods=~Vol2*Method+I(Vol2^2)*Method,random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model8<-rma.mv(yi,vi,mods=~Vol2*Age,random=list(~1|Study,~1|Age),method="ML",data=AGB_samp)
  Model9<-rma.mv(yi,vi,mods=~Vol2*Method+I(Vol2^2),random=list(~1|Study),method="ML",data=AGB_samp)
  Model_AIC<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],#produce AICc values for the models
                               Model3$fit.stats$ML[5],Model4$fit.stats$ML[5],Model5$fit.stats$ML[5],
                               Model6$fit.stats$ML[5],Model7$fit.stats$ML[5],Model8$fit.stats$ML[5],
                               Model9$fit.stats$ML[5]))
  Model_AIC$Vars<-c("Null","Age","Volume","Method","Volume*Method","Volume*Age+Volume*Method", #details of model variables
                    "Voume+Volume^2","Volume*Method+Volume^2*Method","Volume*Age","Volume*Method+Vol^2")
  Model_AIC$logLik<-c(Model0$fit.stats$ML[1],Model1$fit.stats$ML[1],Model2$fit.stats$ML[1],#put logLiklihood in the table
                       Model3$fit.stats$ML[1],Model4$fit.stats$ML[1],Model5$fit.stats$ML[1],
                       Model6$fit.stats$ML[1],Model7$fit.stats$ML[1],Model8$fit.stats$ML[1],
                       Model9$fit.stats$ML[1])
  Null_dev<-deviance(Model0)
  Dev<-c(deviance(Model0),deviance(Model1),deviance(Model2),deviance(Model3),deviance(Model4),deviance(Model5),#calculate deviance of models
           deviance(Model6),deviance(Model7),deviance(Model8),deviance(Model9))
  Model_AIC$R2<-1-(Dev/Null_dev) #calculate pseudo-r squared using model deviance
  Model_AIC$R2<-ifelse(Model_AIC$R2<0,0,Model_AIC$R2)
  Model_AIC<-Model_AIC[order(Model_AIC$AICc),] #reorder from lowest to highest
  Model_AIC$delta<-Model_AIC$AICc-Model_AIC$AICc[1]#calculate AICc delta
  Model_AIC$rel_lik<-exp((Model_AIC$AICc[1]-Model_AIC$AICc)/2)#calculate the relative likelihood of model
  Model_AIC$weight<-Model_AIC$rel_lik/(sum(Model_AIC$rel_lik))
  Model_AIC$Run<-i
  Model_AIC$Rank<-seq(1,10,1) #rank models from 1-10 in terms of parsimony
  Model_AIC_summary<-rbind(Model_AIC,Model_AIC_summary)
}
Model_AIC_summary$Rank1<-ifelse(Model_AIC_summary$Rank==1,1,0)

#summarise the boostrapping routine by giving median values for model statistics - log liklihood, AICc delta AICc, R squared
Model_sel_boot<-ddply(Model_AIC_summary,.(Vars),summarise,Modal_rank=Mode(Rank),Prop_rank=sum(Rank1)/10000,log_liklihood=median(logLik),AICc_med=median(AICc),
                      delta_med=median(delta),R2_med=median(R2))
#write this summary table to file
write.csv(Model_sel_boot,file="Tables/AGB_mod_sel.csv")

#now boostrap the top model to get parameter estimates
Site_unique<-unique(ROM2$MU)
Param_boot<-NULL
for (i in 1:10000){
  print(i)
  AGB_samp<-NULL
  for (j in 1:length(Site_unique)){#use same routine as previously to subsample dataset avoiding pseudo-replication
    AGB_sub<-subset(ROM2,MU==Site_unique[j])
    AGB_sub<-AGB_sub[sample(nrow(AGB_sub), 1), ]
    AGB_samp<-rbind(AGB_sub,AGB_samp)
  }
  Model1_Vol<-rma.mv(yi,vi,mods=~Vol2,random=list(~1|Study,~1|Age),method="REML",data=AGB_samp)
  Param_vals<-data.frame(Parameter=c("Intercept","Vol_slope"),estimate=coef(summary(Model1_Vol))[1],se=coef(summary(Model1_Vol))[2],
                         pval=coef(summary(Model1_Vol))[4],ci_lb=coef(summary(Model1_Vol))[5],ci_ub=coef(summary(Model1_Vol))[6])
  Param_boot<-rbind(Param_vals,Param_boot)
}

#produce summary of parameter estimates
Param_boot_sum<-ddply(Param_boot,.(Parameter),summarise,coef_estimate=median(estimate),lower=median(ci.lb),
                      upper=median(ci.ub),med_pval=median(pval),se=median(se))
#write this table of parameter estimates
write.table(Param_boot_sum,file="Tables/Biomass_parameter_estimates.csv",sep=",")

#create dataframe for predictions
newdat<-data.frame(Vol=seq(5.7,179,length.out=500))
newdat$yi<-Param_boot_sum$coef_estimate[1]+(newdat$Vol*Param_boot_sum$coef_estimate[2])
newdat$UCI<-(Param_boot_sum$upper[1])+(Param_boot_sum$upper[2]*newdat$Vol)
newdat$LCI<-(Param_boot_sum$lower[1])+(Param_boot_sum$lower[2]*newdat$Vol)
#wirte predictions to file
write.csv(newdat,"Data/Preds_AGB.csv",row.names=F)
