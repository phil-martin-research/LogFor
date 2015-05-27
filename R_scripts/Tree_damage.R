#############################################################
#script to predict proportion of trees damaged using volume##
#############################################################

#first load packages
library(ggplot2)
library(MuMIn)
library(lme4)

#clear objects
rm(list=ls())

#import data
Dam<-read.csv("Data/prop_damage.csv")
head(Dam)
colnames(Dam)<-c("Row","ID","Study","Site","Method","Vol","Tree_Ex","Prop_dam","Region","Plot","Replicates")
Dam<-Dam[complete.cases(Dam[,c(6,8)]),]

#######################################################
#now predict proportion of trees damaged using volume##
#######################################################

#create column for Volume squared and log volume
Dam$Vol_sq<-Dam$Vol^2
Dam$Vol_log<-log(Dam$Vol)
Dam$Vol_log2<-(log(Dam$Vol))^2
head(Dam)
Dam$Replicates<-ifelse(is.na(Dam$Replicates),1,Dam$Replicates)
Dam$Plot<-ifelse(is.na(Dam$Plot),median(Dam$Plot,na.rm = T),Dam$Plot)


ggplot(Dam,aes(x=Vol,y=Prop_dam,group=Study,size=Plot*Replicates))+geom_point()+geom_smooth(method=lm,aes(group=Method))





#test for the effects of volume logged, non-linear volume logged (squared and log), differences in method
#and a null model
M1<-lmer(qlogis(Prop_dam)~Vol+(Vol|Study),Dam,weight=I(Plot*Replicates),na.action="na.fail")
M2<-lmer(qlogis(Prop_dam)~Vol_log*Method+(Vol|Study),weight=I(Plot*Replicates),Dam,na.action="na.fail")
M3<-lmer(qlogis(Prop_dam)~Method+(Vol|Study),weight=I(Plot*Replicates),Dam,na.action="na.fail")
M4<-lmer(qlogis(Prop_dam)~Vol*Method+(Vol|Study),weight=I(Plot*Replicates),Dam,na.action="na.fail")
M0<-lmer(qlogis(Prop_dam)~1+(Vol|Study),weight=I(Plot*Replicates),Dam,na.action="na.fail")
AICc(M1,M2,M3,M4,M0)
All_mods<-list(M1,M2,M3,M4,M0)

#diagnostic plots

plot(Dam$Vol,plogis(predict(M2,re.form=NA)))

#model averaging
ms1<-model.sel(object=All_mods,rank="AICc",fit=T,trace=T,subset=dc(Vol2,Vol_sq))
ms1$r2<-c(r.squaredGLMM(M2)[1],r.squaredGLMM(M3)[1],r.squaredGLMM(M1)[1],
          r.squaredGLMM(M0)[1],r.squaredGLMM(M4)[1])
#write this table to csv
write.csv(ms1,"Tables/Tree_damage_models.csv")

# extract coefficients
coefs <- data.frame(coef(summary(M2)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
write.csv(coefs,"Tables/Tree_damage_coefs.csv")

#predict values for plotting
#now create plots of this
newdat<-data.frame(Vol=c(seq(3,164.9,length.out = 500),seq(5,107,length.out = 500)),Method=c(rep("Conventional",500),rep("RIL",500)))
newdat$Vol_log<-log(newdat$Vol)
newdat$Prop_dam<-0
mm <- model.matrix(terms(M2),newdat)
newdat$Prop_dam <- predict(M2,newdat,re.form=NA)
## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(M2),mm))
tvar1 <- pvar1+VarCorr(M2)$Study[1]  ## must be adapted for more complex models
tvar1 <- 
  newdat <- data.frame(
    newdat
    , plo = newdat$Prop_dam-2*sqrt(pvar1)
    , phi = newdat$Prop_dam+2*sqrt(pvar1)
    , tlo = newdat$Prop_dam-2*sqrt(tvar1)
    , thi = newdat$Prop_dam+2*sqrt(tvar1)
  )


#plot these results
theme_set(theme_bw(base_size=12))
Dam_1<-ggplot(Dam,aes(x=Vol,y=Prop_dam,colour=Method,size=Replicates*Plot))+geom_point(shape=1)+geom_line(data=newdat,aes(x=Vol,y=plogis(Prop_dam),colour=Method),size=2)
Dam_2<-Dam_1+geom_ribbon(data=newdat,aes(ymax=plogis(phi),ymin=plogis(plo),fill=Method,size=NULL),colour=NA,alpha=0.2)+ylim(0,0.7)
Dam_3<-Dam_2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Dam_3+xlab("Volume logged per hectare")+ylab("Proportion of residual tree stems damaged")+scale_colour_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+theme(legend.position="none")
ggsave("Figures/Prop_damaged_vol.png",height=6,width=8,dpi=1200)


