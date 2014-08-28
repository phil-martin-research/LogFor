#script to calculate scaling coeffients between different
#measures of logging intensity and damage
#and model damage to proportion of stems


#first load packages
library(nlme)
library(ggplot2)
library(MuMIn)
library(lme4)


#load data
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Prop_dam<-read.csv("prop_damage.csv")
head(Prop_dam)
colnames(Prop_dam)<-c("Study","Site_ID","Age","Method","BA_log","Prop_BA_log","Vol_log","Tree_ex_ha","Dam_tree","Dam_ha","Sev_per_tree","Severe_ha","BA_dam","Prop_ba_dam","Prop_seve_dam","Prop_dam","ID","All","Region","N_logged","Plot","Notes")


#exploratory plots
#volume vs number of trees logged
ggplot(Prop_dam,aes(x=Vol_log,y=Tree_ex_ha,colour=Region,group=Region))+geom_point()+geom_smooth(se=F,method="lm")
#Proportion damaged vs number of trees logged
ggplot(Prop_dam,aes(x=Tree_ex_ha,y=Prop_dam,colour=Method,group=Method))+geom_point()+geom_smooth(se=F,method="lm")
#Proportion damaged vs volume
ggplot(Prop_dam,aes(x=Vol_log,y=Prop_dam,colour=Method,group=Method))+geom_point()+geom_smooth(se=F,method="lm")
#Proportion severe damaged vs damaged
ggplot(Prop_dam,aes(x=Prop_seve_dam,y=Prop_dam))+geom_point()+geom_smooth(se=F,method="lm")+scale_x_continuous(limits=c(0,1))+geom_abline()



#number damaged vs proportion damaged
ggplot(Prop_dam,aes(x=Dam_tree,y=Prop_dam))+geom_point()+geom_smooth(se=F,method="lm")

#proportion BA damaged vs proportion damaged
ggplot(Prop_dam,aes(x=Prop_ba_dam,y=Prop_dam))+geom_point()+geom_smooth(se=F,method="lm")


#model volume as a function of number of trees logged
#create dataset with only complete cases in explanatory variables
Prop_dam_CC<-Prop_dam[complete.cases(Prop_dam[,c(7,8)]),]
Prop_dam_CC$Tree_corr<-Prop_dam_CC$Tree_ex_ha-mean(Prop_dam_CC$Tree_ex_ha)
str(Prop_dam_CC)
M1<-lme(log(Vol_log)~Tree_corr+Region,random=~1|Study,Prop_dam_CC)


plot(Prop_dam_CC$Tree_ex_ha,exp(predict(M1)))

qqnorm(M1)

Tree_pred1<-data.frame(Site_ID=Prop_dam$Site_ID,Tree_corr=Prop_dam$Tree_ex_ha-mean(Prop_dam_CC$Tree_ex_ha),Region=Prop_dam$Region)
Tree_pred1<-Tree_pred1[complete.cases(Tree_pred1),]

#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7)


Tree_pred1$Pred<-exp(predict(avgm,newdata=Tree_pred1,level=0))


M_Damage<-merge(Prop_dam,Tree_pred1,by="Site_ID",all=T)

#put in predicted volume values where there are none
for (i in 1:nrow(M_Damage)){
  M_Damage$Vol_log2[i]<-ifelse(is.na(M_Damage$Vol_log[i])&M_Damage$Tree_ex_ha[i]<25,M_Damage$Pred[i],M_Damage$Vol_log[i])
}



plot(M_Damage$Vol_log2,M_Damage$Prop_dam)

#now predict number of trees logged using volume
M1<-lme(log(Tree_ex_ha)~Vol_log+Region,random=~1|Study,Prop_dam_CC)
plot(M1)
qqnorm(M1)

Tree_pred1<-data.frame(Site_ID=Prop_dam$Site_ID,Vol_log=Prop_dam$Vol_log,Region=Prop_dam$Region)
Tree_pred1<-Tree_pred1[complete.cases(Tree_pred1),]

#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7)

Tree_pred1$Pred<-exp(predict(avgm,newdata=Tree_pred1,level=0))


M_Damage2<-merge(M_Damage,Tree_pred1,by="Site_ID",all=T)

str(M_Damage2)

#put in tree_extracted values where there are none
for (i in 1:nrow(M_Damage2)){
  M_Damage2$Tree_ex_ha2[i]<-ifelse(is.na(M_Damage2$Tree_ex_ha[i]),M_Damage2$Pred.y[i],M_Damage2$Tree_ex_ha[i])
}

summary(M_Damage2$Tree_ex_ha2)

plot(M_Damage2$Tree_ex_ha2,M_Damage2$Prop_dam)
plot(M_Damage2$Vol_log2,M_Damage2$Prop_dam)

(Prop_dam$Vol_log)

ggplot(M_Damage2,aes(x=Tree_ex_ha2,y=Prop_dam,colour=Method))+geom_point()+geom_smooth(se=F, method = 'nls', formula = 'y~a*x^b')+xlim(0,20)


ggplot(M_Damage2,aes(x=Vol_log2,y=Prop_dam,colour=Method))+geom_point()+geom_smooth(se=F, method = 'nls', formula = 'y~a*x^b')


#now remove columns which are no use
keeps <- c("Study","Site_ID","Vol_log2","Tree_ex_ha2","Method","Prop_dam","Prop_seve_dam")
M_Damage3<-M_Damage2[keeps]
M_Damage3<-subset(M_Damage3,Vol_log2<200)
M_Damage3<-subset(M_Damage3,Study!="Guariguata et al 2009")
M_Damage3<-subset(M_Damage3,Prop_dam<0.6)
M_Damage3<-M_Damage3[complete.cases(M_Damage3[,c(3,4,6)]),]

head(M_Damage3)

summary(M_Damage3)

ggplot(M_Damage3,aes(x=Vol_log2,y=Prop_dam,colour=Method))+geom_point()+geom_smooth(se=F, method = 'nls', formula = 'y~a*x^b')



#now predict proportion of trees damaged using volume
M1<-lmer(qlogis(Prop_dam)~Vol_log2*Method+I(Vol_log2^2)*Method+(1|Study),M_Damage3,na.action="na.fail")
summary(M1)

plot(plogis(fitted(M1)),resid(M1))
qqnorm(M1)


xyplot(qlogis(Prop_dam)~ Vol_log2|Study, data=M_Damage3, panel=function(x,y)
{
  panel.xyplot(x,y)
  panel.abline(lm(y~ x), lty=1, col=4)
},
strip= strip.custom( par.strip.text = list(cex=0.75)))


#model averaging
ms1 <- dredge(M1,REML=F,rank =AICc,evaluate = T,trace=T)
ms1
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
avgm <- model.avg(delta7,se.fit=T)


#predict from model averaged stuff
tapply(M_Damage2$Vol_log2,M_Damage2$Method,min)

Damage_pred<-data.frame(Vol_log2=c(seq(10,150,1),seq(5,97,1)),Method=c(rep("Conventional",141),rep("RIL",93)))


Damage_pred2<-predict(avgm,newdata=Damage_pred,level=0,se.fit=T)
Damage_pred2<-cbind(Damage_pred,Damage_pred2)
head(Damage_pred2)

#plot these results
theme_set(theme_bw(base_size=12))
Dam_1<-ggplot(M_Damage3,aes(x=Vol_log2,y=Prop_dam,colour=Method))+geom_point()+geom_line(data=Damage_pred2,aes(x=Vol_log2,y=plogis(fit),colour=Method),size=2)
Dam_2<-Dam_1+geom_line(data=Damage_pred2,aes(x=Vol_log2,y=plogis(fit+(1.96*se.fit)),colour=Method),lty=2)+geom_line(data=Damage_pred2,aes(x=Vol_log2,y=plogis(fit-(1.96*se.fit)),colour=Method),lty=2)
Dam_3<-Dam_2+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Dam_3+xlab("Volume logged per hectare")+ylab("Proportion of total tree stems damaged")

setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Prop_damaged_vol.jpeg",height=6,width=8,dpi=1200)


M1<-lmer(qlogis(Prop_dam)~Vol_log2*Method+I(Vol_log2^2)*Method+(Vol_log2|Study),M_Damage3,na.action="na.fail")

M2<-nlme(qlogis(Prop_dam)~SSasymp(Vol_log2,Asym, R0,lrc),data=M_Damage3,fixed=Asym + R0 + lrc ~ 1,random = R0 ~ 1,groups=Study ~ 1,start = c(Asym = 1, R0 =0, lrc = -3.3))
summary(M2)

plot(predict(M2),resid(M2))

plot(M_Damage3$Method,predict(M2))

