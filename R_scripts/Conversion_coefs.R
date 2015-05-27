#############################################################
#script to calculate scaling coeffients between different####
#measures of logging intensity and damage####################
#############################################################

rm(list=ls())

#first load packages
library(ggplot2)
library(MuMIn)
library(lme4)
library(plyr)
library(reshape2)

#load data
#for Proportional damage
Prop_dam<-read.csv("Data/prop_damage_edit2.csv")
#for data on trees logged vs volume logged relationships
Vol_tree<-read.csv("Data/Vol_trees.csv")
head(Prop_dam)

length(which(!is.na(Prop_dam$Damage_prop)))

colnames(Prop_dam)<-c("Study","Site_ID","ID","Method","Vol","Tree_ex","Dam_ha","Prop_dam","Region","Plot","Replicates")

################################################################
#model of volume as a function of number of trees logged per ha#
################################################################
#code africa and asia as one region
Vol_tree<-unique( Vol_tree[ , 1:4 ] )
OldVar<-Vol_tree$Region
NewVar <- factor(rep(NA, length(OldVar) ), 
                 levels=c("Africa/Asia", "Americas") )
NewVar[ OldVar %in% c("Africa", "Asia")] <- "Africa/Asia"
NewVar[ OldVar %in% c("Americas")]<-"Americas"
Vol_tree$Region2<-NewVar
ggplot(Vol_tree,aes(x=Trees2,y=Vol2,colour=Region2,shape=Region2))+geom_point()+geom_smooth(method=lm,se=F)



#create dataset with only complete cases in explanatory variables
Vol_tree_CC<-Vol_tree[complete.cases(Vol_tree[,c(2,3)]),]
Vol_tree_CC$Tree_corr<-(Vol_tree_CC$Trees2-mean(Vol_tree_CC$Trees2,na.rm = T))/sd(Vol_tree_CC$Trees2,na.rm = T)
Vol_tree_CC<-subset(Vol_tree_CC,Trees2<25)

M1<-lmer(Vol2~Tree_corr*Region2-1+(Tree_corr|Study),Vol_tree_CC,na.action="na.fail",REML=F,)
M2<-lmer(Vol2~Tree_corr*Region2+(1|Study),Vol_tree_CC,na.action="na.fail",REML=F)


#look at diagnostic plots
plot(M1)
plot(Vol_tree_CC$Tree_corr,resid(M1))
plot(Vol_tree_CC$Vol2,predict(M1),col=Vol_tree_CC$Region2)


#model averaging
ms1 <- dredge(M2,REML=F,rank =AICc,evaluate = T)
delta7<- get.models(ms1, subset = cumsum(delta) <= 7)
Model_sel<-model.sel(delta7)
avgm <- model.avg(delta7)
summary(avgm)
r.squaredGLMM(M2)

summary(M2)


Vol_tree_CC$pred<-predict(M2)
#create predictions from model
tapply(Vol_tree_CC$Trees2,Vol_tree_CC$Region2,summary)

Tree_pred1<-data.frame(Trees2=c(seq(0.4,13.06,length.out = 500),seq(1.85,16,length.out = 500)),Region2=c(rep("Africa/Asia",500),rep("Americas",500)))
Tree_pred1$Tree_corr<-(Tree_pred1$Trees2-mean(Vol_tree_CC$Trees2,na.rm = T))/sd(Vol_tree_CC$Trees2,na.rm = T)
Tree_pred1$Pred<-predict(avgm,newdata=Tree_pred1,re.form=NA)


#predict values for plotting
#now create plots of this
newdat<-data.frame(Trees2=c(seq(0.4,13.06,length.out = 500),seq(1.85,16,length.out = 500)),Region2=c(rep("Africa/Asia",500),rep("Americas",500)))
newdat$Tree_corr<-(newdat$Trees2-mean(Vol_tree_CC$Trees2,na.rm = T))/sd(Vol_tree_CC$Trees2,na.rm = T)
newdat$Vol2<-0
mm <- model.matrix(terms(M2),newdat)
newdat$Vol2 <- predict(M2,newdat,re.form=NA)
pvar1 <- diag(mm %*% tcrossprod(vcov(M2),mm))
tvar1 <- pvar1+VarCorr(M2)$Study[1]  ## must be adapted for more complex models
tvar1 <- 
  newdat <- data.frame(
    newdat
    , plo = newdat$Vol2-2*sqrt(pvar1)
    , phi = newdat$Vol2+2*sqrt(pvar1)
    , tlo = newdat$Vol2-2*sqrt(tvar1)
    , thi = newdat$Vol2+2*sqrt(tvar1)
  )

newdat

#plot this prediction
theme_set(theme_bw(base_size=12))
Scaling_plot1<-ggplot(Vol_tree_CC,aes(x=Trees2,y=Vol2,colour=Region2,group=Study))+geom_point(size=3,shape=1,alpha=0.5)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Scaling_plot2<-Scaling_plot1+xlab(expression(paste("Number of trees extracted ",ha^-1)))+ylab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
Scaling_plot2+geom_line(data=newdat,aes(x=Trees2,y=Vol2,group=Region2,colour=Region2),size=2)+
  geom_ribbon(data=newdat,aes(x=Trees2,y=NULL,ymax=thi,ymin=tlo,fill=Region2,group=Region2,colour=NULL),alpha=0.2)+
  scale_colour_brewer("Region",palette="Set1")+scale_fill_brewer("Region",palette="Set1")+theme(legend.position="none")
ggsave("Figures/Volume_trees.png",height=6,width=8,dpi=400)

#################################################
#now use predictions from this model in other ###
#datasets########################################
#################################################

#first the data for proportion of stems damaged
#code africa and asia as one region
OldVar<-Prop_dam$Region
NewVar <- factor(rep(NA, length(OldVar) ), 
                 levels=c("Africa/Asia", "Americas") )
NewVar[ OldVar %in% c("Africa", "Asia")] <- "Africa/Asia"
NewVar[ OldVar %in% c("Americas")]<-"Americas"
Prop_dam$Region2<-NewVar


Tree_pred2<-data.frame(ID=Prop_dam$ID,Tree_corr=(Prop_dam$Tree_ex-mean(Vol_tree_CC$Trees2,na.rm = T))/sd(Vol_tree_CC$Trees2,na.rm = T),Region2=Prop_dam$Region2)
Tree_pred2<-Tree_pred2[complete.cases(Tree_pred2),]
Tree_pred2$Pred<-(predict(avgm,newdata=Tree_pred2,re.form=NA))
M_Damage<-merge(Prop_dam,Tree_pred2,by="ID",all=T)
#put in predicted volume values where there are none where number of treesextracted is <25 per ha
M_Damage$Vol2<-NULL
for (i in 1:nrow(M_Damage)){
  if(is.na(M_Damage$Vol[i])&M_Damage$Region2.x[i]=="Africa/Asia"&M_Damage$Tree_ex[i]<12&!is.na(M_Damage$Tree_ex[i])){
    M_Damage$Vol2[i]<-M_Damage$Pred[i]
  }else if (is.na(M_Damage$Vol[i])&M_Damage$Region2.x[i]=="Americas"&M_Damage$Tree_ex[i]<11&!is.na(M_Damage$Tree_ex[i])){
    M_Damage$Vol2[i]<-M_Damage$Pred[i]
  }else {M_Damage$Vol2[i]<-M_Damage$Vol[i]
  }
}

ggplot(M_Damage,aes(x=Tree_ex,fill=is.na(Vol2)))+geom_histogram()+facet_wrap(~Region2.x)

length(na.omit(M_Damage$Vol2))
length(na.omit(M_Damage$Vol))

str(M_Damage)
keeps<-c("ID","Study","Site_ID","Method","Vol2","Tree_ex","Prop_dam","Region","Plot","Replicates")
M_Damage<-M_Damage[,keeps,drop=FALSE]
write.csv(M_Damage,"Data/Prop_damage.csv")



#now do the same for richness
Rich<-read.csv("Data/Rich_intens.csv")
Rich$Tree_corr<-(Rich$Trees-mean(Vol_tree_CC$Trees2,na.rm = T))/sd(Vol_tree_CC$Trees2,na.rm = T)
OldVar<-Rich$Region
NewVar <- factor(rep(NA, length(OldVar) ), 
                 levels=c("Africa/Asia", "Americas") )
NewVar[ OldVar %in% c("Africa", "Asia")] <- "Africa/Asia"
NewVar[ OldVar %in% c("Americas")]<-"Americas"
Rich$Region2<-NewVar
Rich$Pred<-(predict(avgm,newdata=Rich,re.form=NA))
Rich$Vol2<-NULL
for (i in 1:nrow(Rich)){
  if(is.na(Rich$Vol[i])&Rich$Region2[i]=="Africa/Asia"&Rich$Trees[i]<12&!is.na(Rich$Trees[i])){
    Rich$Vol2[i]<-Rich$Pred[i]
  }else if (is.na(Rich$Vol[i])&Rich$Region2[i]=="Americas"&Rich$Trees[i]<11&!is.na(Rich$Trees[i])){
    Rich$Vol2[i]<-Rich$Pred[i]
  }else {Rich$Vol2[i]<-Rich$Vol[i]
  }
}
head(Rich)

write.csv(Rich,"Data/Richness.csv")


#now do the same for biomass
AGB<-read.csv("Data/AGB_intens.csv")
AGB$Tree_corr<-(AGB$Trees-mean(Vol_tree_CC$Trees2,na.rm = T))/sd(Vol_tree_CC$Trees2,na.rm = T)
OldVar<-AGB$Region
NewVar <- factor(rep(NA, length(OldVar) ), 
                 levels=c("Africa/Asia", "Americas") )
NewVar[ OldVar %in% c("Africa", "Asia")] <- "Africa/Asia"
NewVar[ OldVar %in% c("Americas")]<-"Americas"
AGB$Region2<-NewVar
AGB$Pred<-(predict(avgm,newdata=AGB,re.form=NA))
AGB$Vol2<-NULL
for (i in 1:nrow(AGB)){
  if(is.na(AGB$Vol[i])&AGB$Region2[i]=="Africa/Asia"&AGB$Trees[i]<12&!is.na(AGB$Trees[i])){
    AGB$Vol2[i]<-AGB$Pred[i]
  }else if (is.na(AGB$Vol[i])&AGB$Region2[i]=="Americas"&AGB$Trees[i]<11&!is.na(AGB$Trees[i])){
    AGB$Vol2[i]<-AGB$Pred[i]
  }else {AGB$Vol2[i]<-AGB$Vol[i]
  }
}

tail(AGB)
write.csv(AGB,"Data/AGB.csv")



