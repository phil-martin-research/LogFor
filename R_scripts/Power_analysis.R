#######################################################################################
#Script for analyses associated with reviewers comments ###############################
#######################################################################################

#name: Phil Martin
#date:06/03/14

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(metafor)
library(GGally)
library(multcomp)
library(plyr)


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
AGB_NoSafe<-subset(AGB,N_Logged<2)



#calculate RR
ROM<-escalc(data=AGB,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
ROM_nosafe<-escalc(data=AGB_NoSafe,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)

ME_summary_Age<-rma.mv(yi,vi,random=~(1|Age)+(1|ID),method="REML",data=ROM)
ME_summary_Age_NS<-rma.mv(yi,vi,random=~(1|Age)+(1|ID),method="REML",data=ROM_nosafe)

#test for differences between RIL and conventional

RIL_ME<-rma.mv(yi,vi,mods=~Method-1,random=~(1|Age)+(1|ID),method="REML",data=ROM)
summary(RIL_ME)
RIL_ME_nosafe<-rma.mv(yi,vi,mods=~Method-1,random=~(1|Age)+(1|ID),method="REML",data=ROM_nosafe)
summary(RIL_ME_nosafe)



#now insert one row at a time adding an RIL study based on the mean values of
#previous studies

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

#work out mean values for RIL studies
ROM_RIL<-subset(ROM,Method=="RIL")

AGB2<-AGB[,c(-1:-2,-7,-10,-14,-15,-16,-19)]
AGB2_RIL<-subset(AGB2,Method=="RIL")
AGB2_RIL$SDU/AGB2_RIL$MU
UL_RIL<-runif(100000,min=252,max=420)
ULSD_RIL<-runif(100000,min=0.04,max=.25)*UL_RIL
L_RIL<-UL_RIL-abs(((UL_RIL/100)*rnorm(100000,27,10)))
LSD_RIL<-ULSD_RIL*runif(n=100000,min=min((AGB2_RIL$SDL)/(AGB2_RIL$SDU)),max=max((AGB2_RIL$SDL)/(AGB2_RIL$SDU)))
SSU<-round(runif(n=100000,min=3,max=20))
SSL<-round(SSU*runif(n=100000,min=.7,max=2))
Age<-round((abs(rnorm(100000,mean(AGB2$Age),sd=sd(AGB2$Age)))))
rep_row<-seq(from=1,to=1000)
rep_meta<-seq(from=0,to=50)
new.RIL<-data.frame(Age=Age,Method=as.factor("RIL"),Vol=-9999,MU=UL_RIL,SSU=SSU,
  ML=L_RIL,SSL=SSL,VarT="SD",ID=seq(from=37,to=100036),
  SDU=ULSD_RIL,SDL=LSD_RIL)

ptm <- proc.time()
for (i in 1:1000){
  RIL_sample<-new.RIL[sample(nrow(new.RIL),size=50,replace=F),]
  New_AGB<-rbind(AGB2,RIL_sample)
  ROM<-escalc(data=New_AGB,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)
  Results<-data.frame(estimate=numeric(),ci.lb=numeric(),ci.ub=numeric(),add=numeric(),Method=factor(),rep=as.numeric())
  setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables/Power analysis")
  write.csv(Results,paste("Results_",rep_row[i],".csv",sep=""),row.names=F)
  for (y in 1:51){
    setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables/Power analysis")
    Results<-read.csv(paste("Results_",rep_row[i],".csv",sep=""))
    RIL_ME<-rma.mv(yi,vi,mods=~Method-1,random=~(1|Age)+(1|ID),method="REML",data=ROM[1:(37+rep_meta[y]),])
    newrow<-data.frame(estimate=coef(summary(RIL_ME))$estimate,ci.lb=coef(summary(RIL_ME))$ci.lb,
    ci.ub=coef(summary(RIL_ME))$ci.ub,add=rep(rep_meta[y],2),Method=as.factor(c("Conventional","RIL")),rep=rep_row[i])
    Results2<-rbind(Results,newrow)
    write.csv(Results2,paste("Results_",rep_row[i],".csv",sep=""),row.names=F)
  }
}
proc.time() - ptm
warnings()

  
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Tables/Power analysis")

All_meta<-do.call(rbind.fill,lapply(list.files("."),FUN=read.csv))

BS_MA<-data.frame(med=with(All_meta,tapply(estimate,list(Method,add),median))[1:102],
ci.lb=with(All_meta,tapply(ci.lb,list(Method,add),median))[1:102],
ci.ub=with(All_meta,tapply(ci.ub,list(Method,add),median))[1:102],
Method=rep(c("Conventional","RIL"),times=51),add=rep(0:50,each=2))

theme_set(theme_bw(base_size=14))
a<-ggplot(data=All_meta,aes(x=add,ymin=exp(ci.lb)-1,ymax=exp(ci.ub)-1,group=interaction(Method,rep),fill=Method))+geom_ribbon(alpha=0.005,lty=2)+geom_line(aes(y=exp(estimate)-1),alpha=0.01)
b<-a+geom_line(data=BS_MA,aes(y=exp(ci.ub)-1,x=add,group=Method,colour=Method),lty=2,size=.5,fill=NA,alpha=0.8)
c<-b+geom_line(data=BS_MA,aes(y=exp(ci.lb)-1,x=add,group=Method,colour=Method),lty=2,size=.5,fill=NA,alpha=0.8)
d<-c+geom_line(data=BS_MA,aes(y=exp(med)-1,x=add,group=Method,colour=Method),size=1,alpha=0.8)
e<-d+theme(legend.position="none")+scale_size_continuous(range=c(0.5,1.5))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA),axis.title=element_text(size=14,face="bold"),title=element_text(size=14,face="bold"))
f<-e+ylab("Proportional change following logging")+xlab("Number of additional studies")+geom_hline(x=0,lty=2,size=1)
f+coord_cartesian(xlim=c(0,50))+scale_colour_brewer(palette="Set1")+scale_fill_brewer(palette="Set1")
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave(filename="Power_analysis.pdf",width=8,height=4,units="in",dpi=300)

?scale_colour_brewer()

df<-data.frame

for (i in 1:2){
  
}





