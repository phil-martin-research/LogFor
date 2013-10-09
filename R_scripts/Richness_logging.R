#open packages
library(RODBC)
library(ggplot2)
library(metafor)
library(MuMIn)


#connect to database
log <- odbcConnect("Logging")
sqlTables(log)

#import data
Richness<- sqlFetch(log, "Richness query")
head(Richness)
colnames(Richness)<-c("Study","Site","Age","Method","Vol","Type","MU","VU","SSU","ML","VL","SSL","Vtype","Tax","Rar","Region")

#recalculate SDs

#unlogged
Richness$SDU<-ifelse(Richness$Vtype=="SE",Richness$VU*sqrt(Richness$SSU),Richness$VU)
Richness$SDU<-ifelse(Richness$Vtype=="CI",(Richness$VU/1.96)*sqrt(Richness$SSU),Richness$SDU)

#logged
Richness$SDL<-ifelse(Richness$Vtype=="SE",Richness$VL*sqrt(Richness$SSL),Richness$VL)
Richness$SDL<-ifelse(Richness$Vtype=="CI",(Richness$VL/1.96)*sqrt(Richness$SSL),Richness$SDL)
head(Richness)

#impute missing standard deviation values
Richness2<-subset(Richness,Richness$SDU>0)
head(Richness2)
Imp_U<-(sum(Richness2$SDU))/(sum(Richness2$MU))
Imp_L<-(sum(Richness2$SDL))/(sum(Richness2$ML))

Richness$SDU<-ifelse(Richness$SDU<0,Richness$MU*Imp_U,Richness$SDU)
Richness$SDL<-ifelse(Richness$SDL<0,Richness$ML*Imp_L,Richness$SDL)

#change rarefied column
Richness$Rare<-ifelse(Richness$Rar=="Not rarefied","NR",NA)
Richness$Rare<-ifelse(Richness$Rar=="Rarefied - Area","R",Richness$Rare)
Richness$Rare<-ifelse(Richness$Rar=="Rarefied - Individuals","R",Richness$Rare)
Richness$Rare<-as.factor(Richness$Rare)

#calculate the log ratio
ROM<-escalc(data=Richness,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)

###############################################
#summary analysis##############################
###############################################

#random effects meta-analysis
ROM.ma<-rma.uni(yi,vi,method="REML",data=ROM)
summary(ROM.ma)

exp(-.1275)-1
exp(-.2016)-1
exp(-.0533)-1


#forrest plot of this
theme_set(theme_bw(base_size=10))
forrest_data<-rbind(data.frame(ES=ROM.ma$yi,SE=sqrt(ROM.ma$vi),Type="Site",Study=Richness$Site),data.frame(ES=ROM.ma$b,SE=ROM.ma$se,Type="Summary",Study="Summary"))
forrest_data$Study2<-factor(forrest_data$Study, levels=rev(levels(forrest_data$Study)) )
levels(forrest_data$Study2)
plot1<-ggplot(data=forrest_data,aes(x=Study2,y=exp(ES)-1,ymax=exp(ES+(1.96*SE))-1,ymin=exp(ES-(1.96*SE))-1,size=(1/SE)/10,colour=factor(Type)))+geom_pointrange(shape=15)
plot2<-plot1+coord_flip()+geom_hline(aes(x=0), lty=2,size=1)
plot3<-plot2+xlab("Study")+ylab("Proportional change")+scale_colour_manual(values=c("grey","black"))
plot3+theme(legend.position="none")+scale_size_continuous(range=c(0.5,1.5))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Forest_Richness.pdf",height=4,width=6,dpi=1200)

##################################################################
#model to test for taxa, region and method differences############
##################################################################

Model1<-rma.uni(yi,vi,mods=~Rare+Tax+Region,method="ML",data=ROM)
Model2<-rma.uni(yi,vi,mods=~Rare+Tax,method="ML",data=ROM)
Model3<-rma.uni(yi,vi,mods=~Rare+Region,method="ML",data=ROM)
Model4<-rma.uni(yi,vi,mods=~Rare-1,method="ML",data=ROM)
Model5<-rma.uni(yi,vi,mods=~Region+Tax,method="ML",data=ROM)
Model6<-rma.uni(yi,vi,mods=~Region,method="ML",data=ROM)
Model7<-rma.uni(yi,vi,mods=~Tax,method="ML",data=ROM)
Model8<-rma.uni(yi,vi,mods=~1,method="ML",data=ROM)

#work out model AICc
Model_AIC<-data.frame(AIC=c(AIC(Model1),AIC(Model2),AIC(Model3),AIC(Model4),AIC(Model5),AIC(Model6),AIC(Model7),AIC(Model8)))

Model_AIC$Vars<-c("Rarefied+Region+Taxonomic group",
                   "Rarefied+Taxonomic group",
                   "Rarefied+Region","Rarefied",
                   "Region+Taxonomic group", "Region","Taxonomic group","Null")

#reorder from lowest to highest
Model_AIC<-Model_AIC[order(Model_AIC$AIC),]
#calculate AICc delta
Model_AIC$delta<-Model_AIC$AIC-Model_AIC$AIC[1]
#drop last models with delta >7
AIC_sel<-subset(Model_AICc,delta<=7)
#calculate the relative likelihood of model
AIC_sel$rel_lik<-exp((AIC_sel$AICc[1]-AIC_sel$AICc)/2)
#calculate the AICc weight
AIC_sel$weight<-AIC_sel$rel_lik/(sum(AIC_sel$rel_lik))

#dummy coding for variable importance
AIC_sel$Rare<-c(1,0,1,1,0,0)
AIC_sel$Region<-c(0,0,0,1,1,0)
AIC_sel$Tax<-c(0,0,1,0,0,1)

sum(AIC_sel$Rare*AIC_sel$weight)
sum(AIC_sel$Region*AIC_sel$weight)
sum(AIC_sel$Tax*AIC_sel$weight)

summary(Model4)


AIC_sel$R_squared<-c(1-(Model4$tau2/Model8$tau2),1-(Model8$tau2/Model8$tau2),1-(Model2$tau2/Model8$tau2),1-(Model3$tau2/Model8$tau2),1-(Model6$tau2/Model8$tau2),1-(Model7$tau2/Model8$tau2))

#plot results from rarefied model

Rare_model<-rma.uni(yi,vi,mods=~Rare-1,method="REML",data=ROM)

Coef_Rich<-coef(summary(Rare_model))

#put the rarefied results into data frame
model.rare<-data.frame(Mean=c(Coef_Rich[1,1],Coef_Rich[2,1]),SE=c(Coef_Rich[1,2],Coef_Rich[2,2]),Method=c("Unrarefied","Rarefied"))

#plot the results of differences in rarefied
theme_set(theme_bw(base_size=10))
Rare_plot<-ggplot(data=model.rare,aes(y=Mean,x=Method,ymax=Mean+(SE*1.96),ymin=Mean-(SE*1.96)))+geom_pointrange(size=1.5)
Rare_plot2<-Rare_plot+coord_flip()+geom_hline(y=0,lty=2)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Rare_plot2+xlab("Method of estimating species richness")+ylab("Proportional change in species richness \nfollowing selective logging")
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Richness_Rarefied.pdf",height=4,width=6,dpi=1200)



#meta-analysis looking at variation by region

Model1<-rma.uni(yi,vi,mods=~Region-1+Rare+Tax,method="ML",data=ROM)
summary(Model1)
Model2<-rma.uni(yi,vi,mods=~Region-1,method="ML",data=ROM)
summary(Model2)

AIC(Model1)
AIC(Model2)

#plot graph
Region_rich<-data.frame(coef(summary(Model1)))
Region_rich$Region<-as.factor(c("Africa","C & S America", "SE Asia & Australasia"))

a<-ggplot(data=Region_rich,aes(x=Region,y=estimate,ymin=ci.lb,ymax=ci.ub))+geom_pointrange(size=1)
b<-a+geom_hline(x=0,lty=2,size=1)+ylab("Proportional change after logging")+xlab("Regions")
b+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1,colour="black",fill=NA))

#Do analysis for volume

Rich_vol<-subset(Richness,Vol!=-9999)

#calculate the log ratio
ROM_vol<-escalc(data=Rich_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)

summary(Rich_vol$Vol)

Vol<-seq(10,160,1)

plot(ROM_vol$Vol,ROM_vol$yi,col=ROM_vol$Rare)
lines(Vol,-0.9612+(Vol*0.0051))

plot(Model1_Vol)


Model1_Vol<-rma.uni(yi,vi,mods=~Rare,method="ML",data=ROM_vol)
Model2_Vol<-rma.uni(yi,vi,mods=~Vol,method="ML",data=ROM_vol)


AIC(Model1_Vol)
AIC(Model2_Vol)

summary(Model1_Vol)