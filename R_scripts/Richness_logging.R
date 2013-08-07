#open packages
library(RODBC)
library(ggplot2)
library(metafor)

#connect to database
log <- odbcConnect("Logging")
sqlTables(log)

#import data
Richness<- sqlFetch(log, "Richness query")
head(Richness)
colnames(Richness)<-c("Study","Site","Age","Method","Vol","Type","MU","VU","SSU","ML","VL","SSL","Vtype","Tax","Rar")

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

#log ratio
ROM<-escalc(data=Richness,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)


#this runs a random effects meta-analysis for log ratio data
ROM.ma<-rma.uni(yi,vi,method="REML",data=ROM)
summary(ROM.ma)

exp(-.1273)-1
exp(-.04)-1



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
ggsave("Forrest_Richness.jpeg",height=4,width=6,dpi=1200)


#model to look at differences between rarefied and not rarefied
Model1<-rma.uni(yi,vi,mods=~Method+Rare+Tax,method="ML",data=ROM)
summary(Model1)

Model2<-rma.uni(yi,vi,mods=~Rare+Tax,method="ML",data=ROM)

AIC(Model1)
AIC(Model2)

#model 2 looks to be the best model
#now lets use REML to get unbiased estimates of coefficients
Model2<-rma.uni(yi,vi,mods=~Rare+Tax,method="REML",data=ROM)
Coef_Rich<-coef(summary(Model2))
Coef_Rich
Coef_Rich[1,1]

#put the rarefied results into data frame
model.rare<-data.frame(Mean=c(Coef_Rich[1,1],Coef_Rich[1,1]+Coef_Rich[2,1]),SE=c(Coef_Rich[1,2],Coef_Rich[2,2]),Method=c("Unrarefied","Rarefied"))

#plot the results of differences in rarefied
Rare_plot<-ggplot(data=model.rare,aes(y=Mean,x=Method,ymax=Mean+(SE*1.96),ymin=Mean-(SE*1.96)))+geom_pointrange(size=1.5)
Rare_plot2<-Rare_plot+coord_flip()+geom_hline(y=0,lty=2)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Rare_plot2+xlab("Method of estimating species richness")+ylab("Proportional change in species richness \nfollowing selective logging")
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Richness_Rarefied.jpeg",height=4,width=6,dpi=1200)



#subset to only include studies with data on volume
Rich_vol<-subset(Richness,Vol>-9999)
Rich_vol
Rich_vol<-subset(Rich_vol,Rare=="R")
#calculate the log ratio for these studies
ROM_Rich<-escalc(data=Rich_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)

#model to look at importance of volume
Vol_M1<-rma.uni(yi,vi,mods=~Vol+Tax,method="ML",data=ROM_Rich)
summary(Vol_M1)
plot(Vol_M1)


Vol_M2<-rma.uni(yi,vi,mods=~Vol+Tax,method="ML",data=ROM_Rich)
summary(Vol_M2)
plot(Vol_M2)


Vol_M3<-rma.uni(yi,vi,mods=~I(Vol/MU),method="ML",data=ROM_Rich)
summary(Vol_M3)

plot(Rich_vol$Vol,ROM_Rich$yi,col=Rich_vol$Method)



