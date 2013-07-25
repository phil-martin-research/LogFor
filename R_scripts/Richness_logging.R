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

#forrest plot of this
theme_set(theme_bw(base_size=26))
forrest_data<-rbind(data.frame(ES=ROM.ma$yi,SE=sqrt(ROM.ma$vi),Type="Site",Study=Richness$Site),data.frame(ES=ROM.ma$b,SE=ROM.ma$se,Type="Summary",Study="Summary"))
forrest_data$Study2<-factor(forrest_data$Study, levels=rev(levels(forrest_data$Study)) )
levels(forrest_data$Study2)
plot1<-ggplot(data=forrest_data,aes(x=Study2,y=exp(ES)-1,ymax=exp(ES+(1.96*SE))-1,ymin=exp(ES-(1.96*SE))-1,size=(1/SE)/10,colour=factor(Type)))+geom_pointrange(shape=15)
plot2<-plot1+coord_flip()+geom_hline(aes(x=0), lty=2,size=1)
plot3<-plot2+xlab("Study")+ylab("Proportional change")+scale_colour_manual(values=c("grey","black"))
plot3+theme(legend.position="none")+scale_size_continuous(range=c(1,4))

#model to look at differences between rarefied and not rarefied
Model1<-rma.uni(yi,vi,mods=~Method+Rare+Tax,method="ML",data=ROM)
summary(Model1)

Model2<-rma.uni(yi,vi,mods=~Rare+Tax,method="ML",data=ROM)
summary(Model2)

Model3<-rma.uni(yi,vi,mods=~Rare,method="ML",data=ROM)
summary(Model3)


model.res<-data.frame(Mean=c(-0.2943011,-0.2943011+0.3006),SE=c(0.0528,0.0730),Method=c("Unrarefied","Rarefied"))

#plot the results of differences in rarefied
Rare_plot<-ggplot(data=model.res,aes(y=Mean,x=Method,ymax=Mean+(SE*1.96),ymin=Mean-(SE*1.96)))+geom_pointrange(size=1.5)
Rare_plot+coord_flip()+geom_hline(y=0,lty=2)

#plot the results of differences in taxa
Rare_plot<-ggplot(data=model.res,aes(y=Mean,x=Method,ymax=Mean+(SE*1.96),ymin=Mean-(SE*1.96)))+geom_pointrange(size=1.5)
Rare_plot+coord_flip()+geom_hline(y=0,lty=2)


#subset to only include studies with data on volume
Rich_vol<-subset(Richness,Vol>-9999)
Rich_vol
Rich_vol<-subset(Rich_vol,Rare=="R")
#calculate the log ratio for these studies
ROM_Rich<-escalc(data=Rich_vol,measure="ROM",m2i=MU,sd2i=SDU,n2i=SSU,m1i=ML,sd1i=SDL,n1i=SSL,append=T)


#model to look at importance of volume
Vol_M1<-rma.uni(yi,vi,mods=~Method+Vol,method="ML",data=ROM_Rich)
summary(Vol_M1)

Vol_M2<-rma.uni(yi,vi,mods=~Vol,method="ML",data=ROM_Rich)
summary(Vol_M2)

Vol_M3<-rma.uni(yi,vi,mods=~Rare+I(Vol/MU),method="ML",data=ROM_Rich)
summary(Vol_M3)

coef(Vol_M3)

plot(Rich_vol$Vol,ROM_Rich$yi,col=Rich_vol$Rare)



