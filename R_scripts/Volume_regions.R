#open packages
library(RODBC)
library(ggplot2)
library(metafor)
library(reshape2)

#connect to database
log <- odbcConnect("Logging")
sqlTables(log)

#import data
Loc<- sqlFetch(log, "Location")
head(Loc)
colnames(Loc)<-c("Study","Site","Method","Vol","Lat","Long")

#set contintent
Loc$Cont<-ifelse(Loc$Long<(-30),"Americas",NA)
Loc$Cont<-ifelse(Loc$Long>=-30&Loc$Long<=60,"Africa",Loc$Cont)
Loc$Cont<-ifelse(Loc$Long>60,"SE Asia",Loc$Cont)
Loc$Cont<-as.factor(Loc$Cont)

tapply(Loc2$Vol,Loc2$Cont,mean)



theme_set(theme_bw(base_size=26))
a<-ggplot(data=Loc,aes(x=Long,y=Lat))+borders("world", size=0.1,colour="grey",fill="lightgrey")+theme(panel.grid.major = element_line(colour =NA))+geom_point(size=2,alpha=0.5)
b<-a+coord_cartesian(xlim = c(-130, 180),ylim=c(-30, 30))
b+geom_vline(x=c(-30,60))

#how does intensity vary over the tropics
Loc2<-subset(Loc,Vol>0)

#standard error function
std <- function(x) sd(x)/sqrt(length(x))

#calculate mean

Vols<-cbind(melt(tapply(Loc2$Vol,Loc2$Cont,mean)),(melt(tapply(Loc2$Vol,Loc2$Cont,std))[2]))
colnames(Vols)<-c("Region","Vol","SE")

ggplot(data=Vols,aes(x=Region,y=Vol,ymax=Vol+(1.96*SE),ymin=Vol-(1.96*SE)))+geom_pointrange(size=1.5)


#map of regions
a<-ggplot(data=Loc2,aes(x=Long,y=Lat))+borders("world", size=0.1,colour="grey",fill="lightgrey")+theme(panel.grid.major = element_line(colour =NA))+geom_point(aes(colour=Vol),alpha=0.5)
b<-a+coord_cartesian(xlim = c(-130, 180),ylim=c(-30, 30))+theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())
b+scale_size_area()
