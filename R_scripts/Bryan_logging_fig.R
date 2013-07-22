library(reshape2)
library(ggplot2)
library(gridExtra)

Biomass<-seq(100,400,0.01)
mean_log<-data.frame(Biomass-exp(0.012*Biomass))/Biomass
mean_log$Biomass<-Biomass
mean_log$Type<-"Harvest"
colnames(mean_log)<-c("Mean","Biomass","Type")

mean_coll<-data.frame(Biomass-(exp(0.0265*Biomass)^0.5))/Biomass
mean_coll$Biomass<-Biomass
mean_coll$Type<-"Collateral"
colnames(mean_coll)<-c("Mean","Biomass","Type")

mean_def<-data.frame(Biomass-(exp(0.0268*Biomass)^0.5))/Biomass
mean_def$Biomass<-Biomass
mean_def$Type<-"Deforestation"
colnames(mean_def)<-c("Mean","Biomass","Type")

combined<-rbind(mean_log,mean_coll,mean_def)

as.factor(combined$Type)
1-max(sub_log$mean_coll)
length(rep(Biomass,3))

a<-ggplot(combined,aes(x=Biomass,y=Mean+0.032,colour=factor(Type)))+geom_line(size=1)
a
b<-a+theme_bw()+opts(legend.position = "none")+opts(panel.grid.major = theme_line(colour =NA))+opts(axis.title.x = theme_text(size = 20, colour = 'black'))
c<-b+ylab(NULL)+opts(axis.title.y = theme_text(size = 20, colour = 'black'))+xlab("Biomass prior to logging")
d<-c+theme(axis.text.x = theme_text(size = 12),axis.text.y = theme_text(size = 14),strip.text.x = theme_text(size = 16))
e<-d+coord_cartesian(xlim=c(100,400),ylim=c(0.4,1))
e
#save plot
setwd("C:/Documents and Settings/PMART/My Documents/Dropbox/Publications, Reports and Responsibilities/Chapters/5. Tropical Forest degradation/Analysis")
ggsave(filename="Logging biomass Bryan2.png",height=5,width=3,dpi=500)

#Possible effects of RIL
Median_log<-c(0.7,0.4)
Upper_log<-c(0.8,0.6)
Lower_log<-c(0.5,0.2)
Method<-c("RIL","Conventional")
RIL_comb<-data.frame(Median_log,Upper_log,Lower_log,Method)

RIL_plot<-ggplot(RIL_comb,aes(x=Method, y=Median_log,ymax=Upper_log,ymin=Lower_log,fill=factor(Method)))+geom_crossbar()
RIL_plot2<-RIL_plot+theme_bw()+opts(legend.position = "none")+opts(panel.grid.major = theme_line(colour =NA))+ylab("Proportion of biomass \nremaining after harvesting")+opts(axis.title.y = theme_text(size = 14, colour = 'black'))+xlab("Logging method")+opts(axis.title.x =theme_text(size = 14, colour = 'black'))
RIL_plot3<-RIL_plot2+theme(axis.text.x = theme_text(size = 12),axis.text.y = theme_text(size = 14))+coord_cartesian(,ylim=c(0,1))
RIL_plot3
ggsave(filename="RIL_biomass.png",height=5,width=2,dpi=500)
#effects of harvest intensity
HI<-c(19,28,40)
PL<-c(5,10,20)
Harvest<-data.frame(HI,PL)

H_plot<-ggplot(Harvest,aes(y=1-(PL/100),x=HI))+geom_line(size=1)+theme_bw()+opts(legend.position = "none")+opts(panel.grid.major = theme_line(colour =NA))+ylab(NULL)+opts(axis.title.y = theme_text(size = 14, colour = 'black'))+opts(axis.title.x =theme_text(size = 14, colour = 'black'))
H_plot2<-H_plot+theme(axis.text.x = theme_text(size = 12),axis.text.y = theme_text(size = 14))+coord_cartesian(ylim=c(0,1))
H_plot3<-H_plot2+xlab(expression(paste("Harvest intensity", "(", m^3, ha^-1,")")))+coord_cartesian(ylim=c(0.7,1),xlim=c(20,40))
H_plot3
ggsave(filename="Harvest intensity.png",height=5,width=2,dpi=500)

grid.arrange(RIL_plot3,H_plot3,e,ncol=1)
?grid.arrange
setwd("C:/Documents and Settings/Phil/My Documents/My Dropbox/Publications, Reports and Responsibilities/Chapters/5. Tropical Forest degradation/Analysis")
png("logging.png",height=1800,width=3600,res=300)
grid.arrange(RIL_plot3,H_plot3,e,nrow=1)
dev.off()



