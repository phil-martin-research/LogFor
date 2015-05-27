######################################################
#start here just to plot the figure###################
######################################################

AGB_preds<-read.csv("Data/Preds_AGB.csv")
Rich_preds<-read.csv("Data/Preds_Richness.csv")
Rich<-read.csv("Data/Richness_studies.csv")
AGB<-read.csv("Data/AGB_studies_vol.csv")


library(ggplot2)
library(dplyr)

#organise data into one dataframe for predictions and one for data
Rich_preds$Type<-"Species richness"
AGB_preds$Type<-"Aboveground biomass"
newdat<-rbind(AGB_preds,Rich_preds)

keeps<-c("yi","vi","Method","Vol")
keeps2<-c("yi","vi","Method","Vol2")
Rich2<-Rich[,keeps]
AGB2<-AGB[,keeps]
AGB2$Vol<-AGB2$Vol2
AGB2$Vol2<-NULL
Rich2$Type<-"Species richness"
AGB2$Type<-"Aboveground biomass"
  
Data<-rbind(AGB2,Rich2)
#plot results
#first create x axis labels
label_loc<-data.frame(x=c(10,10),y=c(0.15,0.15),lab=c("(a)","(b)"),Type=c("Aboveground biomass","Species richness"))
Vol_ax<-(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))
theme_set(theme_bw(base_size=12))
vol_plot<-ggplot(newdat,aes(x=Vol,y=exp(yi)-1,ymax=exp(UCI)-1,ymin=exp(LCI)-1))+geom_ribbon(alpha=0.2)+geom_line(size=1)+facet_wrap(~Type,scale="free")
vol_plot2<-vol_plot+geom_point(data=Data,aes(ymax=NULL,ymin=NULL,colour=Method,size=1/vi),shape=1)
vol_plot3<-vol_plot2+ylab("Proportional change in metric following logging")
vol_plot4<-vol_plot3+scale_size_continuous(range=c(5,10))+geom_hline(y=0,lty=2,size=1)
vol_plot5<-vol_plot4+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
vol_plot6<-vol_plot5+xlab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+scale_colour_brewer(palette="Set1")
rich_vol_plot<-vol_plot6+theme(legend.position="none")+scale_colour_brewer(palette="Set1")
rich_vol_plot+geom_text(data = label_loc, aes(x=x,y=y,label=lab,ymax=NULL,ymin=NULL))
ggsave("Figures/AGB_Richness_volume.png",height=4,width=8,dpi=1200)
