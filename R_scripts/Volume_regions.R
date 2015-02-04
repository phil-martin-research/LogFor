#######################################################################################
#Script for plots showing variation in volume by location and for different############
#logging methods#######################################################################
#######################################################################################

#name: Phil Martin
#date:04/02/2015

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(plyr)


###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################


setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Vol<-read.csv("Vol_var.csv")
Vol2<-unique(Vol)

head(Vol2)


#plot variation in volume by region and by method
theme_set(theme_bw(base_size=25))
Vol_plot<-ggplot(Vol2,aes(y=Volume_pred,x=Region,fill=Method))+geom_boxplot()
Vol_plot2<-Vol_plot+ylab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Vol_plot2+scale_fill_brewer(palette="Set1")
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Vol_regions.png",height=6,width=10,dpi=400)

tapply(x=Vol2$Volume_pred,INDEX =Vol2$Region,median,na.rm=TRUE)
v<-c("Volume_pred")
sapply(v, function(i) tapply(Vol2[[i]], Vol2$Region, mean, na.rm=TRUE))
