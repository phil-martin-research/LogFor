#######################################################################################
#Script for plotting location of sites used in logging meta-analysis###############
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:09/10/2013

#clear objects
rm(list=ls())

#open packages
library(RODBC)
library(ggplot2)
library(reshape)
library(plyr)

###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################

setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/Data/Fo analysis")
Location<-read.csv("Logged_locations.csv")

#import locations
Location$LatR<-round(Location$Lat*4,-1)/4
Location$LongR<-round(Location$Long*4,-1)/4
Location

#subset to calculate number of points per grid cell for each map
Loc_Rich<-subset(Location,Type=="Species richness")
Loc_AGB<-subset(Location,Type=="Aboveground biomass")

#count number of points
Count_Rich<-count(Loc_Rich,vars=c("LongR","LatR"))
Count_AGB<-count(Loc_AGB,vars=c("LongR","LatR"))
Count_Rich$Measure<-"Species richness"
Count_AGB$Measure<-"Aboveground biomass"

#merge back together
Loc_comb<-rbind(Count_Rich,Count_AGB)

head(Loc_comb)

#plot locations with points
site_map_a<-ggplot(data=Loc_comb,aes(x=LongR, y=LatR,size=freq))+borders("world", size=0.1,colour="grey",fill="lightgrey")+theme_bw()+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA))+geom_point(alpha=0.5,colour="blue")+coord_equal()
site_map_b<-site_map_a+coord_map(project="rectangular",lat0 = 0)+ coord_cartesian(xlim = c(-130, 160),ylim=c(-40, 40))
site_map_c<-site_map_b+theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())+facet_wrap(~Measure,ncol=1)+scale_area(range = c(2, 10),name="No. of sites")
site_map_c
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Site_locations_point_log.pdf",height=5.5,width=12,dpi=1200)

Location2<-subset(Location,Vol>0)
#plot how volume varies with location
#first get number of sites per region
xlabels <- ddply(Location2, .(Region), summarize, 
                 xlabels = paste(unique(Region), '\n(n = ', length(Region),')'))
theme_set(theme_bw(base_size=12))
a<-ggplot(Location2,aes(x=Region,y=Vol))+geom_boxplot()
b<-a+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
b+ylab(expression(paste("Volume of wood logged (",m^3,ha^-1,")")))+xlab("Region")+scale_x_discrete(labels=xlabels[['xlabels']])
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Region_vol.pdf",height=5,width=8,dpi=1200)
