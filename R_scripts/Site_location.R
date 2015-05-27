#######################################################################################
#Script for plotting location of sites used in logging meta-analysis###############
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:18/05/2015

#clear objects
rm(list=ls())

#open packages
library(ggplot2)
library(reshape)
library(plyr)

###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################

Location<-read.csv("Data/Study_locations.csv")

#import locations
Location$LatR<-round(Location$Lat*4,-1)/4
Location$LongR<-round(Location$Long*4,-1)/4
Location

#subset to calculate number of points per grid cell for each map
Loc_Rich<-subset(Location,Type=="Species richness")
Loc_AGB<-subset(Location,Type=="Aboveground biomass")
Loc_Dam<-subset(Location,Type=="Residual damage")


#count number of points
Count_Rich<-count(Loc_Rich,vars=c("LongR","LatR"))
Count_AGB<-count(Loc_AGB,vars=c("LongR","LatR"))
Count_Dam<-count(Loc_Dam,vars=c("LongR","LatR"))
Count_Rich$Measure<-"Species richness"
Count_AGB$Measure<-"Aboveground biomass"
Count_Dam$Measure<-"Residual damage"

#merge back together
Loc_comb<-rbind(Count_Rich,Count_AGB,Count_Dam)

head(Loc_comb)

#plot locations with points
site_map_a<-ggplot(data=Loc_comb,aes(x=LongR, y=LatR,size=freq))+borders("world", size=0.2,colour="grey",fill="lightgrey")+theme_bw()+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA))+geom_point(alpha=0.5,colour="blue",shape=1)+coord_equal()
site_map_b<-site_map_a+coord_map(project="rectangular",lat0 = 0)+ coord_cartesian(xlim = c(-100, 150),ylim=c(-24, 24))
site_map_c<-site_map_b+theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())+facet_wrap(~Measure,ncol=1)+scale_size_area(name="No. of sites")
site_map_c
ggsave("Figures/Site_locations_point_log.png",height=5.5,width=8,dpi=600)

