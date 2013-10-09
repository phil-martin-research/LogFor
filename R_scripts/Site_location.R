#######################################################################################
#Script for plotting location of sites used in logging meta-analysis###############
#and plots for paper###################################################################
#######################################################################################

#name: Phil Martin
#date:09/10/2013

#open packages
library(RODBC)
library(ggplot2)
library(reshape)
library(plyr)

###########################################################################################
#Organise data before analysis#############################################################
###########################################################################################

#connect to database
log <- odbcConnect("Logging")
sqlTables(log)

#import data
Location<- sqlFetch(log, "Location")
head(Location)
colnames(Location)<-c("Study","Site","Method","Vol","Lat","Long","Measure")

#import locations
Location$LatR<-round(Location$Lat*4,-1)/4
Location$LongR<-round(Location$Long*4,-1)/4
Location

#subset to calculate number of points per grid cell for each map
Loc_Rich<-subset(Location,Measure=="Species richness")
Loc_AGB<-subset(Location,Measure=="Aboveground biomass")

#count number of points
Count_Rich<-count(Loc_Rich,vars=c("LongR","LatR"))
Count_AGB<-count(Loc_AGB,vars=c("LongR","LatR"))
Count_Rich$Measure<-"Species richness"
Count_AGB$Measure<-"Aboveground biomass"

#merge back together
Loc_comb<-rbind(Count_Rich,Count_AGB)

head(Loc_comb)

#plot locations with points
site_map_a<-ggplot(data=Loc_comb,aes(x=LongR, y=LatR,size=freq))+borders("world", size=0.1,colour="grey",fill="lightgrey")+theme_bw()+theme(panel.grid.major = element_line(colour =NA))+geom_point(alpha=0.5,colour="blue")
site_map_b<-site_map_a+coord_map(project="rectangular",lat0 = 0)+ coord_cartesian(xlim = c(-130, 160),ylim=c(-40, 40))
site_map_b+theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank())+facet_wrap(~Measure)+scale_area(range = c(2, 5),name="No. of sites")
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Site_locations_point_log.pdf",height=3,width=8,dpi=300)
