#script to look at the reporting of variables in studies

#name: Phil Martin
#date:10/05/2015

#clear objects
rm(list=ls())
std <- function(x) sd(x)/sqrt(length(x))

#open packages
library(ggplot2)
library(plyr)
library(reshape2)

Methods<-read.csv("Data/Methods_reporting.csv",stringsAsFactors=FALSE)
Methods<-Methods[!duplicated(Methods[1:2]),]



#melt dataframe for use in ggplot
Meth_melt<-melt(Methods,id.vars=c("Study","Method"))
Meth_melt$variable<-gsub(".", " ", Meth_melt$variable, fixed = TRUE)
Meth_melt$variable<-ifelse(Meth_melt$variable=="Planned road building   skidder routes","Planned roads and skidder routes",Meth_melt$variable)
Meth_melt$variable<-ifelse(Meth_melt$variable=="Minimum DBH of trees cut  cm ","Minimum DBH",Meth_melt$variable)
Meth_melt$variable<-ifelse(Meth_melt$variable=="Illegal Legal","Illegal/legal",Meth_melt$variable)
Meth_melt$value<-ifelse(is.na(Meth_melt$value),0,1)
head(Meth_melt)
Method_df<-ddply(Meth_melt,.(variable,Method),summarise,Prop=(sum(value))/length(value))
ddply(Method_df,.(Method),summarise,Mean_prop=mean(Prop),SE=std(Prop))



#plot this result
theme_set(theme_bw(base_size=12))
Plot1<-ggplot(Method_df,aes(x=variable,y=Prop,fill=Method))+geom_bar(stat = "identity",position="dodge")+ theme(axis.text.x = element_text(angle = 90,hjust =1,vjust = 0.5))
Plot2<-Plot1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))                                                                                                                                          
Plot2+scale_fill_brewer(palette="Set1")+xlab("Variable")+ylab("Proportion of sites for which data was reported")
ggsave("Figures/Method_rep.png",dpi=400,height=8,width=6,units="in")
