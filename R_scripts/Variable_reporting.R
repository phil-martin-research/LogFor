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
library(data.table)

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


#look at variation in methods
#Liana extraction
head(Methods)
for (i in 6:14){
Methods[,ncol(Methods)+1]<-ifelse(Methods[,i]=="Yes",1,NA)
Methods[,ncol(Methods)]<-ifelse(Methods[,i]=="No",0,Methods[,ncol(Methods)])
Methods[,ncol(Methods)]<-ifelse(is.na(Methods[,i]),0,Methods[,ncol(Methods)])
Methods[,ncol(Methods)]<-ifelse(is.na(Methods[,ncol(Methods)]),0,Methods[,ncol(Methods)])
}
head(Methods)

setnames(Methods, old = c('V15','V16',"V17","V18","V19","V20","V21","V22","V23"), 
         new = c('Lianas','P_roads',"P_Ext","DF","L_Mach","Train","Super","Silvi","Slopes"))

Methods$All<-ifelse((Methods$Lianas+Methods$P_roads+Methods$P_Ext+Methods$DF+Methods$L_Mach+Methods$Train+
                      Methods$Super+Methods$Slopes)>0,1,0)
Methods$All2<-Methods$Lianas+Methods$P_roads+Methods$P_Ext+Methods$DF+Methods$L_Mach+Methods$Train+
                       Methods$Super+Methods$Slopes


Method_var<-ddply(Methods,.(Method),summarise,
      Prop_liana=round((sum(Lianas,na.rm = T)/sum(!is.na(Lianas)))*100,1),
      Prop_road=round((sum(P_roads,na.rm = T)/sum(!is.na(P_roads)))*100,1),
      Prop_ext=round((sum(P_Ext,na.rm = T)/sum(!is.na(P_Ext)))*100,1),
      Prop_DF=round((sum(DF,na.rm = T)/sum(!is.na(DF)))*100,1),
      Prop_L_Mach=round((sum(L_Mach,na.rm = T)/sum(!is.na(L_Mach)))*100,1),
      Prop_train=round((sum(Train,na.rm = T)/sum(!is.na(Train)))*100,1),
      Prop_super=round((sum(Super,na.rm = T)/sum(!is.na(Super)))*100,1),
      Prop_slopes=round((sum(Slopes,na.rm = T)/sum(!is.na(Slopes)))*100,1)
      #Prop_all=round((sum(All,na.rm = T)/sum(!is.na(All)))*100,1),
      #Mean_all=mean(All2),se_all=std(All2)
      )


#plot these results
Method_var_melt<-melt(Method_var)
Method_var_melt$variable<-as.character(Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_liana","Lianas cut",Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_road","Road/skidder route planning",Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_ext","Planned extraction",Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_DF","Directional felling",Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_L_Mach","Lighter machinery",Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_train","Training of loggers",Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_super","Supervision of loggers",Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_silv","Silvicultural techniques",Method_var_melt$variable)
Method_var_melt$variable<-ifelse(Method_var_melt$variable=="Prop_slopes","Exclusion of steep slopes",Method_var_melt$variable)

theme_set(theme_bw(base_size=12))
Plot1<-ggplot(Method_var_melt,aes(x=variable,y=value,fill=Method))+geom_bar(stat = "identity",position="dodge")+ theme(axis.text.x = element_text(angle = 90,hjust =1,vjust = 0.5))
Plot2<-Plot1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))                                                                                                                                          
Plot2+scale_fill_brewer(palette="Set1")+xlab("Technique")+ylab("Percentage of sites \nwhich used method")
ggsave("Figures/Method_variation.png",dpi=400,height=4,width=6,units="in")


#now plot histogram of the number of techniques used 

ggplot(Methods,aes(x=All2))+geom_histogram(binwidth=1)+facet_wrap(~Method,scale="free_y")
