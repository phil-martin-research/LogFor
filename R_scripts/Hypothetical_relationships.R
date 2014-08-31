#script to produce code for hypothesised RIL and C logging

Future<-data.frame(Study=rep(seq(1,10,by = 1),each=20),Method=as.factor(rep(c("Conventional","RIL"),each = 10,times = 10)),
                   Basal=abs(rnorm(n = 200,mean = 7,sd = 3)),Ref=abs(rnorm(n = 200,mean = 300,sd = 60)))

Future$Post<-Future$Ref-(15*Future$Basal)+rnorm(200,10,4)
Future$Post<-ifelse(Future$Method=="RIL",Future$Post+(abs(rnorm(200,8,2))*Future$Basal),Future$Post)

Future$LnRR<-log(Future$Post)-log(Future$Ref)

plot(Future$Basal,exp(Future$LnRR)-1)

library(ggplot2)


theme_set(theme_bw(base_size=25))
fut_plot<-ggplot(Future,aes(x=Basal,y=exp(LnRR)-1,colour=Method))+geom_point(size=4,alpha=0.5)+geom_smooth(se=F,method="lm",size=2)
fut_plot2<-fut_plot+ylab("Proportional change in \nbiomass following logging")
fut_plot3<-fut_plot2+scale_size_continuous(range=c(5,10))+geom_hline(y=0,lty=2,size=2)
fut_plot4<-fut_plot3+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))
fut_plot5<-fut_plot4+scale_colour_brewer(palette="Set1")+ guides(colour = guide_legend(override.aes = list(size=4)))
fut_plot5+xlab("Basal area per hectare removed during logging")
setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Chapters/5. Tropical forest degradation/LogFor/Figures")
ggsave("Potential_plots.jpeg",height=7,width=10,dpi=1200)
