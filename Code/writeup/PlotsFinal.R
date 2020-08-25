source("../Model/modelling_functions.R"
)
require("tidyverse")
require("ggplot2")


which_best<-function(country){
  return(data_wider_means[which(data_wider_means$country==country),"meanflu"])
}  

which_function<-function(country,bests){
  return(bests$best[which(bests$country==country)])
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

require("cowplot")

library(gridExtra)
library(grid)

data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ_POP.csv")
data_wider_means<-read.csv("../../Data/data_wider_means_POP.csv")

latlong=read.csv("../../Data/latlong/latlong_sel_short.csv")

plot_influenza<-function(parms,sims){
  correlation_df<-read.csv(paste0("../../Results/fromfunction/cors/",sims,parms[["climate_label"]],"correlation_dataframe.csv"))
  correlation_df$Region<-cut(correlation_df$lat,breaks=c(min(correlation_df$lat)-1,-23.5,23.5,max(correlation_df$lat)+1), labels=c("Southern","Tropics","Northern"))
  correlation_df$Region_Combined<-cut(abs(correlation_df$lat),breaks=c(0,23.5,max(abs(correlation_df$lat)+1)), labels=c("Tropics","Temperate"))
  correlation_df$combination<-as.factor(correlation_df$combination)
  correlation_df$mismatch<-as.factor(correlation_df$mismatch)
  len_mismatch=length(unique(correlation_df$mismatch))
  lensims=len_mismatch*sims
  print(parms[["climate_label"]])
  if (parms[["climate_label"]]=="Temperature"){
    correlation_df$meanCl<-rep(data_wider_means_summ$meanT,each=lensims)
  }
  if (parms[["climate_label"]]=="AH"){
    correlation_df$meanCl<-rep(data_wider_means_summ$meanAH,each=lensims)
  }
  if (parms[["climate_label"]]=="RH"){
    correlation_df$meanCl<-rep(data_wider_means_summ$meanRH,each=lensims)
  }
  #correlation_df_means<-as_tibble(correlation_df) %>% group_by(mismatch) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
  correlation_df_means_country<-as_tibble(correlation_df) %>% group_by(mismatch,lat,country,maxs,mins,time_max,pop,Region_Combined,Region,meanCl) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
  #to get weekly values 
    
  model_means<-read.csv(paste0("../../Results/fromfunction/",sims,parms[["climate_label"]],".csv"))
  model_means_combined<-model_means %>% group_by(country,mismatch,week) %>%  summarise(meanI=mean(meanI),meanR0=mean(meanR0),meanCl=mean(meanCl),.groups="keep")
  model_means$lat<-latlong[model_means$country,"V2"]
  model_means$country<-data_wider_means_summ[model_means$country,"country"]
  model_means$Region<-cut(model_means$lat,breaks=c(min(model_means$lat)-1,-23.5,23.5,max(model_means$lat)+1), labels=c("Southern","Tropics","Northern"))
  model_means$Region_Combined<-cut(abs(model_means$lat),breaks=c(0,23.5,max(abs(model_means$lat)+1)), labels=c("Tropics","Temperate"))
  model_means_summary<-model_means %>% group_by(country,Region_Combined,mismatch) %>% summarise(meanmeanR0=mean(meanR0),meanmeanI=mean(meanI),maxmeanR0=max(meanR0),maxmeanI=max(meanI),meanq=mean(meanq),.groups="keep")

    #deee <- model_means %>% group_by(country,Region_Combined) %>% summarise(ranges=max(meanq)-min(meanq))
    #deee <- model_means %>% group_by(country,Region_Combined) %>% summarise(ranges=max(meanCl)-min(meanCl))
    ##deee <- model_means %>% group_by(country,Region_Combined) %>% summarise(ranges=max(meancr)-min(meancr))
    #plot(deee$Region_Combined,deee$ranges)
  #range comparisons

  m<-model_means_summary %>%  group_by(Region_Combined,mismatch) %>% summarise(mI=mean(meanmeanI),mR=mean(meanmeanR0),maI=mean(maxmeanI),maR=mean(maxmeanR0), smI=std(meanmeanI),smR=std(meanmeanR0),smaI=std(maxmeanI),smaR=std(maxmeanR0),.groups="keep")
  #correlation_df_means<-as_tibble(correlation_df) %>% group_by(mismatch) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
  pdf(paste0("../../Writeup/",sims,parms[["climate_label"]],"rangelat.pdf"))
  print(ggplot(data=correlation_df_means_country, aes(x= abs(lat) ,y=maxs-mins))+geom_point()+xlab("Absolute Value of Latitude")+theme_bw()+ylab(paste0(parms[["climate_label_long"]]," Range"))+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15))) 
  graphics.off()
  
  pdf(paste0("../../Writeup/boxandwhisker_facet",sims,parms[["climate_label"]],".pdf"),width= 10,height= 5)
  print(ggplot(data=correlation_df, aes(x= as.factor(mismatch),y= corsI)) +geom_boxplot()  +theme_bw()+xlab("Mismatch") +ylab("Mean Correlation")+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)) +facet_wrap(~Region_Combined))
  graphics.off() 
  
  
  maxR<-max(m$maR)
  maxI<-max(m$maI)
  
  
  
  if (parms[["climate_label"]]=="Temperature"){
  
  bests<-correlation_df %>% group_by(country,lat,maxs,mins,Region,Region_Combined,meanCl,combination) %>% summarise(best=mismatch[which.max(corsI)],.groups="keep")
  bests<-bests %>% group_by(country,Region_Combined) %>% summarise(best=getmode(best))
  write.csv(bests,paste0("../../Results/bestmismatchtemp",sims,".csv"))
  plot(bests$Region_Combined,bests$best)  
  tab<-table(bests$best,bests$Region_Combined)
  write.csv(tab/matrix(rep(colSums(tab),5),nrow=5,byrow = T),"../../Writeup/mismatch.csv")
  pdf(paste0("../../Writeup/mismatchseverity_meanmeanI",parms[["climate_label"]],sims,parms[["climate_label"]],".pdf"),width=7.5,height= 5)
  print(ggplot(data=m,aes(mismatch,mI,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Mean Proportion Infected")+xlab("Mismatch")+
          geom_errorbar(aes(ymin=mI-smI,ymax=mI+smI))+scale_y_continuous()+scale_x_continuous(limits = c(-0.125,1.125), breaks = seq(0,1,by=1/(len_mismatch-1)))+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15))) +labs(colour="Region")
  graphics.off()
  pdf(paste0("../../Writeup/mismatchseverity_meanmeanro",parms[["climate_label"]],sims,parms[["climate_label"]],".pdf"),width=7.5,height= 5)
  print(ggplot(data=m,aes(mismatch,mR,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Mean R0 Infected")+xlab("Mismatch")+
          geom_errorbar(aes(ymin=mR-smR,ymax=mR+smR))+scale_y_continuous()+scale_x_continuous(limits = c(-0.125,1.125), breaks = seq(0,1,by=1/(len_mismatch-1)))+labs(colour="Region")+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15))) 
  graphics.off()
  pdf(paste0("../../Writeup/mismatchseverity_maxmeanI",parms[["climate_label"]],sims,parms[["climate_label"]],".pdf"),width=7.5,height= 5)
  print(ggplot(data=m,aes(mismatch,maI,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Maximum Proportion Infected")+xlab("Mismatch")+
          geom_errorbar(aes(ymin=maI-smaI,ymax=maI+smaI))+scale_y_continuous()+scale_x_continuous(limits = c(-0.125,1.125), breaks = seq(0,1,by=1/(len_mismatch-1)))+labs(colour="Region")+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15))) 
  graphics.off()
  
  pdf(paste0("../../Writeup/mismatchseverity_maxmeanro",parms[["climate_label"]],sims,parms[["climate_label"]],".pdf"),width=7.5,height= 5)
  print(ggplot(data=m,aes(mismatch,maR,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Maximum R0 Infected")+xlab("Mismatch")+
          geom_errorbar(aes(ymin=maR-smaR,ymax=maR+smaR))+scale_y_continuous()+scale_x_continuous(limits = c(-0.125,1.125), breaks = seq(0,1,by=1/(len_mismatch-1)))+labs(colour="Region")+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15))) 
  graphics.off()

  #   #for some reason mean AH is a good predictor of mismatch????????    

  
  
    #could the reason be incorreact asumptions? 
    #get AH and contact rates for each country
    
    #ggplot(model_means_2,aes(meanCl,meanR0,group=country))+geom_line() # 
    #ggplot(model_means_2,aes(meanCl,meancr,group=country,col=best))+geom_line() # 
    #ggplot(model_means_2,aes(meanCl,meanq,group=country,col=best))+geom_line() # 
  if(sims==1){
    ranges_comparison<-correlation_df[which(correlation_df$mismatch==0),]
    ranges_comparison<-ranges_comparison %>% group_by(Region_Combined) %>% summarise(ranges=maxs-mins,.groups="keep")  
    ranges_comparison_min<-ranges_comparison %>% group_by(Region_Combined) %>% summarise(means=mean(ranges),sterr=std(ranges),.groups="keep")
    
    write.csv(ranges_comparison_min,"../../Writeup/tempranges.csv")
    ranges_comparison<- model_means%>% group_by(country,Region_Combined) %>% summarise(meanq=mean(meanq),.groups="keep") 
    ranges_comparison_min<- ranges_comparison %>% group_by(Region_Combined) %>% summarise(means=mean(meanq),sterr=std(meanq),.groups="keep")
    write.csv(ranges_comparison_min,"../../Writeup/survivalmeans.csv")
    
    #ranges_comparison_min<-ranges_comparison %>% group_by(Region_Combined) %>% summarise(means=mean(ranges),sterr=std(ranges),.groups="keep")
    ranges_comparison<- model_means %>% group_by(country,Region_Combined) %>% summarise(diff=max(meanq)-min(meanq),.groups="keep") 
    ranges_comparison_min<- ranges_comparison %>% group_by(Region_Combined) %>% summarise(means=mean(diff),sterr=std(diff),.groups="keep")
  
    write.csv(ranges_comparison_min,"../../Writeup/survivalranges.csv")
    
    bests<-correlation_df %>% group_by(country,lat,maxs,mins,Region,Region_Combined,meanCl) %>% summarise(best=mismatch[which.max(corsI)],.groups="keep")
    #bests<-correlation_df %>% group_by(country,lat,maxs,mins,Region) %>% summarise(best=which.max(corsI),.groups="keep")
    model_means$best<-sapply(model_means$country,which_function,bests=bests)
   model_means$mismatch<-as.factor(model_means$mismatch)
  #model_means_2<-model_means[which(model_means$mismatch==1),] 
    model_means_2<-model_means[which(model_means$mismatch==model_means$best),] 
    model_means_2$data=NA
    for (coun in unique(model_means_2$country)){
      data=which_best(coun)
      data_space=which(model_means_2$country==coun)
      data_space=data_space[1:length(data)]
      model_means_2[data_space,"data"]<-data
    }
    cbPalette <- c( "#000000","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")
    
    #make a function that finds how much to shift 
    # do for seperate then pic best
    #then go back to hump thing- look for temp= nothing
    #for (country in unique(model_means_2$country)){
    pdf(paste0("../../Writeup/",parms[["climate_label"]],"sel_I_TIME_DATA.pdf"),width=15,height=5)
    
    model_means_3<-model_means_2[which(model_means_2$mismatch==1),]
    # if (parms[["climate_label"]]=="AH"){
    #   model_means_3<-model_means_2[which(model_means_2$country %in% c("El Salvador", "France", "Costa Rica",  "Panama","Paraguay")),] 
    #   model_means_3$c_order = factor(model_means_3$country, levels=c( "France","Paraguay","El Salvador","Costa Rica", "Panama"))  
    #   mult<-50
      
    model_means_3<-model_means_2[which(model_means_2$country %in% c("El Salvador", "France", "Dominican Republic",  "Senegal","Thailand")),] 
    model_means_3$c_order = factor(model_means_3$country, levels=c( "France", "Senegal","Thailand", "Dominican Republic","El Salvador"))  
    model_means_3$mismatch_cntry = paste0(model_means_3$mismatch,": ",model_means_3$country)
    mult=50
  
    p<-ggplot(data=model_means_3, aes(week, meanI,group=country,col=country))+geom_line()+facet_wrap(~mismatch_cntry,nrow=1)
    p<-ggplot(data=model_means_3, aes(week, meanI,colour="Model"))+geom_line()+facet_wrap(~mismatch_cntry,nrow=1)
    
    p <- p + geom_line(aes(y = mult*data, colour = "Data"))
    
    p <- p + scale_y_continuous(sec.axis = sec_axis(~./mult, name = "per Capita Flu Cases from Data"))

    p <- p + labs(x = "Week",
                 y = "Proportion of Infected Individuals in Model", colour="") +theme_bw()
    #p <- p + labs(x = "Week",
    #              y = "Model Number of Infectious") +theme_bw()
    
    p<-p+  scale_colour_manual(values=cbPalette)
    p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)) 
    
    print(p)
    graphics.off()
    # for (country in unique(model_means_2$country)){
    #   model_means_4<-model_means_2[which(model_means_2$country==country),]
    #   png(paste0("../../Writeup/loop",country,parms[["climate_label"]],"sel_I_TIME_DATA.png"),width=50*10,height=50*5)
    #   mult<-max(model_means_4$meanI)/max(model_means_4$data,na.rm = T)
    #   
    #   p<-ggplot(data=model_means_4, aes(week, meanI,group=mismatch,col=mismatch))+geom_line()
    #   p <- p + geom_line(aes(y = mult*data, colour = "Data"))
    #   p <- p + scale_y_continuous(sec.axis = sec_axis(~./mult, name = "Data Flu Cases"))
    #   p <- p + labs(x = "Week",
    #                 y = "Model Number of Infectious", colour="Mismatch") +theme_bw()
    #   p<-p+  scale_colour_manual(values=cbPalette)
    #   p<-p+theme(text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)) 
    #   print(p)
    #   graphics.off()
    # }
 

  }
  
  }
}  
# 



parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
              g=0.085,q0=-9.079,Climate_Variables=NA,climate_label_long="Temperature")
sims=1

plot_influenza(parms,sims=1)

plot_influenza(parms,sims=100)

#covid_plots(parms)
parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="RH",
              g=0.085,q0=-9.079,Climate_Variables=NA,climate_label_long="Relative Humidity")

plot_influenza(parms,sims=1)

plot_influenza(parms,sims=100)
#sims=1
parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="AH",extra="Absolute Humidity",
              g=0.062,q0=-30.162,Climate_Variables=NA)
plot_influenza(parms,sims=1)
plot_influenza(parms,sims=100)#need t oredo rh with proportion


pdf(paste0("../../Writeup/AHtemp.pdf"))
ggplot(data_wider_means,aes(meantemp,meanAH ))+geom_point()+theme_bw()+ylab("Absolute Humidity")+xlab("Temperature")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)) 
graphics.off()#covid_plots(parms)

pdf(paste0("../../Writeup/RHtemp.pdf"))
ggplot(data_wider_means,aes(meantemp,meanRH ))+geom_point()+theme_bw()+ylab("Relative Humidity")+xlab("Temperature")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)) 
graphics.off()#covid_plots(parms)


parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
              g=0.085,q0=-9.079,Climate_Variables=NA,climate_label_long="Temperature")
sims=1
IT<-data_wider_means[data_wider_means$country=="Italy",]

C<-seq(min(IT$meantemp),max(IT$meantemp),length.out = 1000)
dummy_dataframe<-as.data.frame(matrix(ncol=5,nrow=length(C)*5))
colnames(dummy_dataframe)<-c("Temperature","Contact_Rate","Survival_Rate","Max_cr","Mismatch")
dummy_dataframe$Climate<-rep(C,5)
range<-seq(min(C),max(C),length.out = 5)
dummy_dataframe$Max_cr<-rep(range,each=length(C))
dummy_dataframe$Mismatch<-(dummy_dataframe$Max_cr-min(C))/(max(C)-min(C))
for (k in unique(dummy_dataframe$Max_cr)){
  dummy_dataframe[which(dummy_dataframe$Max_cr==k),"Contact_Rate"]<-cr_climate(Max_Coordinates_cr = c(k,30),range_C=c(min(C),max(C)),Climate =C)
}
dummy_dataframe$Survival_Rate<-q_climate(q0 = parms[["q0"]],g=parms[["g"]],Climate = dummy_dataframe$Climate)

add<-30-max(dummy_dataframe$Survival_Rate,na.rm = T)
pdf(paste0("../../Writeup/",parms[["climate_label_long"]],"mismatch_vis.pdf"),width=15,height=5)

p<-ggplot(data=dummy_dataframe, aes(Climate,Contact_Rate,col="Contact Rate"))+geom_line()+
  geom_line(aes(y = Survival_Rate+add, colour = "Growth Rate"))+scale_x_continuous(breaks = c(5,10,15,20,25))+
  scale_y_continuous(breaks=c(0,10,20,30),sec.axis = sec_axis(~.-add, name = "Virus Growth Rate (1/day)"))+ ylab("Human Contact Rate (1/day)")+ labs(colour="")+
  scale_colour_manual(values = c("blue", "red"))+facet_wrap(~Mismatch,nrow=1)+xlab("Temperature")+
  theme_bw()    +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15))

print(p)
graphics.off()

pdf(paste0("../../Writeup/",parms[["climate_label_long"]],"mismatch_vis_mini.pdf"),width=6,height=5)
dummy_dataframe_mini<-dummy_dataframe[which(dummy_dataframe$Mismatch %in% c(0,1)),]
p<-ggplot(data=dummy_dataframe_mini, aes(Climate,Contact_Rate,col="Contact Rate"))+geom_line()+
  geom_line(aes(y = Survival_Rate+add, colour = "Growth Rate"))+scale_x_continuous(breaks = c(5,10,15,20,25))+
  scale_y_continuous(breaks=c(0,10,20,30),sec.axis = sec_axis(~.-add, name = "Virus Growth Rate (1/day)"))+ ylab("Human Contact Rate (1/day)")+ labs(colour="")+
scale_colour_manual(values = c("blue", "red"))+facet_wrap(~Mismatch,nrow=1)+xlab("Temperature")+
  theme_bw()    +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15))

print(p)
graphics.off()




#COVID

parms = list( climate_label="Temperature",Climate_Variables=NA,climate_label_long="Temperature")
sims=1
 covid_means<-read.csv(paste0("../../Results/fromfunction/covid/1temperature_shift.csv"))
 covid_means$lat<-latlong[covid_means$country,"V2"]
 covid_means$country<-data_wider_means_summ[covid_means$country,"country"]
 covid_means$Region<-cut(covid_means$lat,breaks=c(min(covid_means$lat)-1,-23.5,23.5,max(covid_means$lat)+1), labels=c("Southern","Tropics","Northern"))
# shift_vec<-c(0,0.25,0.5,0.75,1)
# covid_means$shift<-shift_vec[covid_means$shift+1]
 covid_means$Region_Combined<-cut(abs(covid_means$lat),breaks=c(0,23.5,max(abs(covid_means$lat)+1)), labels=c("Tropics","Temperate"))
 #covid_means$mismatch<-as.factor(covid_means$mismatch)
 bests<-read.csv("../../Results/bestmismatchtemp100.csv")
 #rownames(bests)=bests[,2]
 covid_means$best_influenza<-sapply(covid_means$country,which_function,bests)
 covid_means$shift<-covid_means$mismatch-covid_means$best_influenza
# 
# 
 covid_means_summary<-covid_means %>% group_by(country,Region_Combined,combination,shift) %>% summarise(meanmeanR0=mean(meanR0),meanmeanI=mean(meanI),maxmeanR0=max(meanR0),maxmeanI=max(meanI),.groups="keep")
# 
# 
# 
 m<-covid_means_summary %>%  group_by(Region_Combined,shift) %>% summarise(mI=mean(meanmeanI),mR=mean(meanmeanR0),maI=mean(maxmeanI),maR=mean(maxmeanR0), smI=std(meanmeanI),smR=std(meanmeanR0),smaI=std(maxmeanI),smaR=std(maxmeanR0))
# correlation_df_means<-as_tibble(covid_means) %>% group_by(mismatch) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
 maxR<-max(m$maR)
 maxI<-max(m$maI)
 
 pdf(paste0("../../Writeup/shiftseverity_meanmeanI",sims,".pdf"),width=7.5,height= 5)
 print(ggplot(data=m,aes(shift,mI,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Mean Proportion Infected")+xlab("Shift")+
         geom_errorbar(aes(ymin=mI-smI,ymax=mI+smI),width=0.25)+scale_y_continuous()+scale_x_continuous()+labs(colour="Region")+
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)))
 graphics.off()
 pdf(paste0("../../Writeup/shiftseverity_meanmeanro",sims,".pdf"),width=7.5,height= 5)
 print(ggplot(data=m,aes(shift,mR,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Mean R0")+xlab("Shift")+
         geom_errorbar(aes(ymin=mR-smR,ymax=mR+smR))+scale_y_continuous()+scale_x_continuous()+labs(colour="Region")+
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)))
 graphics.off()
 pdf(paste0("../../Writeup/shiftseverity_maxmeanI",sims,".pdf"),width=7.5,height= 5)
 print(ggplot(data=m,aes(shift,maI,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Maximum Proportion Infected")+xlab("Shift")+
         geom_errorbar(aes(ymin=maI-smaI,ymax=maI+smaI))+scale_y_continuous()+scale_x_continuous()+labs(colour="Region")+
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)))
 graphics.off()
 
 pdf(paste0("../../Writeup/shiftseverity_maxmeanro",sims,".pdf"),width=7.5,height= 5)
 print(ggplot(data=m,aes(shift,maR,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Maximum R0")+xlab("Shift")+
         geom_errorbar(aes(ymin=maR-smaR,ymax=maR+smaR))+scale_y_continuous()+scale_x_continuous()+labs(colour="Region")+
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)))
 graphics.off()
# 
 leg<-ggplot(data=m,aes(shift,maR,group=Region_Combined, colour = Region_Combined))+geom_point()+geom_line()+theme_bw()+ylab("Maximum R0")+xlab("Shift")+
   geom_errorbar(aes(ymin=maR-smaR,ymax=maR+smaR))+scale_y_continuous()+scale_x_continuous()+labs(colour="Region")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)) 
 legend <- get_legend(leg)
 pdf(paste0("../../Writeup/legend.pdf"),width=6,height= 1)
 grid.newpage()
 grid.draw(legend) 
 graphics.off()
 
 trop<-covid_means[which(covid_means$Region=="Southern"),]
trop<-trop[which(trop$shift==0),]
# 
 noshift<-covid_means[which(covid_means$shift==0),]
 noshift$Country<-noshift$country
 selected_countries<-noshift[which(noshift$country %in% c("United Kingdom", "Cambodia","India","Australia","Niger")),]
 pdf(paste0("../../Writeup/selectedseriesR0",".pdf"),width=10,height=5)
 print(ggplot(data=selected_countries,aes(week,meanR0,group=Country,color=Region))+geom_line()+theme_bw()+ylab("R0")+xlab("Week")+
         theme(legend.position="bottom",text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15))+facet_wrap(c("Country"),nrow=1)+
         theme(strip.background = element_blank(), strip.text.x = element_blank())+
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)))
 graphics.off()
 
# 
# 
 covid_means_summ<-covid_means %>% group_by(Region,country,shift,combination) %>% summarise(peakweek=which.max(meanR0),peakI=which.max(meanI))
 covid_means_summ<-covid_means_summ[which(covid_means_summ$shift==0),]
 pdf(paste0("../../Writeup/covidtime_R0_peekweek_",sims,".pdf"),width=10,height=5)
 print(ggplot(data=covid_means_summ,aes(y=peakweek,x=Region))+geom_boxplot()+theme_bw()+ylab("Week of Maximal R0")+xlab("Region")+
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)))
 graphics.off()
 pdf(paste0("../../Writeup/covidtime_I_peekweek_",sims,".pdf"))
 print(ggplot(data=covid_means_summ,aes(y=peakI,x=Region))+geom_boxplot()+theme_bw()+ylab("Week of Maximal I")+xlab("Region")+
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)))
 graphics.off()
 
 summs<-covid_means_summ %>% group_by(Region) %>% summarise(meanweek=mean(peakweek),errweek=std(peakweek),.groups="keep")
 write.csv(summs,paste0("../../Writeup/covidpeekweek",sims,".csv"))

#Variance
 sims=1
 parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
               N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
               g=0.085,q0=-9.079,Climate_Variables=NA,climate_label_long="Temperature")
correlation_df<-read.csv("../../Results/fromfunction/cors/var1Temperaturecorrelation_dataframe.csv")
correlation_df$Region<-cut(correlation_df$lat,breaks=c(min(correlation_df$lat)-1,-23.5,23.5,max(correlation_df$lat)+1), labels=c("Southern","Tropics","Northern"))
correlation_df$Region_Combined<-cut(abs(correlation_df$lat),breaks=c(0,23.5,max(abs(correlation_df$lat)+1)), labels=c("Tropics","Temperate"))
correlation_df$combination<-as.factor(correlation_df$combination)
correlation_df$mismatch<-as.factor(correlation_df$mismatch)
correlation_df$v_scal<-as.factor(correlation_df$v_scal)
correlation_df_TR<-correlation_df[which(correlation_df$Region_Combined=="Tropics"),]
correlation_df_TE<-correlation_df[which(correlation_df$Region_Combined=="Temperate"),]

pdf(paste0("../../Writeup/boxandwhisker_facet_TR",sims,parms[["climate_label"]],".pdf"),width= 10,height= 5)
print(ggplot(data=correlation_df_TR, aes(x= as.factor(mismatch),y= corsI)) +geom_boxplot()  +theme_bw()+xlab("Mismatch") +ylab("Mean Correlation")+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)) +facet_wrap(~v_scal))
graphics.off() 

pdf(paste0("../../Writeup/boxandwhisker_facet_TE",sims,parms[["climate_label"]],".pdf"),width= 10,height= 5)
print(ggplot(data=correlation_df_TE, aes(x= as.factor(mismatch),y= corsI)) +geom_boxplot()  +theme_bw()+xlab("Mismatch") +ylab("Mean Correlation")+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 15),axis.text = element_text(size=15),axis.title = element_text(size=15)) +facet_wrap(~v_scal))
graphics.off() 


