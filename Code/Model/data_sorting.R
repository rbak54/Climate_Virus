
require("tidyverse")


###############################################################
#initial data tidying- obtain and widen data
data_influenza<-read.csv("../../Data/deyle_data.txt")
#data are already per capita
#get week number  
data_influenza_order<-dplyr::group_by(data_influenza,variable,country,year) %>% mutate(week=row_number())
#keep only full years
data_influenza_order_max<- data_influenza_order %>% group_by(variable,country,year) %>% mutate(no=max(week))
data_influenza_full_years<-data_influenza_order_max[which(data_influenza_order_max$no>51),]
#widen
data_edited_wide<-data_influenza_full_years %>% spread(variable,value)
#save
for (i in unique(data_edited_wide$country)){
  png(paste0("../../Results/Plots/climate/",i,".png"))
      toplot<-data_edited_wide[which(data_edited_wide$country==i),]
      plot(toplot$RH,toplot$T)
      graphics.off()
}
#remove iraq because extremely high RH values- likely to be errors
DATA<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/population.csv")
RECENT<-DATA[,c("Country.Name","X2019")] #or midway??

data_edited_wide <-data_edited_wide[which(data_edited_wide$country!="Iraq"),]

write.csv(data_edited_wide,"../../Data/deyle_edited_wide.csv")
#countries missing:
#Egypt
#French Guiana
#Kyrgyzstan
#Slovakia

# so change recent
#Slovak Republic
#Egypt, Arab Rep.
#Kyrgyz Republic
countries_deyle<-data_edited_wide
#RH extremely high - eror likely


countries_deyle <- data.frame(lapply(countries_deyle, function(x) {gsub("Egypt", "Egypt, Arab Rep.", x)}))
countries_deyle <- data.frame(lapply(countries_deyle, function(x) {gsub("Kyrgyzstan", "Kyrgyz Republic", x)}))
countries_deyle <- data.frame(lapply(countries_deyle, function(x) {gsub("Slovakia", "Slovak Republic", x)}))
countries_deyle <-countries_deyle[which(countries_deyle$country!="French Guiana"),]
data_edited_wide<-countries_deyle
write.csv(data_edited_wide,"../../Data/deyle_edited_wide_POP.csv")
#H added here

###############################################################

###############################################################
#open data
data_wider<-read.csv("../../Data/deyle_edited_wide.csv")
################################################################

################################################################
#data tidying and reformatting 
#remove missing flu and temp values
data_wider<-data_wider[which(data_wider$T!="NAN"),]
data_wider<-data_wider[which(data_wider$flu!="NAN"),]
data_wider<-data_wider[which(data_wider$RH!="NAN"),]
data_wider<-data_wider[which(data_wider$RH!="NaN"),]
data_wider<-data_wider[which(data_wider$T!="NaN"),]
data_wider<-data_wider[which(data_wider$flu!="NaN"),]

data_wider<-data_wider[which(data_wider$AH!="NAN"),]
data_wider<-data_wider[which(data_wider$AH!="NaN"),]
#convert to celcius
data_wider$T<-(data_wider$T-32)*(5/9)
#find weekly mean flu and temp over the years
data_wider_means<-data_wider %>% group_by(country,week) %>% summarise(meanflu=mean(flu),meantemp=mean(T),.groups="keep")
data_wider_means<-data_wider %>% group_by(country,week) %>% summarise(meanflu=mean(flu),meantemp=mean(T), meanRH=mean(RH),meanAH=mean(AH),.groups="keep")

#convert week to day
data_wider_means$day<-(data_wider_means$week-1)*7+3.5

# find summary for each country- week and values for maxima and mimima
data_wider_means_summ<-data_wider_means %>% group_by(country) %>% summarise(minflu=min(meanflu),maxflu=max(meanflu),peakflu=week[which.max(meanflu)],troughflu=week[which.min(meanflu)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)],.groups="keep")
data_wider_means_summ<-data_wider_means %>% group_by(country) %>% summarise(minflu=min(meanflu),maxflu=max(meanflu),peakflu=week[which.max(meanflu)],troughflu=week[which.min(meanflu)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)],minRH=min(meanRH),maxRH=max(meanRH),peakRH=week[which.max(meanRH)], troughRH=week[which.min(meanRH)],.groups="keep")
#,minAH=min(meanAH),maxAH=max(meanAH),peakAH=week[which.max(meanAH)], troughAH=week[which.min(meanAH)]
# convert to data frame
data_wider_means_summ<-as.data.frame(data_wider_means_summ)
#summary for each year for each country
data_wider_summ<-data_wider %>% group_by(country, year) %>% summarise(minflu=min(flu),maxflu=max(flu),peakflu=week[which.max(flu)],troughflu=week[which.min(flu)],minT=min(T),maxT=max(T),peakT=week[which.max(T)], troughT=week[which.min(T)],.groups="keep")
data_wider_summ<-data_wider %>% group_by(country, year) %>% summarise(minflu=min(flu),maxflu=max(flu),peakflu=week[which.max(flu)],troughflu=week[which.min(flu)],minT=min(T),maxT=max(T),peakT=week[which.max(T)], troughT=week[which.min(T)],minRH=min(RH),maxRH=max(RH),peakRH=week[which.max(RH)],troughRH=week[which.min(RH)],.groups="keep")
#minAH=min(AH),maxAH=max(AH),peakAH=week[which.max(AH)],troughAH=week[which.min(AH)]
write.csv(data_wider_means,"../../Data/data_wider_means.csv")
write.csv(data_wider_means_summ,"../../Data/data_wider_means_summ.csv")
write.csv(data_wider_summ,"../../Data/data_wider_summ.csv")

#################################################################
#summary of pop shortened dataset
data_wider<-read.csv("../../Data/deyle_edited_wide_POP.csv")

data_wider<-data_wider[which(data_wider$T!="NAN"),]
data_wider<-data_wider[which(data_wider$flu!="NAN"),]
data_wider<-data_wider[which(data_wider$T!="NaN"),]
data_wider<-data_wider[which(data_wider$flu!="NaN"),]
data_wider<-data_wider[which(data_wider$RH!="NAN"),]
data_wider<-data_wider[which(data_wider$RH!="NaN"),]
#data_wider<-data_wider[which(data_wider$AH!="NAN"),]
#data_wider<-data_wider[which(data_wider$AH!="NaN"),]
#convert to celcius
data_wider$T<-(data_wider$T-32)*(5/9)
#find weekly mean flu and temp over the years
#data_wider_means<-data_wider %>% group_by(country,week) %>% summarise(meanflu=mean(flu),meantemp=mean(T),.groups="keep")
data_wider_means<-data_wider %>% group_by(country,week) %>% summarise(meanflu=mean(flu),meantemp=mean(T), meanRH=mean(RH), meanAH=mean(AH),varflu=var(flu),vartemp=var(T),.groups="keep")

#convert week to day
data_wider_means$day<-(data_wider_means$week-1)*7+3.5

# find summary for each country- week and values for maxima and mimima
data_wider_means_summ<-data_wider_means %>% group_by(country) %>% summarise(minflu=min(meanflu),maxflu=max(meanflu),peakflu=week[which.max(meanflu)],
                            troughflu=week[which.min(meanflu)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)],.groups="keep")
data_wider_means_summ<-data_wider_means %>% group_by(country) %>% summarise(minflu=min(meanflu),maxflu=max(meanflu),peakflu=week[which.max(meanflu)],
                            troughflu=week[which.min(meanflu)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)],minRH=min(meanRH),
                            maxRH=max(meanRH),peakRH=week[which.max(meanRH)], troughRH=week[which.min(meanRH)],meanvarflu=mean(varflu),meanvartemp=mean(vartemp), varvarflu=var(varflu), varvartemp=var(vartemp), meanmeanflu=mean(meanflu),
                            minAH=min(meanAH),maxAH=max(meanAH),peakAH=week[which.max(meanAH)], troughAH=week[which.min(meanAH)],meanAH=mean(meanAH),meanT=mean(meantemp),meanRH=mean(meanRH),.groups="keep")

# convert to data frame
data_wider_means_summ<-as.data.frame(data_wider_means_summ)
#summary for each year for each country
data_wider_summ<-data_wider %>% group_by(country, year) %>% summarise(minflu=min(flu),maxflu=max(flu),peakflu=week[which.max(flu)],
                                      troughflu=week[which.min(flu)],minT=min(T),maxT=max(T),peakT=week[which.max(T)], troughT=week[which.min(T)],.groups="keep")
data_wider_summ<-data_wider %>% group_by(country, year) %>% summarise(minflu=min(flu),maxflu=max(flu),peakflu=week[which.max(flu)],
              troughflu=week[which.min(flu)],minT=min(T),maxT=max(T),peakT=week[which.max(T)], troughT=week[which.min(T)],minRH=min(RH),maxRH=max(RH),peakRH=week[which.max(RH)],
              troughRH=week[which.min(RH)],
minAH=min(AH),maxAH=max(AH),peakAH=week[which.max(AH)],troughAH=week[which.min(AH)],st.groups="keep")
write.csv(data_wider_means,"../../Data/data_wider_means_POP.csv")
write.csv(data_wider_means_summ,"../../Data/data_wider_means_summ_POP.csv")
write.csv(data_wider_summ,"../../Data/data_wider_summ_POP.csv")

#################################################################
