#https://data.worldbank.org/indicator/SP.POP.TOTL
require(tidyverse)
DATA<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/population.csv")
RECENT<-DATA[,c("Country.Name","X2019")]
RECENT
countries_deyle<-read.csv("../../Data/deyle_edited_wide_POP.csv")
countries_deyle<-unique(countries_deyle$country)

populations_deyle<-as.data.frame(matrix(nrow=length(countries_deyle),ncol=2))
populations_deyle[,1]<-countries_deyle
for (i in 1:nrow(populations_deyle)){
  if (populations_deyle[i,1] %in% RECENT$Country.Name ){
    #value<-RECENT[which(RECENT$Country.Name==as.character(populations_deyle[i,1])),"X2019"]
    populations_deyle[i,2]<-RECENT[which(RECENT$Country.Name==as.character(populations_deyle[i,1])),"X2019"]
  }
}
write.csv(populations_deyle,"../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/populations_sel.csv")


#https://www.kaggle.com/paultimothymooney/latitude-and-longitude-for-every-country-and-state/data
latlong<-read.csv("../../Data/latlong/552239_1006003_bundle_archive/world_country_and_usa_states_latitude_and_longitude_values.csv")
latlong<-latlong[,c(1:4)]

countries_deyle<-read.csv("../../Data/deyle_edited_wide.csv")
countries_deyle<-unique(countries_deyle$country)
latlongs_deyle<-as.data.frame(matrix(nrow=length(countries_deyle), ncol=3))
latlongs_deyle[,1]<-countries_deyle
for (i in 1:nrow(latlongs_deyle)){
  if (latlongs_deyle[i,1] %in% latlong$country){
    #value<-RECENT[which(RECENT$Country.Name==as.character(latlongs_deyle[i,1])),"X2019"]
    latlongs_deyle[i,2]<-latlong[which(latlong$country==as.character(latlongs_deyle[i,1])),"latitude"]
    latlongs_deyle[i,3]<-latlong[which(latlong$country==as.character(latlongs_deyle[i,1])),"longitude"]
    
  }
}
latlongs_deyle
write.csv(latlongs_deyle,"../../Data/latlong/latlong_sel.csv")
###shorter
latlongs_deyle
latlongs_deyle<-latlongs_deyle[which(latlongs_deyle$V1!="French Guiana"),]

latlongs_deyle <- data.frame(lapply(latlongs_deyle, function(x) {gsub("Egypt", "Egypt, Arab Rep.", x)}))
latlongs_deyle <- data.frame(lapply(latlongs_deyle, function(x) {gsub("Kyrgyzstan", "Kyrgyz Republic", x)}))
latlongs_deyle <- data.frame(lapply(latlongs_deyle, function(x) {gsub("Slovakia", "Slovak Republic", x)}))
write.csv(latlongs_deyle,"../../Data/latlong/latlong_sel_short.csv")
