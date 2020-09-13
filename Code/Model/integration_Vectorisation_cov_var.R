sims_range=c(1,10)

#data and functions 
source("modelling_functions.R")

require("tidyverse")
require("ggplot2")
source("latincube.R")
print("latin cube done")
source("data_sorting.R")
print("data sorting done")
source("gazeteer.R")
print("gazeteer run")

data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ_POP.csv")
data_wider_summ<-read.csv("../../Data/data_wider_summ_POP.csv")
data_wider_means<-read.csv("../../Data/data_wider_means_POP.csv")
#integration  N changing
population_deyle<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/populations_sel.csv")
latlong=read.csv("../../Data/latlong/latlong_sel_short.csv")
pop<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/populations_sel.csv")


ranges_all<-c(min(data_wider_means_summ$minT),max(data_wider_means_summ$minT))

integration_general<-function(parms,sims,time){
  if (sims>1){
    Y<-make_lhs(n=sims,parms = parms)
    write.csv(Y,"../../Data/latincube.csv")
  }else{
    Y<-as.data.frame(matrix(c(parms[["epsilon"]],parms[["h"]],parms[["mu"]],parms[["f"]]),nrow=1))
    colnames(Y)<-c("epsilon","h","mu","f")
    
    print("estimated parameters used")
  }
  
  extra_cols<-12
  model_means=matrix(ncol =13)
  names(model_means)<-c("country","week","mismatch","v_scal","combination","meanI","varI","meanR0","meancr","meanCl","varCl","meanq","varR0")
  
  parms_temp<-parms
  parms_temp[["epsilon"]]<-NA
  parms_temp[["h"]]<-NA
  parms_temp[["mu"]]<-NA
  parms_temp[["f"]]<-NA
  parms_temp[["N"]]<-NA
  
  for (combination in 1:nrow(Y)){
    parms_temp[["epsilon"]]<-Y$epsilon[combination]
    parms_temp[["h"]]<-Y$h[combination]
    parms_temp[["mu"]]<-Y$mu[combination]
    parms_temp[["f"]]<-Y$f[combination]
    for (v_scal in parms[["v_scal_range"]]){
    for (location_index in 1:nrow(data_wider_means_summ)){
      parms_temp[["N"]]=population_deyle[location_index,3]
      start = c(S = (1-1e-4)*parms_temp[["N"]],
                E = 0.00*parms_temp[["N"]], 
                I =(1e-4)*parms_temp[["N"]],
                R = 0*parms_temp[["N"]])
      
      if(parms[["climate_label"]]=="Temperature"){
        peak_contact_seq<-seq(data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"],length.out=5)
        for (i in peak_contact_seq){
          #difference between temperature where contact rate is highest and lowest temperature in range (i.e virus does best survivasl or virus does best contact)
          #mismatch= 0 is when contact rate is highest at low temp 
          parms_temp[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakT"]*7,range_C=c(data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"]),Max_Climate_cr=i)
          mismatch=(i-parms_temp[["Climate_Variables"]][["range_C"]][1])/(parms_temp[["Climate_Variables"]][["range_C"]][2]-parms_temp[["Climate_Variables"]][["range_C"]][1])
          mismatch<-round(mismatch,2)
          if (abs(i-parms_temp[["Climate_Variables"]][["range_C"]][1])>abs(i-parms_temp[["Climate_Variables"]][["range_C"]][2])){
            s = v_scal *(-(parms_temp[["Climate_Variables"]][["range_C"]][1]- i) / 1.96)
          }else{
            s = v_scal * ((parms_temp[["Climate_Variables"]][["range_C"]][2] - i) / 1.96)
          }
          parms_temp[["Climate_Variables"]][["s"]]=s
          temp = ode(    y = start,    time = time,    func = SEIR_model_quicker,    parms = parms_temp)
          temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"],data_wider_means_summ[location_index,"peakT"],mismatch,NA,NA,NA,NA,NA,combination,v_scal),nrow(temp)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
          colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","T","R0","cr","q","week","combination","v_scal")
          temp_extra[,"T"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms_temp[["Climate_Variables"]][["range_C"]][1],max=parms_temp[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms_temp[["Climate_Variables"]][["time_at_peak"]] )
          temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"T"],parms=parms_temp, Climate_Variables_Temp=parms_temp[["Climate_Variables"]], max_R0_Req=F)
          temp_extra[,"cr"]<-cr_climate_quick(Max_Coordinates_cr = c(parms_temp[["Climate_Variables"]][["Max_Climate_cr"]],parms_temp[["Max_cr"]]),Climate=temp_extra[,"T"],s=s)
          
          temp_extra[,"q"]<-q_climate(parms_temp[["q0"]],parms_temp[["g"]],temp_extra[,"T"])
          
          year<-ceiling((temp[,"time"])/365)
          day<-(temp[,"time"]-(year-1)*365)
          temp_extra[,"week"]<-ceiling(day/7)
          temp_extra[,"country"]<-location_index
          temp<-cbind(temp,temp_extra)
          temp<-as_tibble(temp)
          temp<-temp[which(year>1),]
          model_means_temp<- temp %>% group_by(country,week,mismatch,combination,v_scal) %>% summarise(meanI=mean(I/(S+E+I+R)),varI=var(I/(S+E+I+R)),meanR0=mean(R0),meancr=mean(cr),
                                                                                                meanCl=mean(T),varCl=var(T),meanq=mean(q),varR0=var(R0),.groups="keep")
          model_means<-bind_rows(model_means,model_means_temp)
        }
      }
    
      
      if(parms[["climate_label"]]=="RH"){
        peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=5)
        for (i in peak_contact_seq){
          #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
          #mismatch= 0 is when contact rate is highest at low RH 
          parms_temp[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
          mismatch=(i-parms_temp[["Climate_Variables"]][["range_C"]][1])/(parms_temp[["Climate_Variables"]][["range_C"]][2]-parms_temp[["Climate_Variables"]][["range_C"]][1])
          mismatch<-round(mismatch,2)
          if (abs(i-parms_temp[["Climate_Variables"]][["range_C"]][1])>abs(i-parms_temp[["Climate_Variables"]][["range_C"]][2])){
            s = v_scal*(-(parms_temp[["Climate_Variables"]][["range_C"]][1]- i) / 1.96)
          }else{
            s = v_scal*((parms_temp[["Climate_Variables"]][["range_C"]][2] - i) / 1.96)
          }
          parms_temp[["Climate_Variables"]][["s"]]=s
          temp = ode(    y = start,    time = time,    func = SEIR_model_quicker,    parms = parms_temp)
          temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA,NA,NA,NA,combination,v_scal)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
          colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","RH","R0","cr","q","week","combination","v_scal")
          temp_extra[,"RH"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms_temp[["Climate_Variables"]][["range_C"]][1],max=parms_temp[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms_temp[["Climate_Variables"]][["time_at_peak"]] )
          temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"RH"],parms=parms_temp, Climate_Variables_Temp=parms_temp[["Climate_Variables"]], max_R0_Req=F)
          temp_extra[,"cr"]<-cr_climate_quick(Max_Coordinates_cr = c(parms_temp[["Climate_Variables"]][["Max_Climate_cr"]],parms_temp[["Max_cr"]]),Climate=temp_extra[,"RH"],s=s)
          
          temp_extra[,"q"]<-q_climate(parms_temp[["q0"]],parms_temp[["g"]],temp_extra[,"RH"])
          
          year<-ceiling((temp[,"time"])/365)
          day<-(temp[,"time"]-(year-1)*365)
          temp_extra[,"week"]<-ceiling(day/7)
          temp_extra[,"country"]<-location_index
          temp<-cbind(temp,temp_extra)
          temp<-as_tibble(temp)
          temp<-temp[which(year>1),]
          model_means_temp<- temp %>% group_by(country,week,mismatch,combination,v_scal) %>% summarise(meanI=mean(I/(S+E+I+R)),varI=var(I/(S+E+I+R)),meanR0=mean(R0),meancr=mean(cr),
                                                                                                meanCl=mean(RH),varCl=var(RH),meanq=mean(q),varR0=var(R0),.groups="keep")
      
          
          model_means<-bind_rows(model_means,model_means_temp)

        }
      }
      
      
      if(parms[["climate_label"]]=="AH"){
        peak_contact_seq<-seq(data_wider_means_summ[location_index,"minAH"],data_wider_means_summ[location_index,"maxAH"],length.out=5)
        for (i in peak_contact_seq){
          #difference between Absolute_humidity where contact rate is highest and lowest Absolute_humidity in range (i.e virus does best survivasl or virus does best contact)
          #mismatch= 0 is when contact rate is highest at low AH 
          parms_temp[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakAH"]*7,range_C=c(data_wider_means_summ[location_index,"minAH"],data_wider_means_summ[location_index,"maxAH"]),Max_Climate_cr=i)
          mismatch=(i-parms_temp[["Climate_Variables"]][["range_C"]][1])/(parms_temp[["Climate_Variables"]][["range_C"]][2]-parms_temp[["Climate_Variables"]][["range_C"]][1])
          #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
          mismatch<-round(mismatch,2)
          if (abs(i-parms_temp[["Climate_Variables"]][["range_C"]][1])>abs(i-parms_temp[["Climate_Variables"]][["range_C"]][2])){
            s = v_scal*(-(parms_temp[["Climate_Variables"]][["range_C"]][1]- i) / 1.96)
          }else{
            s = v_scal*((parms_temp[["Climate_Variables"]][["range_C"]][2] - i) / 1.96)
          }
          parms_temp[["Climate_Variables"]][["s"]]=s
          temp = ode(    y = start,    time = time,    func = SEIR_model_quicker,    parms = parms_temp)
          #  png(paste0("../../Results/Plots/model_series/AH",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
          #  plottime(temp)
          # graphics.off()
          temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minAH"],data_wider_means_summ[location_index,"maxAH"],data_wider_means_summ[location_index,"peakAH"],mismatch,NA,NA,NA,NA,NA,combination,v_scal)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
          colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","AH","R0","cr","q","week","combination","v_scal")
          temp_extra[,"AH"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms_temp[["Climate_Variables"]][["range_C"]][1],max=parms_temp[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms_temp[["Climate_Variables"]][["time_at_peak"]] )
          temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"AH"],parms=parms_temp, Climate_Variables_Temp=parms_temp[["Climate_Variables"]], max_R0_Req=F)
          temp_extra[,"cr"]<-cr_climate_quick(Max_Coordinates_cr = c(parms_temp[["Climate_Variables"]][["Max_Climate_cr"]],parms_temp[["Max_cr"]]),Climate=temp_extra[,"AH"],s=s)
          temp_extra[,"q"]<-q_climate(parms_temp[["q0"]],parms_temp[["g"]],temp_extra[,"AH"])
          
          year<-ceiling((temp[,"time"])/365)
          day<-(temp[,"time"]-(year-1)*365)
          temp_extra[,"week"]<-ceiling(day/7)
          temp_extra[,"country"]<-location_index
          temp<-cbind(temp,temp_extra)
          temp<-as_tibble(temp)
          temp<-temp[which(year>1),]
          
          model_means_temp<- temp %>% group_by(country,week,mismatch,combination,v_scal) %>% summarise(meanI=mean(I/(S+E+I+R)),varI=var(I/(S+E+I+R)),meanR0=mean(R0),meancr=mean(cr),
                                                                                                meanCl=mean(AH),varCl=var(AH),meanq=mean(q),varR0=var(R0),.groups="keep")
        
          model_means<-bind_rows(model_means,model_means_temp)

          
        }
        }
      }
    }
  }
  model_means<-model_means[-1,]
  return(model_means)
}


correlations<-function(mean,country){
  data_sub<-data_wider_means[which(data_wider_means$country==country),"meanflu"]
  return(unname(cor.test(mean[1:52],data_sub[1:52])[4][[1]]))
  
}

correlation_function<-function(model_means,parms){
  model_means$country<-data_wider_means_summ[model_means$country,"country"]
  correlation_df<- as_tibble(model_means)  %>% group_by(country,combination,mismatch,v_scal) 
  #  correlation_df<-na.omit(correlation_df$)
  if (parms[["climate_label"]]=="Temperature"){
    correlation_df <- correlation_df %>%  summarise(corsI=correlations(meanI,unique(country)),corsR=correlations(meanR0,unique(country)),maxs=max(meanCl),mins=min(meanCl),time_max=which.max(meanCl),means=mean(meanCl),.groups="keep")
  }
  if (parms[["climate_label"]]=="RH"){
    correlation_df <- correlation_df %>%  summarise(corsI=correlations(meanI,unique(country)),corsR=correlations(meanR0,unique(country)),maxs=max(meanCl),mins=min(meanCl),time_max=which.max(meanCl),means=mean(meanCl),.groups="keep")
  }
  if (parms[["climate_label"]]=="AH"){
    correlation_df <- correlation_df %>%  summarise(corsI=correlations(meanI,unique(country)),corsR=correlations(meanR0,unique(country)),maxs=max(meanCl),mins=min(meanCl),time_max=which.max(meanCl),means=mean(meanCl),.groups="keep")
  }

  #addidng latitude and pop info
  
  matrix_extra<-as.data.frame(matrix(nrow=nrow(correlation_df),ncol=3))
  colnames(matrix_extra)<-c("lat","long","pop")
  for (i in 1:nrow(matrix_extra)){
    matrix_extra[i,c(1:2)]<-latlong[which(latlong$V1==correlation_df$country[i]),c(3:4)]
    matrix_extra[i,3]<-pop[which(pop$V1==correlation_df$country[i]),"V2"]
    
  }
  correlation_df<-as.data.frame(correlation_df)
  correlation_df<-cbind(correlation_df,matrix_extra)
  return(correlation_df)
}

run_integration<-function(parms,sims_range){
  for (sss in sims_range){
    
    sims=sss
    time = seq(1,365*11, by=1)
    extra<-parms[["extra"]]
    model_means<-integration_general(parms,sims,time)
    write.csv(model_means,paste0("../../Results/fromfunction/",sims,parms[["climate_label"]],extra,".csv"))
    correlation_df<-correlation_function(model_means,parms)
    write.csv(correlation_df,paste0("../../Results/fromfunction/cors/",sims,parms[["climate_label"]],"correlation_dataframe",extra,".csv"))
    
    bests<-correlation_df %>% group_by(country) %>% summarise(best=mismatch[which.max(corsI)],.groups="keep")
    
    correlation_df_means<-as_tibble(correlation_df) %>% group_by(mismatch) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
    correlation_df_means_country<-as_tibble(correlation_df) %>% group_by(mismatch,lat,maxs,mins,time_max,pop) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
    
  
  }
}

run_integration_covid<-function(parms,sims_range){
  for (sss in sims_range){
    
    sims=sss
    time = seq(1,365*11, by=1)
    extra<-parms[["extra"]]
    model_means<-integration_general(parms,sims,time)
    write.csv(model_means,paste0("../../Results/fromfunction/covid",sims,parms[["climate_label"]],extra,".csv"))
  }
}
run_integration_var<-function(parms,sims_range){
  for (sss in sims_range){
    sims=sss
    time = seq(1,365*11, by=1)
    extra<-parms[["extra"]]
    model_means<-integration_general(parms,sims,time)
    write.csv(model_means,paste0("../../Results/fromfunction/var",sims,parms[["climate_label"]],extra,".csv"))

      correlation_df<-correlation_function(model_means,parms)
    write.csv(correlation_df,paste0("../../Results/fromfunction/cors/var",sims,parms[["climate_label"]],"correlation_dataframe",extra,".csv"))
      }
}


parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=26.97,climate_label="Temperature",extra="",
              g=0.085,q0=-9.079,Climate_Variables=NA,v_scal_range=c(1))
a<-Sys.time()
run_integration(parms,sims_range)
print(Sys.time()-a)
parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=26.97,climate_label="RH",
              g=0.0209,q0=-21.98,Climate_Variables=NA,v_scal_range=c(1)) 
a<-Sys.time()
run_integration(parms,sims_range)
print(Sys.time()-a)

parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=26.97,climate_label="AH",extra="",
              g=0.062,q0=-30.162,Climate_Variables=NA,v_scal_range=c(1))
a<-Sys.time()
run_integration(parms,sims_range)
print(Sys.time()-a)


#COVID

parms = list( mu = 2.06e-5,sigma = 1/5 ,p = 0.01, gamma =1/21,f=0.05,
              N = NA, nu = 5.07e-5, h=0.25/24 ,epsilon= 0.01, d=4/24,Max_cr=26.97,climate_label="Temperature",
              g=0.07841,q0=-0.2557,Climate_Variables=NA,v_scal_range=c(1))
sims_range=c(1)
run_integration_covid(parms,sims_range)
#VAR
sims_range<-c(1)

parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=26.97,climate_label="Temperature",extra="",
              g=0.085,q0=-9.079,Climate_Variables=NA,v_scal_range=c(0.75,1,1.5,2,2.5,3))
a<-Sys.time()
run_integration_var(parms,sims_range)
print(Sys.time()-a)

