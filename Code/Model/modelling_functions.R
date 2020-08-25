#copied from notebook 
require("deSolve")
require("ggplot2")


beta <-function(c_r,q,d,h,epsilon){
  return(c_r * epsilon *  (( (-1 / q) * (1 - exp( q * d)))/((- 1 / q) * (1 - exp( q * d)) + h)))
}

cr_normal <- function(Climate, m , c, stype) {
  #parameters
  #m is the cooredinates where the contact rate is at a maximum 
  #c is 95 % confidence interval for duration
  #t range_C is the T values the function is calculated for
  #s type determines whether s is calculated using the upper or lower confidence interval
  c_u = c[2]
  c_l = c[1]
  m_C = m[1]
  m_cr = m[2]
  #find s
  if (stype==1){
    s = (c_u - m_C) / 1.96
  }else{
    s = -(c_l - m_C) / 1.96
  }
  #find y axis of normal distribution given s   
  scaling = sqrt(2 * pi) * s * m_cr
  
  return(scaling*(1 / (sqrt(2 * pi) * s))*exp((-(Climate - m_C) ^ 2) / (2 * s ^ 2)))
}


#this is where quantile_95 is used 
cr_climate <- function(Max_Coordinates_cr, range_C, Climate) {
  #finds the equation fiven temp_at_max and range_C. then finds S for given T
  #uses where temp at max relates to range_C to determine stype (to determine whether lower or upper value of range_C is used to find SI)
  if (abs(Max_Coordinates_cr[1]-range_C[1])>abs(Max_Coordinates_cr[1]-range_C[2])){
    stype=0
  }else{
    stype=1
  }
  #finds d at T values for given parameters (contact finds equation using duration_normal_scaled)
  #t_range_C not strictly needed 
  #return(mean_duration(T = T,m = c(temp_at_max, max_mean_contact), c = quantile_95, T_range_C = range_C, stype=stype))
  return(cr_normal(Climate = Climate,m = Max_Coordinates_cr, c = range_C, stype=stype))
  
}

#this is where quantile_95 is used 
cr_climate <- function(Max_Coordinates_cr, range_C, Climate) {

  #finds the equation fiven temp_at_max and range_C. then finds S for given T
  #uses where temp at max relates to range_C to determine stype (to determine whether lower or upper value of range_C is used to find SI)
      #parameters
    #m is the cooredinates where the contact rate is at a maximum 
    #c is 95 % confidence interval for duration
    #t range_C is the T values the function is calculated for
    #s type determines whether s is calculated using the upper or lower confidence interval
    #range_C=c(-22.774,28.865)
    #find s  
  # c_u = range_C[2]
  # c_l = range_C[1]
  # m_C = Max_Coordinates_cr[1]
  # m_cr = Max_Coordinates_cr[2] 
  # 
  # if (abs(Max_Coordinates_cr[1]-range_C[1])>abs(Max_Coordinates_cr[1]-range_C[2])){
  #   s = (-(c_l - m_C) / 1.96)
  # }else{
  #   s = ((c_u - m_C) / 1.96)
  # }
  # #finds d at T values for given parameters (contact finds equation using duration_normal_scaled)
  # #t_range_C not strictly needed 
  # #return(mean_duration(T = T,m = c(temp_at_max, max_mean_contact), c = quantile_95, T_range_C = range_C, stype=stype))
  # return(sqrt(2 * pi) * s * m_cr * (1 / (sqrt(2 * pi) * s))*exp((-(Climate - m_C) ^ 2) / (2 * s ^ 2)))

  if (abs(Max_Coordinates_cr[1]-range_C[1])>abs(Max_Coordinates_cr[1]-range_C[2])){
    s = (-(range_C[1] - Max_Coordinates_cr[1]) / 1.96)
  }else{
    s = ((range_C[2] - Max_Coordinates_cr[1]) / 1.96)
  }
  #finds d at T values for given parameters (contact finds equation using duration_normal_scaled)
  #t_range_C not strictly needed 
  #return(mean_duration(T = T,m = c(temp_at_max, max_mean_contact), c = quantile_95, T_range_C = range_C, stype=stype))
  return(sqrt(2 * pi) * s * Max_Coordinates_cr[2] * (1 / (sqrt(2 * pi) * s))*exp((-(Climate - Max_Coordinates_cr[1]) ^ 2) / (2 * s ^ 2)))
}

cr_climate_quick <- function(Max_Coordinates_cr, s, Climate) {

  
 # if (abs(Max_Coordinates_cr[1]-range_C[1])>abs(Max_Coordinates_cr[1]-range_C[2])){
#    s = (-(range_C[1] - Max_Coordinates_cr[1]) / 1.96)
#  }else{
 #   s = ((range_C[2] - Max_Coordinates_cr[1]) / 1.96)
  #}
  #finds d at T values for given parameters (contact finds equation using duration_normal_scaled)
  #t_range_C not strictly needed 
  #return(mean_duration(T = T,m = c(temp_at_max, max_mean_contact), c = quantile_95, T_range_C = range_C, stype=stype))
  return(sqrt(2 * pi) * s * Max_Coordinates_cr[2] * (1 / (sqrt(2 * pi) * s))*exp((-(Climate - Max_Coordinates_cr[1]) ^ 2) / (2 * s ^ 2)))
}



q_climate <- function(q0,g, Climate) {
  return(q0*exp(g*Climate))
}



Climate_Time_Function<- function(time,min,max,time_at_peak) {
  # here, climate is varying with time at a yearly cycle.
  return((max-min)/2 * cos((2 * pi / 365)* (time - time_at_peak ) ) + (max+min)/2)
}


#SEIR_model <- function(time, values, parms) {
  #function where c changes with the climate and humidity at each timepoint
  #initial values
#  S <- values[1]
#  E <- values[2]
 # I <- values[3]
 # R <- values[4]
#  N <- values[5]
  #parameters
  
  #for (parameter in names(parms)){
  #  assign(parameter,parms[[parameter]])
  #}
  
  #for (variable in names(Climate_Variables)){
  #  assign(variable,Climate_Variables[[variable]])
  #}
 # mu=parms[["mu"]]
#  sigma=parms[["sigma"]]
 # p=parms[["p"]]
#  gamma=parms[["gamma"]]
 # f=parms[["f"]]
  #N=parms[["N"]]
#  nu=parms[["nu"]]
#  h=parms[["h"]]
  #epsilon=parms[["epsilon"]]
  #d=parms[["d"]]
##  Max_cr=parms[["Max_cr"]]
#  #climate_label=parms[["mu"]]
#  g=parms[["g"]]
#  q0=parms[["q0"]]
 # time_at_peak=parms[["Climate_Variables"]][["time_at_peak"]]
#  #range_C=parms[["Climate_Variables"]][["range_C"]]
 # Max_Climate_cr=parms[["Climate_Variables"]][["Max_Climate_cr"]]
  
  
  
#  #find climate given time
#  Climate <- Climate_Time_Function(time=time ,min=range_C[1],max=range_C[2],time_at_peak = time_at_peak)   
#  
#  cr_value=cr_climate(Climate =  Climate,Max_Coordinates_cr = c(Max_Climate_cr,Max_cr), range_C=range_C)
#  q_value= q_climate(q0=q0,g=g, Climate=Climate)
#  beta_value = beta(c_r = cr_value,q =q_value ,d = d,h =h ,epsilon =epsilon )    
#  #seir model
#  dS = nu * N - beta_value * I * S / N - mu * S + f * R
#  dE = beta_value * S * I / N - (sigma + mu) * E
#  dI = sigma * E - (mu + gamma) * I * (1/(1-p))
#  dR = gamma * I - mu * R - f * R
#  
#  dN = nu * N - (p/(1-p)) *(mu + gamma) * I - mu * (S + E + I + R)
  
#  
#  list(c(dS, dE, dI, dR, dN))
#}

SEIR_model <- function(time, values, parms) {
  #function where c changes with the climate and humidity at each timepoint
  #initial values
  S <- values[1]
  E <- values[2]
  I <- values[3]
  R <- values[4]
  #N <- values[5]
  #parameters
  
  #for (parameter in names(parms)){
  #  assign(parameter,parms[[parameter]])
  #}
  
  #for (variable in names(Climate_Variables)){
  #  assign(variable,Climate_Variables[[variable]])
  #}
  mu=parms[["mu"]]
  sigma=parms[["sigma"]]
  p=parms[["p"]]
  gamma=parms[["gamma"]]
  f=parms[["f"]]
  #N=parms[["N"]]
  nu=parms[["nu"]]
  h=parms[["h"]]
  epsilon=parms[["epsilon"]]
  d=parms[["d"]]
  Max_cr=parms[["Max_cr"]]
  #climate_label=parms[["mu"]]
  g=parms[["g"]]
  q0=parms[["q0"]]
  time_at_peak=parms[["Climate_Variables"]][["time_at_peak"]]
  range_C=parms[["Climate_Variables"]][["range_C"]]
  Max_Climate_cr=parms[["Climate_Variables"]][["Max_Climate_cr"]]
  
  s=parms[["Climate_Variables"]][["s"]]
  
  #find climate given time
  Climate <- Climate_Time_Function(time=time ,min=range_C[1],max=range_C[2],time_at_peak = time_at_peak)   
  
  #cr_value=cr_climate(Climate =  Climate,Max_Coordinates_cr = c(Max_Climate_cr,Max_cr), range_C=range_C)
  cr_value=cr_climate_quick(Climate =  Climate,Max_Coordinates_cr = c(Max_Climate_cr,Max_cr), s=s)
  
    q_value= q_climate(q0=q0,g=g, Climate=Climate)
  beta_value = beta(c_r = cr_value,q =q_value ,d = d,h =h ,epsilon =epsilon )    
  currentpop<-S+E+I+R
  #seir model
  infectives<-beta_value * I * S / (currentpop)
    dS = nu * (currentpop) - infectives - mu * S + f * R
    dE = infectives - (sigma + mu) * E
    dI = sigma * E - (mu + gamma) * I * (1/(1-p))
    dR = gamma * I - mu * R - f * R
  #dS = nu * (currentpop) - beta_value * I * S / (currentpop) - mu * S + f * R
  #dE = beta_value * S * I / (currentpop) - (sigma + mu) * E
  #dI = sigma * E - (mu + gamma) * I * (1/(1-p))
  #dR = gamma * I - mu * R - f * R
  
  #dN = nu * N - (p/(1-p)) *(mu + gamma) * I - mu * (S + E + I + R)
  
  
  list(c(dS, dE, dI, dR))
}

SEIR_model_quicker <- function(time, values, parms) {
  #function where c changes with the climate and humidity at each timepoint
  #initial values
  
  
  #find climate given time
  Climate <- Climate_Time_Function(time=time ,min=parms[["Climate_Variables"]][["range_C"]][1],max=parms[["Climate_Variables"]][["range_C"]][2],time_at_peak = parms[["Climate_Variables"]][["time_at_peak"]])   
  
  cr_value=cr_climate_quick(Climate =  Climate,Max_Coordinates_cr = c(parms[["Climate_Variables"]][["Max_Climate_cr"]],parms[["Max_cr"]]), s=parms[["Climate_Variables"]][["s"]])
  
  q_value= q_climate(q0=parms[["q0"]],g=parms[["g"]], Climate=Climate)
  beta_value = beta(c_r = cr_value,q =q_value ,d = parms[["d"]],h =parms[["h"]] ,epsilon =parms[["epsilon"]] )    
  currentpop<-values[1]+values[2]+values[3]+values[4]
  #seir model
  infectives<-beta_value * values[3] * values[1] / (currentpop)
  dS = parms[["nu"]] * (currentpop) - infectives - parms[["mu"]] * values[1] + parms[["f"]] * values[4]
  dE = infectives - (parms[["sigma"]] + parms[["mu"]]) * values[2]
  dI = parms[["sigma"]] * values[2] - (parms[["mu"]] + parms[["gamma"]]) * values[3] * (1/(1-parms[["p"]]))
  dR = parms[["gamma"]] * values[3] - parms[["mu"]] * values[4] - parms[["f"]] * values[4]
  
  
  list(c(dS, dE, dI, dR))
}


#plottime <- function(out) {
  #plotting seir components over time
 # #choose <- parms[["choose"]]
#  plot(x = out[, "time"], y= out[, "S"], ylab = "Number", xlab = "Time",
#       type = "l", xlim=c(min(time), max(time)),ylim=c(0,  max(out[, "N"],na.rm = T)))
#  lines(x = out[, "time"],y= out[, "E"], col = "orange")
#  lines(x = out[, "time"],y= out[, "I"], col = "red",lwd=5)
#  lines(x = out[, "time"],y= out[, "R"], col = "green")
#  lines(x = out[, "time"],y= out[, "N"], col = "purple")  
#  abline(h = max(out[,"I"]),lty=3,col="red")
#  legend( "topright",legend = c("S", "E", "I", "R","N","Max I"), 
#          col = c("black", "orange", "red", "green","purple","red"), lty = c(1,1,1,1,1,3), cex = 0.8)
#  abline(h=seq(min(out[,"time"]),max(out[,"time"],by=365),col="cyan"))
#  maxpercent=max(( out[, "I"]/ out[, "N"])*100)
 # maximum= paste("Maximum I is:",round(maxpercent,digits = 2),"%")
#  text(max(time)/2, max(out[, "N"]),maximum)
  
#}

plottime <- function(out) {
  #plotting seir components over time
  #choose <- parms[["choose"]]
  plot(x = out[, "time"], y= out[, "S"], ylab = "Number", xlab = "Time",
       type = "l", xlim=c(min(time), max(time)),ylim=c(0,  max((out[, "S"]+out[, "E"]+out[, "I"]+out[, "R"]),na.rm = T)))
  lines(x = out[, "time"],y= out[, "E"], col = "orange")
  lines(x = out[, "time"],y= out[, "I"], col = "red",lwd=5)
  lines(x = out[, "time"],y= out[, "R"], col = "green")
  #lines(x = out[, "time"],y= out[, "N"], col = "purple")  
  abline(h = max(out[,"I"]),lty=3,col="red")
  abline(v=seq(min(out[,"time"]),max(out[,"time"]),by=365),col="cyan")
  legend( "topright",legend = c("S", "E", "I", "R","Max I"), 
          col = c("black", "orange", "red", "green","red"), lty = c(1,1,1,1,1,3), cex = 0.8)
  maxpercent=max(( out[, "I"]/ (out[, "S"]+out[, "E"]+out[, "I"]+out[, "R"]))*100)
  maximum= paste("Maximum I is:",round(maxpercent,digits = 2),"%")
  text(max(time)/2, parms[["N"]],maximum)
  
}


find_R0_function<-function(Climate,parms, Climate_Variables_Temp, max_R0_Req){
  #finds r0 given temperatures where peak occurs. can use this to find maximum if max_r0_req is true
  #take variables from list to individual variables
  for (parameter in names(parms)){
    assign(parameter,parms[[parameter]])
  }
  
  for (variable in names(Climate_Variables_Temp)){
    assign(variable,Climate_Variables_Temp[[variable]])
  }
  # max_R0_req- find all or just max?
  #d and cr as functions of time
  cr_value=cr_climate(Climate =  Climate,Max_Coordinates_cr = c(Max_Climate_cr,Max_cr), range_C=range_C)
  q_value= q_climate(q0=q0,g=g, Climate=Climate)
  beta_value = beta(c_r = cr_value,q =q_value ,d = d,h =h ,epsilon =epsilon ) 
  
  
  #below finds the R0 at maximum if required
  #Climate[which.max((sigma / (sigma + mu)) * beta_value/(mu + alpha + gamma))]
  if (max_R0_Req =="True"){
    return(max((sigma / (sigma + mu)) * beta_value / ((mu + gamma)*( 1 / ( 1 - p)))))  
    
  }else{
    return(c((sigma / (sigma + mu)) * beta_value / ((mu + gamma)*( 1 / ( 1 - p )))))  
  }
}


std <- function(x) sd(x)/sqrt(length(x))

