#finding b
require("minpack.lm")
shrunk_data<-as.data.frame(matrix(nrow=8,ncol=6))
colnames(shrunk_data)<-c("Temp","Hum","Half","sd_Half","b","sd_b")
shrunk_data$Temp<-c(24,24,24,24,28,35,35,35)
shrunk_data$Hum<-c(20,40,60,80,40,20,40,60)
shrunk_data$Half<-c(15.33,11.52,9.15 ,8.33 ,6.11,7.33 ,7.52,  2.26)/24
shrunk_data$sd_Half<-c(2.75,1.72,3.39,1.80,3.02,1.33,1.22,1.42)/24
shrunk_data$b<-log(2)/shrunk_data$Half
shrunk_data$st_b<-log(2)/shrunk_data$sd_Half
plot(shrunk_data$Temp,shrunk_data$b)



plot(shrunk_data$Hum,shrunk_data$b)


Temps<-seq(-10,40,length.out = 1000)
Temp_Model<-nlsLM(b~b0*exp(Temp*g),data=shrunk_data,start = list(g=0,b0=0))

g<-summary(Temp_Model)$coefficients[1,1]
b0<-summary(Temp_Model)$coefficients[2,1]

plot(Temps,b0*exp(g*Temps),"s")
points(shrunk_data$Temp,shrunk_data$b)
shrunk_data_repeats

#result of all this is that equation for how b responds to temp is 
b0
0.256
g
0.078

#Finding q

shrunk_data<-as.data.frame(matrix(nrow=8,ncol=6))
colnames(shrunk_data)<-c("Temp","Hum","Half","sd_Half","b","sd_b")
shrunk_data$Temp<-c(24,24,24,24,28,35,35,35)
shrunk_data$Hum<-c(20,40,60,80,40,20,40,60)
shrunk_data$Half<-c(15.33,11.52,9.15 ,8.33 ,6.11,7.33 ,7.52,  2.26)/24
shrunk_data$sd_Half<-c(2.75,1.72,3.39,1.80,3.02,1.33,1.22,1.42)/24
shrunk_data$b<-log(1/2)/shrunk_data$Half
shrunk_data$st_b<-log(1/2)/shrunk_data$sd_Half
plot(shrunk_data$Temp,shrunk_data$b)



plot(shrunk_data$Hum,shrunk_data$b)


Temps<-seq(-10,40,length.out = 1000)
Temp_Model<-nlsLM(b~b0*exp(Temp*g),data=shrunk_data,start = list(g=0,b0=0))

g<-summary(Temp_Model)$coefficients[1,1]
b0<-summary(Temp_Model)$coefficients[2,1]

plot(Temps,b0*exp(g*Temps),"s")
points(shrunk_data$Temp,shrunk_data$b)
shrunk_data_repeats

#result of all this is that equation for how b responds to temp is 
b0
#-0.256
g
#0.078



