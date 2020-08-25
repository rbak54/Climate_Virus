#Harper

require("minpack.lm")

data<-read.csv("../../Data/converted harper 1961.csv",header =F)
# in this dataset nil becomes 0, trace or n.d becomes NA,and the midpoint of the range is used
#where there is a range I have used the midpoint
data<-as.data.frame(data)


viability<-as.vector(as.matrix((data[4:10,])))
formatted_data<-as.data.frame(matrix(nrow=length(viability),ncol = 9))
formatted_data[,5]<-viability
formatted_data[,4]<-rep(c(1/60/60/24,1/12/24,1/2/24,1/24,4/24,6/24,23/24),length(viability)/7)
formatted_data[,6]<-rep(1:11,each=7)

formatted_data[,1]<-rep(as.vector(as.matrix((data[1,]))),each=7)
formatted_data[,2]<-rep(as.vector(as.matrix((data[2,]))),each=7)
formatted_data[,3]<-rep(as.vector(as.matrix((data[3,]))),each=7)
colnames(formatted_data)<-c("Temp","Hum","Rep","Time","Viability","Experiment","b1","v0","AH")

est0=6.11
lrv= 2256000/461.15
formatted_data$AH<-est0*exp(lrv*((1/273.15)-(1/(formatted_data$Temp+273.15))))*formatted_data$Hum/100


times<-seq(0,1,length.out = 100)
par(mfrow=c(3,4))
for (i in (unique(formatted_data$Experiment))){
 data_subset<-formatted_data[which(formatted_data[,6]==i),]

 model<-nlsLM(Viability~v0*exp(-b*Time),data=data_subset,start=list(v0=100,b=0),lower = c(v0=0  ,b=0),upper = c(v0=100  ,b=5000))

 formatted_data[which(formatted_data[,6]==i),7]<-summary(model)$parameters[2,1]
 formatted_data[which(formatted_data[,6]==i),8]<-summary(model)$parameters[1,1] 
 plot(times,summary(model)$parameters[1,1]*exp(-times*summary(model)$parameters[2,1]),"s",ylim=c(0,100))
 points(data_subset$Time,data_subset$Viability,col="blue",cex=1)
 }

shrunk_data<-formatted_data[seq(1,77,by=7),c(1,2,3,4,5,7,8,9)]


par(mfrow=c(2,3))
plot(shrunk_data$Temp,shrunk_data$b1)
plot(shrunk_data$Temp,shrunk_data$v0)
plot(shrunk_data$Hum,shrunk_data$b1)
plot(shrunk_data$Hum,shrunk_data$v0)
plot(shrunk_data$AH,shrunk_data$b1)
plot(shrunk_data$AH,shrunk_data$v0)


##TEMPERATURE




shrunk_data_repeats<-as.data.frame(matrix(NA,ncol=ncol(shrunk_data),nrow = sum(shrunk_data$Rep)))
colnames(shrunk_data_repeats)=colnames(shrunk_data)
index<-0
for (i in 1:nrow(shrunk_data)) {
 shrunk_data_repeats[(index+1):(index+shrunk_data$Rep[i]) , ] <-shrunk_data[i,]
 index<-index+shrunk_data$Rep[i]
}

Temps<-seq(-10,40,length.out = 1000)
Temp_Model<-nlsLM(b1~b0*exp(Temp*g),data=shrunk_data_repeats,start = list(g=0,b0=0))
g<-summary(Temp_Model)$coefficients[1,1]
b0<-summary(Temp_Model)$coefficients[2,1]
  
plot(Temps,b0*exp(g*Temps),"s")
points(shrunk_data_repeats$Temp,shrunk_data_repeats$b1)

#result of all this is that equation for how b responds to temp is 
#(b0 * exp (g * Temps))
b0
#9.079
g
#0.085
#q=-b so q0=-b0

#HUMIDITY

Hums<-seq(0,100,length.out = 1000)
Hum_Model<-nlsLM(b1~b0*exp(Hum*g),data=shrunk_data_repeats,start = list(g=0,b0=0))
g<-summary(Hum_Model)$coefficients[1,1]
b0<-summary(Hum_Model)$coefficients[2,1]

plot(Hums,b0*exp(g*Hums),"s")
points(shrunk_data_repeats$Hum,shrunk_data_repeats$b1)

#result of all this is that equation for how b responds to temp is 
#(b0 * exp (g * Temps))
b0
#21.98
g
#0.0209
###AH

AH<-seq(min(formatted_data$AH),max(formatted_data$AH),length.out = 1000)
AH_Model<-nlsLM(b1~b0*exp(AH*g),data=shrunk_data_repeats,start = list(g=0,b0=0.000000))

AH_Model<-nlsLM(b1~b0*exp(AH*g),data=shrunk_data_repeats,start = list(g=0.1,b0=0.0000001))

g<-summary(AH_Model)$coefficients[1,1]
b0<-summary(AH_Model)$coefficients[2,1]

g 
#0.06
b0
#30.162

#plot(shrunk_data_repeats$AH,shrunk_data_repeats$q1)
plot(AH,b0*exp(g*AH),"s")
points(shrunk_data_repeats$AH,shrunk_data_repeats$b1)

