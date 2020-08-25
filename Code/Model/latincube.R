
make_lhs<-function(n, parms){
require("lhs")
X<-randomLHS(n,4)
#mu
#epsilon
#f
#h
#eg epsilon between 0.2 and 1
#epsilon
colnames(X)<-c("epsilon","h","mu","f")
Y<-matrix(nrow=nrow(X), ncol=4)
colnames(Y)<-colnames(X)
# parms = list( mu = 2.06e-5,epsilon = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
#                    N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
# #                    g=0.085,q0=-9.079,Climate_Variables=NA)
# Y[,1]<-qunif(X[,1],0.4,1)
# Y[,2]<-qunif(X[,2],01/60/24,6/24)
# Y[,3]<-qunif(X[,3],0.5e-5,3.5e-5)
# Y[,3]<-qunif(X[,3],0.5e-5,2.5e-5)
# Y[,4]<-qunif(X[,4],0.01,0.5)

#when you vary by 50* either way pattern looks the same 
Y[,1]<-qunif(p=X[,1],max=parms[["epsilon"]]+0.5*parms[["epsilon"]],min=parms[["epsilon"]]-parms[["epsilon"]]*0.5)
Y[,2]<-qunif(p=X[,2],max=parms[["h"]]+0.5*parms[["h"]],min=parms[["h"]]-parms[["h"]]*0.5)
Y[,3]<-qunif(p=X[,3],max=parms[["mu"]]+0.5*parms[["mu"]],min=parms[["mu"]]-parms[["mu"]]*0.5)
Y[,4]<-qunif(p=X[,4],max=parms[["f"]]+0.5*parms[["f"]],min=parms[["f"]]-parms[["f"]]*0.5)

# 
# parms_temp = list( mu = 2.06e-5,epsilon = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
#                    N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
#                    g=0.085,q0=-9.079,Climate_Variables=NA)

# head(Y)
Y<-as.data.frame(Y)

# parms_temp<-list( mu = NA,epsilon = NA ,p = 0.001, gamma =0.25,f=NA,
#                            N = 1, nu = 5.07e-5, h=NA ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
#                            g=0.085,q0=-9.079,Climate_Variables=NA)
return(Y)
}
#for (i in 1:nrow(Y)){
#  parms_temp[["epsilon"]]<-Y$epsilon[i]
#  parms_temp[["h"]]<-Y$h[i]
#  parms_temp[["mu"]]<-Y$mu[i]
#  parms_temp[["f"]]<-Y$f[i]
  
  #run odes
  
  #if quicker store out
  #or find means for  each week for each mismatch for each country, then find correlation for each mismatch for each country with data
  
  #matrix with row for country mismatch parms (either simulation no or actual values) and correlation
  #then could do a linear model
