
make_lhs<-function(n, parms){
require("lhs")
X<-randomLHS(n,4)

colnames(X)<-c("epsilon","h","mu","f")
Y<-matrix(nrow=nrow(X), ncol=4)
colnames(Y)<-colnames(X)

Y[,1]<-qunif(p=X[,1],max=parms[["epsilon"]]+0.5*parms[["epsilon"]],min=parms[["epsilon"]]-parms[["epsilon"]]*0.5)
Y[,2]<-qunif(p=X[,2],max=parms[["h"]]+0.5*parms[["h"]],min=parms[["h"]]-parms[["h"]]*0.5)
Y[,3]<-qunif(p=X[,3],max=parms[["mu"]]+0.5*parms[["mu"]],min=parms[["mu"]]-parms[["mu"]]*0.5)
Y[,4]<-qunif(p=X[,4],max=parms[["f"]]+0.5*parms[["f"]],min=parms[["f"]]-parms[["f"]]*0.5)


Y<-as.data.frame(Y)

return(Y)
}
