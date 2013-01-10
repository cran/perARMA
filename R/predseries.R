predseries <-function(real,x,T,p,...)
{
 predseries_full<-function(real,x,T,p,predcol,realcol){

 k=length(real)/T

pm<-permest(x,T, 0.05, NaN,'series', pp=0)
pmean=pm$pmean
xd=pm$xd
n=length(xd)
pred<-predictperYW(xd,T,p,NaN, k, predcol, realcol)
new=pred$new

forecast<-matrix(0,T*k,1)

for(l in 1:k){
for(i in 1:T){
forecast[(l-1)*T +i]=new[(l-1)*T+i]+pmean[i]
}}


plot (forecast, xlab="time",ylab= "values of the series", type="l", lwd=1, lty=1,col=predcol)
lines(real, type="l", lwd=1, lty=1,col=realcol)     
title(main="Prediction of the series", sub= paste("No. periods to predict =", k))
legend("bottomright", c("real series","predicted series"), fill=c(realcol,predcol),ncol=2,title="legend")
}

L<-modifyList(list(predcol="red",realcol="blue"), list(real=real,x = x, T=T, p=p,...))

 do.call(predseries_full,L)

 }