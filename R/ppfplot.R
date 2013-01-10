ppfplot <-
function(ppf,nsamp,alpha,datastr){

     T=nrow(ppf)
     nc=ncol(ppf)
     thr=qnorm(1-alpha/2,0,1/sqrt(nsamp))               
     bfalpha=alpha/(T*nc)                             
     bfthr=qnorm(1-bfalpha/2,0,1/sqrt(nsamp))
       
      ylab.name=expression( pi(t,n+1))
      matplot(seq(0,nc-1),t(ppf), xlab="n = no. samples ", ylab=ylab.name,type="l",lwd=1, ylim=c(-1,1))
     
      lines(seq(0,nc-1),-thr%*%matrix(1,1,nc),type="l", col="red",lwd=1)
      lines(seq(0,nc-1),thr%*%matrix(1,1,nc),type="l",  col="red",lwd=1)
      lines(seq(0,nc-1),-bfthr%*%matrix(1,1,nc),type="l",  col="black",lwd=1)
      lines(seq(0,nc-1),bfthr%*%matrix(1,1,nc),type="l", col="black",lwd=1)

        ccs<-matrix(0,1,T)
      
        for (i in 1:T)  {
                        ccs[i]=i}
                        ccs=as.character(ccs)

       par(xpd=NA,oma=c(3,0,0,0))
        legend(nc-2,1.5,ccs, col=ccs, lty=1,title="t = lags")
      title(paste('perpacf for ', datastr))

       colmeans<-matrix(0,1,nc)
       for (i in 1:nc)  { colmeans[i]=mean(ppf[,i]) }

      nrows=nc
      cat(paste('p-values for pi(t,n+1)=0 for t=0,...,T-1','\n'))
      cat(paste(' n  pi(:,n+1)=0   pi(:,n+1:nmax)=0, for nmax=',nc,'\n'))
    
        pv1save<-0
      for (j in 1:nc)
        {  mean_T=colmeans[j]
           mzscore1=mean_T%*%sqrt(nsamp)%*%sqrt(T)            
           pv1=2*(1-pnorm(abs(mzscore1),0,1))
           pv1save=cbind(pv1save,pv1)
           mzscore2=mean(colmeans[j:nc])%*%sqrt(nsamp)%*%sqrt(T)%*%sqrt(nrows)
           pv2=2*(1-pnorm(abs(mzscore2),0,1))
           if (j==1) {pv2save=pv2}
           nrows=nrows-1
           cat(paste(j-1,pv1,pv2,'\n'))}

       pv1save=pv1save[2:length(pv1save)]
       cat(paste('test for all pi(t,n+1)=0 (1st row, 2nd col) pv=',pv2save,'\n'))
       cat(paste('min of first col, bf corrected for',nc,' rows:',nc*min(pv1save),'\n'))
 }

