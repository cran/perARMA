ppfcoeffab <-
function(ppf,nsamp,printflg,datastr){

     T=nrow(ppf)
     nc=ncol(ppf)                        
     maxn=nc-1                          
     nout= floor((T+1)/2) 
     ppf[is.nan(ppf)]=0

      pkab=phth2ab(ppf) 
      pkab=pkab$a                
      pkab=as.matrix(pkab)

      r=nrow(pkab)
      c=ncol(pkab)
                                   
      pkpv=matrix(1,r,c)
      colsamp=matlab::sum(nsamp,FALSE)                          
      sigma0=1/sqrt(colsamp)                          
      pkpv[1,]=2*(1-pnorm(abs(pkab[1,]),0,sigma0))   
                                               
     sigma=sqrt(2)*sigma0 
       for (j in 2:T)
           {if (j==T & matlab::rem(T,2)==0)
              {sigma=sigma0}
               pkpv[j,]=2*(1-pnorm(abs(pkab[j,]),0,sigma))
           }
      if (printflg)
         {for (k in 1:nc)
            {cat(paste('pi coeffs in ab form for',datastr, 'n= ',k-1,'\n'))
             cat(paste('k  pihat_k   pv','\n'))
          for (kk in 1:T)
           { cat(paste('',kk-1,pkab[kk,k],pkpv[kk,k],'\n'))}
        }    
    }
      
      result = list(pkab=pkab,pkpv=pkpv)  
      class(result) = "ppfcoeffab"
      result
}

