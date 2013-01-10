predictperYW <-function(x,T,p,missval, k,...)
{
   predictperYW_full <- function(x, T, p, missval, k, realcol, 
        predcol) {
        numpts = length(x)
        nper = floor(length(x)/T)
        estimators <- perYW(x, T, p, missval) 	
        phi = estimators$phi			
        phi = as.matrix(phi)
        del = estimators$del						
        del = as.matrix(del)
        estimators = cbind(phi, del)		
        exestimators <- rbind(estimators, estimators) 	
							       # for the external 2 periods 
        for (i in 1:(nper + k - 2)) 		 
		{
            	exestimators <- rbind(exestimators, estimators)
        	}
        for (m in 1:k)				               #periods of predictions re end 
		{
            	xdpred <- matrix(0, T, 1) 	
            	for (j in 1:T)			               #j is time in base period
			{
                	for (i in 1:p) 		               # i is the lag of the X used for prediction 
				{
                  		xdpred[j] = xdpred[j] + exestimators[T * (nper + 
                    		m - 1) + j, i] * x[T * (nper + m - 2) + j - 
                    		i]
				                               # standing at time j in the current period m;
				                               # exestimators[T * (nper + m - 1) + j, i] coeff for lag i
				                               # multiply X at time [T * (nper + m - 2) + j - i]
				                               # which is ith time back  
                		}
            		}
            		x = c(x, xdpred)	
        	}
	
        length = length(x)
        plot(seq(1, numpts), x[1:numpts], xlab = "time", ylab = "Prediction of the series", 
            type = "l", lwd = 1, lty = 1, col = realcol)
        lines(seq((numpts), length), x[(numpts):length], xlim = c((numpts + 
            1), length), type = "l", lwd = 2, lty = 1, col = predcol)
        title(main = "Prediction of series after removing periodic mean", 
            sub = paste("No. periods to predict =", k))
        legend("bottomright", c(expression(real), expression(prediction)), 
            fill = c(realcol, predcol), ncol = 2, title = "legend")
        new = x[(length - k * T + 1):length]
        xp = x
        result = list(xp = xp, new = new)
        class(result) = "predictperYW"
        result
    }
    L <- modifyList(list(realcol = "blue", predcol = "red"), 
        list(x = x, T = T, p = p, missval=missval, k = k, ...))
    do.call(predictperYW_full, L)
}