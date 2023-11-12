#' Subset AR Graph for "Selectmodel" Object
#' 
#' A graphical depiction is given of the output from `SelectModel`.
#' 
#' @usage plot(x, ...) ## S3 method for class 'Selectmodel'
#' @param x out from `SelectModel`.
#' @param ... optional arguments.
#' @details
#' The relative plausibility of Model A vs the best Model B, is defined as 
#' \eqn{R=e^{(AIC_B-AIC_A)/2}}. Values of R less than 1 R is defined similarly 
#' if the BIC/UBIC criterion is used.
#' @returns No value. Plot produced as side-effect.
#' @author A.I. McLeod.
#' @seealso [SelectModel()].
#' @examples
#' # takes about 10 seconds
#' ## Not run: 
#' out<-SelectModel(log(Willamette),lag.max=150,ARModel="AR",Best=5,Criterion="AIC")
#' plot(out)
#' ## End(Not run)
#' 
#' @export
plot.Selectmodel <-
  function(x, ...){
    ans<-x #for clarity we rename object
    par(mar = c(5,4,4,5)+0.1)
    on.exit(par(mar=c(5, 4, 4, 2) + 0.1)) #default
    if (is.list(ans)){    
      X<-Y<-a<-numeric(0)
      for (i in 1:length(ans)){
        p<-ans[[i]]$p
        X<-c(X,p)
        aic<-ans[[i]][[2]]
        Y<-c(Y,rep(aic,length(p)))
        a<-c(a,aic)
      }
      ylab<-names(ans[[1]])[2]
      plot(X,Y,xlab="lags", ylab=ylab,pch=15, cex=2)
      abline(h=a, lty=2)
      title(main=paste(attr(ans,"model"),"model selection"))
      ytic<-a
      r<-round(exp((min(Y)-ytic)/2),3)
      axis(4,at=ytic, labels=r)
      mtext(side=4, line=2.7, text="Relative Plausibility")
    }
    if (is.matrix(ans)){
      X<-ans[,1]
      Y<-ans[,2]
      ylab<-dimnames(ans)[[2]][2]
      plot(X,Y,ylab=ylab,pch=15,cex=2, xlab="p")
      title(main="AR(p) model selection")
      ytic<-Y
      r<-round(exp((min(Y)-ytic)/2),3)
      axis(4,at=Y,labels=r)
      mtext(side=4, line=2.7, text="Relative Plausibility")
    }
    invisible()
  }
