#' Multi-Panel or Single-Panel Time Series Plot with Aspect-Ratio Control
#' 
#' Cleveland (1993) pointed out that the aspect-ratio is important in graphically 
#' showing the rate-of-change or shape information. For many time series, it is 
#' preferably to set this ratio to 0.25 than the default. In general, Cleveland (1993) 
#' shows that the best choice of aspect-ratio is often obtained by if the average 
#' apparent absolute slope in the graph is about 45 deg. But for many stationary 
#' time series, this would result in an aspect-ratio which would be too small. As a 
#' comprise we have chosen a default of 0.25 but the user can select other choices. 
#' 
#' @usage TimeSeriesPlot(z, SubLength = Inf, aspect = 0.25, type="l", 
#'   xlab = "Observation Number", ylab=NULL, main=NULL, ...)
#' @param z ts object or vector, time series data.
#' @param SubLength maximum number of data points per panel. 
#'   Default SubLength=Inf and regular graphics. For trellis graphics, 
#'   set SubLength to a finite value.
#' @param aspect optional setting for the aspect-ratio.
#' @param type plot type, default type="l" join points with lines.
#' @param xlab label for horizontal axis.
#' @param ylab optional label for vertical axis.
#' @param main optional title.
#' @param ... optional arguments passed to `xyplot`.
#' @details
#' If z has attribute "title" containing a character string, this is used on the plot. 
#' Time series input using the function `Readts` always have this attribute set.
#' @returns If SubLength is finite, the lattice package is used and a graphic 
#'   object of class trellis is produced. Otherwise, the standard R graphics 
#'   system is used and the plot is produced as a side-effect and there is no output.
#' @note Requires `lattice` library.   
#' @author A.I. McLeod.
#' @references W.S. Cleveland (1993), Visualizing Data.
#' @seealso [plot.ts()], [Readts()].
#' @examples
#' # from built-in datasets
#' TimeSeriesPlot(AirPassengers)
#' title(main="Monthly number of trans-Atlantic airline passengers")
#' 
#' # compare plots for lynx series
#' plot(lynx)
#' TimeSeriesPlot(lynx, type="o", pch=16, ylab="# pelts", main="Lynx Trappings")
#' 
#' # lattice style plot
#' data(Ninemile)
#' TimeSeriesPlot(Ninemile, SubLength=200)
#' 
#' @export
TimeSeriesPlot <-
  function(z, SubLength=Inf, aspect=0.25, type="l", xlab="Observation Number", ylab=NULL, main=NULL, ...){
    if (SubLength==Inf) {
      pin <- par("pin")
      default.aspect <- pin[2]/pin[1]
      if(aspect < default.aspect)
        pin.new <- c(pin[1], pin[1] * aspect)
      else 
        pin.new <- c(pin[2]/aspect, pin[2])
      par(pin = pin.new)
      yl <- ylab
      ti <- main
      if (is.null(ti))
        ti <- attr(z, "title")
      plot(z, type=type,  main=ti, ylab=yl, xlab = xlab, ...)
      par(pin=pin)
    }
    else { #this part requires lattice library
      n<-length(z)
      nblocks<-ceiling(n/SubLength)
      ti <- main
      if (is.null(ti))
        ti <- attr(z, "title")
      if (nblocks ==1)
        xyplot(z~(1:n) , aspect=aspect, type=type, ylab=ylab, xlab = xlab, main=ti, ...)
      else {
        y<-x<-numeric(0)
        u<-1:n
        if (n>SubLength){
          LastBlockz<-z[(n-SubLength+1):n]
          LastBlockx<-u[(n-SubLength+1):n]
          for (i in 1:(nblocks-1)) {
            ii<-seq((i-1)*SubLength+1,SubLength*i)
            y<-c(y,z[ii])
            x<-c(x,u[ii])
          }
          y<-c(y,LastBlockz)
          x<-c(x,LastBlockx)
        }
        epoch<-ordered(rep(1:nblocks,rep(SubLength,nblocks)))
        xyplot(y~x | epoch, aspect=aspect, type=type, ylab=ylab, xlab=xlab, scales=list(x="sliced",y="same"),strip=FALSE,main=ti, ...)
      }
    }
  }
