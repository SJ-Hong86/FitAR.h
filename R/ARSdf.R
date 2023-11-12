#' Autoregressive Spectral Density Function
#' 
#' Spectral density function of AR(p) is computed.
#' 
#' @usage ARSdf(phi, pFFT = 8)
#' @param phi	 AR Coefficient vector.
#' @param pFFT	 FFT with $2^pFFF$ FFF frequencies, default 8.
#' @details 
#' The Fast Fourier Transform (FFT) is used to compute the spectral density function.
#' @returns A vector of the density function values, \eqn{(f(1),\dots,f(2^p FFF))}.
#' @author A.I. McLeod and Y. Zhang.
#' @seealso [spectrum()], [spec.pgram()], [spec.ar()].
#' @examples 
#' ARSdf(0.8)
#' ARSdf(c(0.1,0.2))
#' 
#' @export
ARSdf <-
  function(phi, pFFT=8){
    pext<-c(1,-phi,rep(0,2^(1+pFFT) -1 -length(phi)))
    ft<-fft(pext)
    1/(Re(ft*Conj(ft)))[1:2^pFFT]
  }
