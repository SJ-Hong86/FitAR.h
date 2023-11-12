#' Converts a Matrix from Symmetric Storage Mode to Regular Format
#' 
#' Utility function.
#' 
#' @usage FromSymmetricStorageUpper(x)
#' @param x a vector which represents a matrix in upper triangular form.
#' @returns symmetric matrix.
#' @author A.I. McLeod.
#' @examples 
#' FromSymmetricStorageUpper(1:5)
#' 
#' @export
FromSymmetricStorageUpper <-
  function(x){
    n<-floor((-1+sqrt(1+8*length(x)))/2)
    z<-matrix(numeric(n^2), nrow=n)
    i<-as.vector(lower.tri(z,diag=TRUE))
    z[i]<-x
    ztranspose<-t(z)
    diag(ztranspose)<-0
    z+ztranspose
  }
