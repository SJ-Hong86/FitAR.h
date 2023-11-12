#' Input a Time Series
#' 
#' This function inputs time series stored in ASCII in a format that the first 
#' line in the file is a title, next few lines, beginning with a \#, are 
#' comments, and the remaining lines contain the data. Here is an example:
#' 
#' Changes In Global Temperature, Annual, 1880-1985 #Surface air 
#' "temperature change" for the globe, 1880-1985. #Degrees Celsius. 
#' "Temperature change" actually means temperature #Surface Air Temperature", 
#' `Journal of Geophysical Research`, Vol. 92, -.40 -.37 -.43 
#' .... ............... .27 .42 .02 .30 .09 .05.
#' 
#' @usage Readts(file = "", freq = 1, start = 1, VerboseQ=TRUE)
#' @param file location for input file.
#' @param freq tsp parameter, =1, annual, =12 monthly etc.
#' @param start tsp parameter.
#' @param VerboseQ normally prompt for arguments but set VerboseQ=FALSE to automate.
#' @returns ts object with attribute 'title'.
#' @author A.I. McLeod.
#' @seealso [scan()], [ts()].
#' @examples
#' # You will need to change save the data given above in a file
#' # and change the directory as appropriate
#' # z<-Readts(file="d:/datasets/mhsets/annual/globtp.1", start=1880, VerboseQ=FALSE)
#' 
#' @export
Readts <-
  function(file = "", freq = 1, start = 1,  VerboseQ=TRUE)
  {
    cat(title <- scan(file, what = "", sep = "\n", n = 1))
    commentcount <- 1
    while("#" == substring(scan(file, what = "", sep = "\n", n = 1, skip = 
                                commentcount), first = 1, last = 1)) {
      commentcount <- commentcount + 1
    }
    if (VerboseQ) {
      cat("\n start = ")
      start <- as.numeric(eval(parse(text = readline())))
    }
    if (VerboseQ) {
      cat("\n frequency = ")
      freq <- as.integer(readline())
    }
    x <- scan(file = file, skip = commentcount)
    zts <- ts(x, start = start, frequency = freq)
    i <- nchar(title)
    while(substring(title, first = i, last = i) == " ") {
      i <- i - 1
    }
    title2 <- substring(title, first = 1, last = i)
    title2
    attr(zts, "title") <- title2
    zts
  }
