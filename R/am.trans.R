am.trans <- function(y) {
  n <- ncol(y)
  if (n < 2) {
    stop("There are no replicated arrays!")
  } 
  A <- c()
  M <- c()
  cc <- permute(1:n)
  for (i in 1:(n-1)) {
    A <- c(A, c((y + y[,cc[i,]])/2), recursive=TRUE)
    M <- c(M, c(y - y[,cc[i,]]), recursive=TRUE)
  }
  return(cbind(A,M))
}

# The above function transforms the replicated arrays
# in the (A,M) format and duplicates the data
# eg: (Y1+Y2, Y1-Y2) and (Y2+Y1, Y2-Y1) are the 
# four columns returned for input column Y1 and Y2
# i.e. duplicate arrays

