
# Finding quartile-range (iqr)
# Default is 50%, i.e. inter-quartile-range
quan.norm <- function(x,percent=50) {
  low <- 0.5*(100 - percent)/100
  high <- 0.5*(100 + percent)/100
  difference <- as.vector(diff(quantile(x, probs=c(low,high), na.rm = TRUE)))
  return(difference)
}

