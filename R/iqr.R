iqr <- function(x) diff(quantile(x,c(0.25,0.75),na.rm=TRUE))

