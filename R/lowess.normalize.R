

# Lowess normalization
lowess.normalize <- function(x,y)
{
   # x = log(cy3 or chip1) and y = log(cy5 or chip2)
   na.point <- (1:length(x))[!is.na(x) & !is.na(y)]
   x <- x[na.point]; y <- y[na.point] 
   fit <- lowess(x+y, y-x)

   diff.fit <- approx(fit$x,fit$y,x+y,ties=mean)
   out <- y - diff.fit$y
   return(out)  
}

