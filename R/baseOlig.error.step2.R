baseOlig.error.step2 <- function(y, baseOlig.error.step1.res, df=10, stats=median, min.genes.int=10, div.factor=1) {
 AM <- am.trans(y)
  A <- AM[, 1]
  M <- AM[, 2]
  median.y <- apply(y,1,stats)

var.genes.ord <- baseOlig.error.step1.res	
genes.sub.int <- n.genes.adaptive.int(var.genes.ord,min.genes.int=min.genes.int, div.factor=div.factor)

  ## Re-calculating the baseline distribution based on new adaptive intervals

  j.start <- 1
  j.end <- 0
  var.M.adap <- rep(NA, length = length(genes.sub.int))
  medianAs.adap <- rep(NA, length = length(genes.sub.int))
    
  for (i in 2:(length(genes.sub.int)+1)) {
        j.start <- j.end + 1
        j.end <- j.start+genes.sub.int[i-1]-1
	
	vect.temp <- (A > var.genes.ord[j.start,1] & 
			A <= var.genes.ord[j.end,1])
	n.i <- length(!is.na(M[vect.temp]))
        mult.factor <- 0.5 * ((n.i - 0.5)/(n.i - 1))
        var.M.adap[i - 1] <- mult.factor * var(M[vect.temp], na.rm = TRUE)
        medianAs.adap[i - 1] <- median(A[vect.temp], na.rm = TRUE)
  }

  var.M.adap[1:which(var.M.adap == max(var.M.adap))] <- max(var.M.adap)

  base.var.adap <- cbind(A.adap = medianAs.adap, var.M.adap = var.M.adap)
  sm.spline.adap <- smooth.spline(base.var.adap[, 1], base.var.adap[, 2], df = df)
  min.Var <- min(base.var.adap[,2])
  var.genes.adap <- fixbounds.predict.smooth.spline(sm.spline.adap, median.y)$y
  var.genes.adap[var.genes.adap < min.Var] <- min.Var
  basevar.all.adap <- cbind(A = median.y, var.M = var.genes.adap)
  
 return(basevar.all.adap)

#  return(cbind(A = medianAs, var.M = var.M))
}

 

