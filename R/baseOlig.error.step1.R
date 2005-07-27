baseOlig.error.step1 <- function(y, stats=median, q=0.01, df=10) {

  AM <- am.trans(y)
  A <- AM[, 1]
  M <- AM[, 2]
  median.y <- apply(y,1,stats)
	
  ## baseline is calculated in two steps
  ## In the 1st step, total number os subintervals are chosen to be 100.
  quantile.A <- quantile(A, probs = seq(0, 1, q), na.rm = TRUE)
  quan.n <- length(quantile.A) - 1
  var.M <- rep(NA, length = quan.n)
  medianAs <- rep(NA, length = quan.n)
  
  # If minimum expression value(s) is more than q percentage of total
  # data, then exclude those minimum numbers for calculating quantiles

  if(sum(A == min(A)) > (q * length(A)))
    {
      tmpA <- A[!(A == min(A))]

      quantile.A <- c(min(A), quantile(tmpA, probs = seq(q, 1, q),
                                       na.rm = TRUE))
    }

  for (i in 2:(quan.n + 1))
    {
      n.i <- length(!is.na(M[A > quantile.A[i - 1] & A <= quantile.A[i]]))
      
      if (n.i >1)
        {
          mult.factor <- 0.5*((n.i - 0.5)/(n.i -1))
          
          var.M[i - 1] <- mult.factor * var(M[A > quantile.A[i - 1] & A <=
                                              quantile.A[i]], na.rm=TRUE)
          
          medianAs[i - 1] <- median(A[A > quantile.A[i - 1] & A <=
                                      quantile.A[i]], na.rm = TRUE)
        }
      
    } ## End of "for (i in 2:...)" loop
  
  
  ## In some intervals, if there are many thresholded points, then a
  ## few consecutive quantiles of A are same, causing n.i as zero!!
  ## For such intervals, assign the interval variance to the average
  ## variance of the previous and next interval if previous interval
  ## has non-NA variance, else assign it to that of next
  ## interval. Logic is that: at the high intensity region, the genes
  ## are sparse, i.e. it is unlikely that the consecutive quantiles in
  ## the max intensity region are same.
  
  if (any(is.na(var.M)) )
    {
      for (i in (quan.n-1):1 )
        {
          if (is.na(var.M[i]))
            {
              var.M[i] <- ifelse(!is.na(var.M[i-1]),
                                 mean(var.M[i+1], var.M[i-1]),
                                 var.M[i+1])
            }
        }
    }

  ## Correct the mathematical artifact of low intensity region
  var.M[1:which(var.M==max(var.M))] <- max(var.M)

  base.var <- cbind(A = medianAs, var.M = var.M)
  sm.spline <- smooth.spline(base.var[, 1], base.var[, 2], df = df)
  min.Var <- min(base.var[,2])
  var.genes <- fixbounds.predict.smooth.spline(sm.spline, median.y)$y
  if (any(var.genes <min.Var))
    var.genes[var.genes < min.Var] <- min.Var
  basevar.step1 <- cbind(A = median.y, var.M = var.genes)
  
  ord.median <- order(basevar.step1[,1])
  var.genes.ord <- basevar.step1[ord.median,]
  
  ## END of step 1 of basline caculation
  ## Here number of sub-intervals were assumed to be 100
  
  return(var.genes.ord)

}

