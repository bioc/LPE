
######################################
##
## Local Pooled Error (LPE) test for microarray data with
## a small number of replicates.
##
######################################


#.First.lib <- function(lib, pkg) {
#  #cat("LPE library by Nitin Jain, Michael O'Connell and Jae K. Lee\n") 
#  cat("Version 1.1.2 (2003-10-29)\n") 
#  library.dynam("LPE", pkg, lib)
#  invisible()
#
#if(.Platform$OS.type=="windows" && require(Biobase) && interactive() && .Platform$GUI=="Rgui"){
#   addPDF2Vig("LPE")
#}
#
#}


permute <- function (a) {
  aa <- matrix(NA, length(a)-1, length(a))
  for (i in 1:(length(a)-1)) {
    aa[i,] <- a[c((i+1):length(a), 1:i)]
  }
  return(aa)
}

# The above function computes the all possible 
# combinations of a vector a. For example,
# for a <- 1:3, the result will be a matrix
# of 2*3:
# 2 3 1
# 3 1 2


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

fixbounds.predict.smooth.spline  <- function(object, x, deriv = 0) {
  if(missing(x)) {
    if(deriv == 0) {
      return(object[c("x", "y")])
    } else {
      x <- object$x
    }
  }
  if(is.null(object)) {
    stop("not a valid smooth.spline object")
  } else {
    out <- predict(object, x, deriv)
    maxpredY <- object$y[object$x == max(object$x)]
    out$y[out$x > max(object$x)] <- maxpredY
    minpredY <- object$y[object$x == min(object$x)]
    out$y[out$x < min(object$x)] <- minpredY
    invisible(out)
  }
}

# Above function is called by lpe function and makes sure
# that predicted values don't go negative.


# Finding quartile-range (iqr)
# Default is 50%, i.e. inter-quartile-range
quan.norm <- function(x,percent=50) {
  low <- 0.5*(100 - percent)/100
  high <- 0.5*(100 + percent)/100
  difference <- as.vector(diff(quantile(x, probs=c(low,high), na.rm = TRUE)))
  return(difference)
}

# Normalizing based on quartile range
# "percent" is the percentage of genes unchanged,
# default is 50, so the data becomes iqr normalized
quartile.normalize <- function(x,percent=50) {
  quartile <- apply(x,2,quan.norm,percent=percent)
  max.quartile <- max(quartile)
  ratio <- (quartile/max.quartile) 
  ratio.vect <- rep(ratio,nrow(x))
  adjusted <- matrix(ratio.vect, nrow=nrow(x),byrow=TRUE)
  normalized <- data.frame(x/adjusted)
  return(normalized)
}
 
# The above functions are generalized form of iqr function (below) 
iqr <- function(x) diff(quantile(x,c(0.25,0.75),na.rm=TRUE))

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



preprocess <- function(x, data.type = "MAS5", threshold=1,LOWESS=TRUE) {
  
  # Removing NA values
  x <- as.matrix(na.exclude(x))
  
  # IQR normalization
  if (data.type =="MAS4" || data.type == "MAS5") {
    x <- quartile.normalize(x, percent=50) 
  }

  # Thresholding to 'threshold' (default = 1)
  if (data.type == "MAS4" || data.type =="MAS5"|| data.type == "dChip") {
    if (length(x[x<threshold]) !=0) {
      x[x<threshold] <- threshold 
    }
  } 

  # Log based 2 transformation
  x <- logb(x,2)

  # Loess normalization of all the chips w.r.t. first one
  if (LOWESS) {
  y <- matrix(NA, nrow=nrow(x), ncol=ncol(x))
  y[,1] <- x[,1]
    for (i in 2:ncol(x)) {
      y[,i] <- lowess.normalize(x[,1],x[,i])
    }
  x <- y
  }
return(x) 
}

# Above function preprocesses the data from MAS4/5 and dchip.
# First, IQR (inter-quartile normalization) is applied to the data 
# from MAS 4/5. Then for MAS4/dChip data thresholding is applied
# at 1 and for MAS5 data, thresholding is applied at 0.1
# Finally, the data is log transformed to the base 2. 

n.genes.adaptive.int <- function(baseOlig.error.step1.res, 
	min.genes.int = 10, div.factor =1){
 # a) min # genes in each interval = 10
  # b) max # genes in each interval = (total genes)/100
  # c) genes are chosen s.t. they are within [median + fraction(S.D.)] from
  #	the starting gene of each interval, if above coonditions are met.
  
  y <- baseOlig.error.step1.res
  max.genes.int <- round(nrow(y)/100, digits=0)
  genes.sub.int <- c()
  i <- 1  
  
  while(i < nrow(y)){
    temp.median <- y[i,1]
    temp.sd <- sqrt(y[i,2])/div.factor
    n.genes.int <- length(which(y[,1] <= temp.median + temp.sd))-(i-1)
    if (n.genes.int < min.genes.int) n.genes.int <- min.genes.int
    if (n.genes.int > max.genes.int) n.genes.int <- max.genes.int
    genes.sub.int <- c(genes.sub.int,n.genes.int)
   i <- i+n.genes.int
  } 
  
  # Checking if "extra" genes have been assigned to the last interval  
  if(sum(genes.sub.int)!=nrow(y)){
    len.temp <- length(genes.sub.int)
    temp1 <- genes.sub.int[1:(len.temp-1)]
    extra.genes <- nrow(y)-sum(temp1) 
    temp1[len.temp-1] <- temp1[len.temp-1] + extra.genes
    genes.sub.int <- temp1
    rm(temp1,len.temp,extra.genes)
  }
  
 return(genes.sub.int)

}


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
  # data, then exclude those minimum numbers for
  # calculating quantiles

  if(sum(A == min(A)) > (q * length(A))) {
    tmpA <- A[!(A == min(A))]
    quantile.A <- c(min(A), quantile(tmpA, probs = seq(q, 1, q), na.rm = TRUE))
  }

  for(i in 2:(quan.n + 1)) {
    n.i <- length(!is.na(M[A > quantile.A[i - 1] & A <= quantile.A[i]]))
    mult.factor <- 0.5*((n.i - 0.5)/(n.i -1))
    var.M[i - 1] <- mult.factor * var(M[A > quantile.A[i - 1] & A <= quantile.A[i]], na.rm=TRUE) 
    medianAs[i - 1] <- median(A[A > quantile.A[i - 1] & A <= quantile.A[i]], na.rm = TRUE)
  }
  var.M[1:which(var.M==max(var.M))] <- max(var.M)
  
  base.var <- cbind(A = medianAs, var.M = var.M)
  sm.spline <- smooth.spline(base.var[, 1], base.var[, 2], df = df)
  var.genes <- fixbounds.predict.smooth.spline(sm.spline, median.y)$y
  basevar.step1 <- cbind(A = median.y, var.M = var.genes)
  
  ord.median <- order(basevar.step1[,1])
  var.genes.ord <- basevar.step1[ord.median,]
  
  ## END of step 1 of basline caculation
  ## Here number of sub-intervals were assumed to be 100
  
  return(var.genes.ord)

}

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
  var.genes.adap <- fixbounds.predict.smooth.spline(sm.spline.adap, median.y)$y
  basevar.all.adap <- cbind(A = median.y, var.M = var.genes.adap)
  
 return(basevar.all.adap)

#  return(cbind(A = medianAs, var.M = var.M))
}

 

baseOlig.error <- function(y, stats=median, q = 0.01, min.genes.int =10, div.factor=1) {
 # Calls baseOlig.error.step1 and baseOlig.error.step2 functions
 
 baseline.step1 <- baseOlig.error.step1(y, stats=stats)
 baseline.step2 <- baseOlig.error.step2(y, stats=stats, baseline.step1, min.genes.int=min.genes.int, div.factor=div.factor)
 return(baseline.step2) 
  
}

# The above function evaluates baseline distribution of M at percentile intervals of A.
# y is (log transformed intensity) of replicated Oligo arrays after normalization and
# q = quantile width 

lpe <- function(x, y, basevar.x, basevar.y, df = 10, array.type = "olig",
                probe.set.name = "OLIG.probe.name",
	        trim.percent = 5) {

  # LPE significance evaluation with replicates
  # x and y are two array samples with n1 and n2 replicates
  # basevar.x (basevar.y): n.quantile x 2 matrix of LPE baseline error of x (y)
  # array.type: "olig" for Affymetrix oligo array and "cDNA" for cDNA array
  # subset rows, removing rows with any na's

  # trim.size <- round((trim.percent/100) * nrow(basevar.x), digits=0)
  # basevar.x <- basevar.x[(trim.size + 1): nrow(basevar.x),  ]
  # basevar.y <- basevar.y[(trim.size + 1): nrow(basevar.y),  ]
  # express.df <- cbind(x, y)
  # express.df <- na.exclude(express.df)

  # x <- as.matrix(express.df[, 1:ncol(x)])
  # y <- as.matrix(express.df[, (ncol(x) + 1):ncol(express.df)] )

  # If some rows containing NAs were removed then, remove corresponding IDs too

  # pickoff <- attr(express.df, "na.action")
  # if(!is.null(pickoff)) {
  #   probe.set.name <- probe.set.name[-pickoff]
  # }
 
  n1 <- ncol(x)
  n2 <- ncol(y)
  ngenes <- nrow(x)

  if (n1 < 2 | n2 < 2) {
    stop("No replicated arrays!")
  }
  if (n1 > 2 |n2 >2){
    # median.x <- apply(x, 1, median)
    # median.y <- apply(y, 1, median)
    # sf.x <- smooth.spline(basevar.x[, 1], basevar.x[, 2], df = df)
    # var.x <- fixbounds.predict.smooth.spline(sf.x, median.x)$y
    # sf.y <- smooth.spline(basevar.y[, 1], basevar.y[, 2], df = df)
    # var.y <- fixbounds.predict.smooth.spline(sf.y, median.y)$y
    var.x <- basevar.x[,2]
    var.y <- basevar.y[,2]
    median.x <- basevar.x[,1]
    median.y <- basevar.y[,1]
    median.diff <- median.x - median.y

    std.dev <- sqrt(1.57 * ((var.x/n1) + (var.y/n2)))
    z.stats <- median.diff/std.dev
    
    data.out <- data.frame(x=x, median.1 = median.x, std.dev.1 = sqrt(var.x),
			   y=y, median.2 = median.y, std.dev.2 = sqrt(var.y),
       		           median.diff = median.diff, pooled.std.dev = std.dev,
			   z.stats=z.stats)
    row.names(data.out) <- probe.set.name
    
    return(data.out)
  } 
  if (n1 ==2 & n2 ==2) {
    # median.x <- (x[, 1] + x[, 2])/2
    # median.y <- (y[, 1] + y[, 2])/2
    
    # sf.x <- smooth.spline(basevar.x[, 1], basevar.x[, 2], df = df)
    # var.x <- fixbounds.predict.smooth.spline(sf.x, median.x)$y
    # sf.y <- smooth.spline(basevar.y[, 1], basevar.y[, 2], df = df)
    # var.y <- fixbounds.predict.smooth.spline(sf.y, median.y)$y
    var.x <- basevar.x[,2]
    var.y <- basevar.y[,2]
    median.x <- basevar.x[,1]
    median.y <- basevar.y[,1]
    median.diff <- median.x- median.y

    std.dev <- sqrt((var.x/n1) + (var.y/n2))
    z.stats <- median.diff/std.dev
    pnorm.diff <- pnorm(median.diff, mean = 0, sd = std.dev)
    p.out <- 2 * apply(cbind(pnorm.diff, 1 - pnorm.diff), 1, min)
   
    # Outlier checking 
    var.xo <- var.yo <- matrix(NA, ngenes, 2)
    for(i in 1:2) {
      # prediction of variance for each x[,i] and each y[,i]
      # rather than for the median/mean of these
      sf.xi <- smooth.spline(basevar.x[, 1], basevar.x[, 2], df = df)
      var.xo[, i] <- fixbounds.predict.smooth.spline(sf.xi, x[, i])$y
      sf.yi <- smooth.spline(basevar.y[, 1], basevar.y[, 2], df = df)
      var.yo[, i] <- fixbounds.predict.smooth.spline(sf.yi, y[, i])$y
    }

    p.val <- matrix(NA, ngenes, 2)
    var.diff <- (var.xo[, 1] + var.xo[, 2] + var.yo[, 1] + var.yo[, 2])/4
    diff.xy <- x - y
    diff.xy <- (diff.xy[, 1] - diff.xy[, 2])
    ### 
    p.val[, 1] <- pnorm(diff.xy, mean = 0, sd = sqrt(var.diff))
    p.val[, 1] <- apply(cbind(p.val[, 1], 1 - p.val[, 1]), 1, min)

    diff.xy <- x - y[, 2:1]
    diff.xy <- (diff.xy[, 1] - diff.xy[, 2])

    p.val[, 2] <- pnorm(diff.xy, mean = 0, sd = sqrt(var.diff))
    p.val[, 2] <- apply(cbind(p.val[, 2], 1 - p.val[, 2]), 1, min)

    p.outlier <- apply(p.val, 1, min)
    flag <- rep(".", ngenes)
    flag[p.outlier < 0.05] <- "*"
    flag[p.outlier < 0.01] <- "**"
    flag[p.outlier < 0.001] <- "***"

    data.out <- data.frame(x=x, median.1=median.x, std.dev.1 = sqrt(var.x),
			   y=y, median.2=median.y, std.dev.2 = sqrt(var.y),
        	           median.diff = median.diff, pooled.std.dev=std.dev,
   			  z.stats=z.stats, flag, p.outlier)
    row.names(data.out) <- probe.set.name
 
    return(data.out)
  }
}

### NEW FUNCTION
fdr.adjust <- function(lpe.result,adjp="resamp",target.fdr=c(10^-3 ,seq(0.01,0.10,0.01), 0.15, 0.20, 0.50), iterations=5, ALL=FALSE){

if(adjp=="resamp"){
    x.location <- grep("^x",names(lpe.result))
    y.location <- grep("^y",names(lpe.result))
    x <- lpe.result[,x.location]
    y <- lpe.result[,y.location]

    z.null.resamp <- resamp.adj(x, y, q=0.01, iterations=iterations)
    z.null.iterations <- abs(z.null.resamp)

    num.genes <- nrow(x)
    # Computing FDR for all z-values
    z.real <- sort(abs(lpe.result$z.stats), decreasing = TRUE)
    num.all <-  length(z.null.iterations)

    min.fdr <- 1.0/num.all
#   pi0 <- length(which(z.real <= median(z.null.iterations)))/(length(z.real)/2.0)
    pi0 <- length(which(z.real <= quantile(z.null.iterations,0.90)))/
		 (length(z.real)*0.9)

    pi0 <- min(1,pi0)

      FDR <- rep(NA, num.genes)
      cat("Computing FDR...\n")
      n.col.iterations <-  ncol(z.null.iterations)
     
    for (i in 1:num.genes) {
      n.false.pos <- length(which(z.null.iterations >= z.real[i]))
      avg.n.false.pos <- n.false.pos/n.col.iterations 
      n.total.pos <- length(which(z.real >= z.real[i]))
      temp <- avg.n.false.pos/n.total.pos
      FDR[i] <- max(min(pi0*temp, 1),min.fdr)
    }

     for (j in num.genes:1) {
      FDR[j] <- min(FDR[j:num.genes]) # enforcing monotonicity
    }
   
    if(!ALL) {
      # Computing FDR for only desired values (target.fdr)
      z.critical <- rep(NA, length(target.fdr))
      for (j in 1:length(target.fdr)) {
        temp <- (target.fdr[j] >= FDR)
        genes.signif <- (1:num.genes)[temp]
        num.genes.signif <- length(genes.signif)
        z.critical[j] <- z.real[num.genes.signif]
      }  # end of ' for j in target.fdr..' loop
      data.out <- cbind(target.fdr,z.critical)
      } else { # end of 'if(!ALL..' loop
        data.out <- cbind(FDR,z.real)
      } # end of else loop

      return(data.out)
    }  # end of 'adjp==resamp' loop
    
   if (adjp == "BH" || adjp=="BY") {
      x.location <- grep("^x",names(lpe.result))
      y.location <- grep("^y",names(lpe.result))
      x <- lpe.result[,x.location]
      y <- lpe.result[,y.location]
      pnorm.diff <- pnorm(lpe.result$median.diff, mean = 0, 
	  sd = lpe.result$pooled.std.dev)
      p.out <- 2 * apply(cbind(pnorm.diff, 1 - pnorm.diff), 1, min)
      p.adj <- mt.rawp2adjp(p.out, proc=adjp)
      data.out <- data.frame(x=x, median.1 = lpe.result$median.1, 
		std.dev.1 = lpe.result$std.dev.1, y=y, 
		median.2 = lpe.result$median.2,
		std.dev.2 = lpe.result$std.dev.2, 
		median.diff = lpe.result$median.diff, 
		pooled.std.dev= lpe.result$pooled.std.dev,
		abs.z.stats=abs(lpe.result$z.stats),p.adj=p.adj)
        col.id <- grep(adjp, colnames(data.out))
 	aa <- cbind(FDR=data.out[,col.id],
		 z.real=data.out$abs.z.stats)
	aa <- aa[order(aa[,2],decreasing=TRUE),]
       return(aa)
    } # End of adjp=="BH" || adjp=="BY"

if (adjp == "mix.all") { # this is like SAM strategy
      x.location <- grep("^x", names(lpe.result))
      y.location <- grep("^y", names(lpe.result))
      x <- lpe.result[, x.location]
      y <- lpe.result[, y.location]
      orig.data <- as.matrix(cbind(x,y))
      # Sampling all genes to generate NULL data
      z.stats.null <- matrix(NA,nrow(x), ncol=iterations)
      for (m in 1:iterations) {
        null.data <- matrix(sample(orig.data), nrow=nrow(orig.data),
			ncol=ncol(orig.data))
        basevar.x.null <- baseOlig.error(null.data[, 1:ncol(x)])
        basevar.y.null <- baseOlig.error(null.data[, (1 + ncol(x)):(ncol(x) + 
            ncol(y))])
        lpe.null <- lpe(null.data[, 1:ncol(x)], null.data[, (1 + 
            ncol(x)):(ncol(x) + ncol(y))], basevar.x.null, 
	    basevar.y.null, probe.set.name = as.character(1:nrow(x)))
        z.stats.null[,m] <- (lpe.null$z.stats)
        cat("iteration number",m," finished \n")
      }
      
      z.stats.null.abs <- abs(z.stats.null)
      z.real <- sort(abs(lpe.result$z.stats), decreasing = TRUE)
      num.all <- length(z.stats.null) # the number of genes in null
      z.null.iterations <- z.stats.null.abs

     min.fdr <- 1/num.all
#     pi0 <- length(which(z.real <= median(z.null.iterations)))/(length(z.real)/2)
 pi0 <- length(which(z.real <= quantile(z.null.iterations,0.90)))/(length(z.real)*0.9)    
        pi0 <- min(1, pi0)
        num.genes <- length(z.real)

        FDR <- rep(NA, num.genes)
        cat("Computing FDR...\n")
        n.col.iterations <-  ncol(z.null.iterations)

        for (i in 1:num.genes) {
          n.false.pos <- length(which(z.null.iterations >= z.real[i]))
          avg.n.false.pos <- n.false.pos/n.col.iterations 
          n.total.pos <- length(which(z.real >= z.real[i]))
          temp <- avg.n.false.pos/n.total.pos
          FDR[i] <- max(min(pi0 * temp, 1), min.fdr)
          FDR[i] <- max(FDR[1:i]) # enforcing monotonicity
        } # end of 'for i in num.genes..' loop

        if (!ALL) {
            z.critical <- rep(NA, length(target.fdr))
            for (j in 1:length(target.fdr)) {
                temp <- (target.fdr[j] >= FDR)
                genes.signif <- (1:num.genes)[temp]
                num.genes.signif <- length(genes.signif)
                z.critical[j] <- z.real[num.genes.signif]
            }
            data.out <- cbind(target.fdr, z.critical)
        } else { # end of 'if(!ALL..'
          data.out <- cbind(FDR, z.real)
        } # end of 'else ..'
       return(data.out)
     } # END mix.all
} # END fdr.adjust

        

######################
resamp.adj <- function(x, y, q=0.01, iterations=5, min.genes.int=10) {
 median.x <- apply(x,1,median)
 rank.x <- rank(median.x)
 median.y <- apply(y,1,median)
 rank.y <- rank(median.y)
 diff.rank <- abs(rank.x-rank.y)

 # finding overall median of genes across all conditions
 data.all <- cbind(x,y)
 median.all <- apply(data.all,1, median) 
 data.median <- cbind(median.all, diff.rank)
 # sorting the median data
 data.median.ord <- data.median[order(median.all),]

 # Finding the quantiles of overall median
 #q.all <- quantile(median.all, probs=seq(0,1,q))
 #n.gene <- nrow(x)
 
 # Finding the baseline of overall data
 baseline.all.step1 <- baseOlig.error.step1(data.all)

 ## Obtain the number of intervals (like in step 2 of baseOlig.error function

  genes.sub.int <- n.genes.adaptive.int(baseline.all.step1, 
      min.genes.int = min.genes.int, div.factor = 1)

   n.genes.invar.int <- paste("n.genes.invar.int", 1:length(genes.sub.int), 
        sep = "")
    invar.genes.all <- c()
    keep.insig.genes <- 0.90
   
 start <- 1
    for (i in 1:length(genes.sub.int)) {
        temp1 <- genes.sub.int[i]
        temp <- start:(start + temp1 - 1)
        data.interval <- data.median.ord[temp, ]
        len.invar.genes <- round(keep.insig.genes* nrow(data.interval),
		 digits = 0)
        assign(n.genes.invar.int[i], len.invar.genes)
        rank.diff.interval <- order(data.interval[, 2])
        invar.genes.interval <- rank.diff.interval[1:len.invar.genes]
        invar.genes.all <- c(invar.genes.all, data.interval[invar.genes.interval, 
            1])
        start <- start + temp1
    } # end of 'for (i in... ' loop 

  genes.invar <- (1:nrow(x))[row.names(x) %in% names(invar.genes.all)]
    n.genes.total <- nrow(x)
   
    n.genes.invar <- length(genes.invar)
    orig.data <- data.frame(x, y)
   
 invar.data <- orig.data[genes.invar, ]
    invar.median <- apply(invar.data, 1, median)
    invar.median.order <- order(invar.median)
    invar.data.ord <- invar.data[invar.median.order, ]
    null.data <- data.frame(matrix(NA, nrow = n.genes.total, 
        ncol = (ncol(x) + ncol(y))))
    
 intervals <- length(genes.sub.int)
    genes.sub.int.invar <- c()
    for (i in 1:intervals) {
        genes.sub.int.invar <- c(genes.sub.int.invar, (get(n.genes.invar.int[i])))
    }

    # sum(genes.sub.int.invar)
    z.stats.null.iterations <- matrix(NA, nrow = nrow(x), ncol = iterations)

 for (m in 1:iterations) {
   cat("iteration number", m, "is in progress\n")

### Generating null distribution 
   start <- 1 # genes in interval
   start1 <- 1 # invar genes
   finish <- genes.sub.int[1] # num. of all genes
   finish1 <- genes.sub.int.invar[1] # num. of invar.genes
   new.intervals <- length(genes.sub.int.invar)

     for (i in 1:new.intervals) { 
      temp1 <- start1:(start1+finish1-1)
      invar.data.int <- invar.data.ord[temp1,]
      invar.data.x <- as.vector(as.matrix(invar.data.int[,1:ncol(x)]))
      invar.data.y <- as.vector(as.matrix(invar.data.int[,(1+ncol(x)):
			(ncol(x)+ncol(y))]))
      pooled.data <- c(invar.data.x, invar.data.y)
        null.data[start:(start+finish-1), ] <-
    	 sample(pooled.data, size=finish*(ncol(x)+ncol(y)), replace=TRUE)
        start <- start+ finish				
        start1 <- start1+ finish1	# invar			
        finish1 <- genes.sub.int.invar[i+1]
        finish <- genes.sub.int[i+1]
      } # end of 'for i in new.intervals..' loop

 # Calculating the z-statistics of 'NULL' distribution:
      # Significant genes in 'NULL' distribution are "FALSE-POSITIVES"

      basevar.x <- baseOlig.error(null.data[, 1:ncol(x)])
      basevar.y <- baseOlig.error(null.data[, (1 + ncol(x)):(ncol(x)+ncol(y))])
      lpe.null <- lpe(null.data[, 1:ncol(x)], null.data[, (1 +ncol(x)):
	     (ncol(x) + ncol(y))], basevar.x, basevar.y, 
              probe.set.name = as.character(1:nrow(x)))

      z.stats.null <- lpe.null$z.stats
      #bitmap("z.null.resamp.png", type="png256")
      #plot(z.stats.null, pch=16, cex=0.4)
      #title("Distribution of z-null Resampling based technique")
      #dev.off()
      z.stats.null.iterations[,m] <- z.stats.null
      cat("iteration number", m, "finished\n")
    } # end of 'for m in iterations..' loop

    return(z.stats.null.iterations)
}

#####

mt.rawp2adjp <-  function (rawp, proc = c("Bonferroni", "Holm", 
            		   "Hochberg", "SidakSS", "SidakSD", "BH", "BY")) {
  m <- length(rawp)
  n <- length(proc)

  index <- order(rawp)
  spval <- rawp[index]

  adjp <- matrix(0, m, n + 1)
  dimnames(adjp) <- list(NULL, c("rawp", proc))
  adjp[, 1] <- spval
  if (is.element("Bonferroni", proc)) {
    tmp <- m * spval
    tmp[tmp > 1] <- 1
    adjp[, "Bonferroni"] <- tmp
  }
  if (is.element("Holm", proc)) {
    tmp <- spval
    tmp[1] <- min(m * spval[1], 1)
    for (i in 2:m) {
      tmp[i] <- max(tmp[i - 1], min((m - i + 1) * spval[i], 1))
    } 
    adjp[, "Holm"] <- tmp
  }
  if (is.element("Hochberg", proc)) {
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((m - i + 1) * spval[i], 1))
    }
    adjp[, "Hochberg"] <- tmp
  }
  if (is.element("SidakSS", proc)) {
    adjp[, "SidakSS"] <- 1 - (1 - spval)^m
  }
  if (is.element("SidakSD", proc)) {
    tmp <- spval
    tmp[1] <- 1 - (1 - spval[1])^m
    for (i in 2:m) {
      tmp[i] <- max(tmp[i - 1], 1 - (1 - spval[i])^(m - i + 1))
    }
    adjp[, "SidakSD"] <- tmp
  }
  if (is.element("BH", proc)) {
    tmp <- spval
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((m/i) * spval[i], 1))
    }
    adjp[, "BH"] <- tmp
  }
  if (is.element("BY", proc)) {
    tmp <- spval
    a <- sum(1/(1:m))
    tmp[m] <- min(a * spval[m], 1)
    for (i in (m - 1):1) {
      tmp[i] <- min(tmp[i + 1], min((m * a/i) * spval[i], 1))
    }
    adjp[, "BY"] <- tmp
  }

  # inversing the sort
  ndx <- order(index) 
  adjp <- adjp[ndx,]
  
  list(adjp = adjp, index = index)
}

# The above function "mt.rawp2adjp" was taken from Bioconductor (Dudoit et al)
