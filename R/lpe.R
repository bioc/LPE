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

preprocess <- function(x, data.type = "MAS5") {
  x <- as.matrix(na.exclude(x))
  if (data.type =="MAS4" || data.type == "MAS5") {
    x <- quartile.normalize(x, percent=50)
    #adj <- matrix(rep(apply(x,2,iqr)/max(apply(x,2,iqr)),
    #			nrow(x)),nrow=nrow(x), byrow=TRUE)
    #x <- data.frame(x/adj)
    #rm(adj)
  }
  if (data.type == "MAS4" || data.type == "dChip") {
    if (length(x[x<1]) !=0) {
      x[x<1] <- 1 
    }
  } 
  if (data.type == "MAS5") {
    if (length(x[x <0.1]) !=0) {
      x[x < 0.1] <- 0.1 
    }
  }
  x <- logb(x,2)
}

# Above function preprocesses the data from MAS4/5 and dchip.
# First, IQR (inter-quartile normalization) is applied to the data 
# from MAS 4/5. Then for MAS4/dChip data thresholding is applied
# at 1 and for MAS5 data, thresholding is applied at 0.1
# Finally, the data is log transformed to the base 2. 

baseOlig.error <- function(y, q = 0.01) {

  AM <- am.trans(y)
  A <- AM[, 1]
  M <- AM[, 2]
	
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
  return(cbind(A = medianAs, var.M = var.M))
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

  trim.size <- round((trim.percent/100) * nrow(basevar.x), digits=0)
  basevar.x <- basevar.x[(trim.size + 1): nrow(basevar.x),  ]
  basevar.y <- basevar.y[(trim.size + 1): nrow(basevar.y),  ]
  express.df <- cbind(x, y)
  express.df <- na.exclude(express.df)

  x <- as.matrix(express.df[, 1:ncol(x)])
  y <- as.matrix(express.df[, (ncol(x) + 1):ncol(express.df)] )

  # If some rows containing NAs were removed then, remove corresponding IDs too

  pickoff <- attr(express.df, "na.action")
  if(!is.null(pickoff)) {
    probe.set.name <- probe.set.name[-pickoff]
  }
 
  n1 <- ncol(x)
  n2 <- ncol(y)
  ngenes <- nrow(x)

  if (n1 < 2 | n2 < 2) {
    stop("No replicated arrays!")
  }
  if (n1 > 2 |n2 >2){
    median.x <- apply(x, 1, median)
    median.y <- apply(y, 1, median)
    median.diff <- median.x - median.y
    sf.x <- smooth.spline(basevar.x[, 1], basevar.x[, 2], df = df)
    var.x <- fixbounds.predict.smooth.spline(sf.x, median.x)$y
    sf.y <- smooth.spline(basevar.y[, 1], basevar.y[, 2], df = df)
    var.y <- fixbounds.predict.smooth.spline(sf.y, median.y)$y
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
    median.x <- (x[, 1] + x[, 2])/2
    median.y <- (y[, 1] + y[, 2])/2
    median.diff <- median.x- median.y

    sf.x <- smooth.spline(basevar.x[, 1], basevar.x[, 2], df = df)
    var.x <- fixbounds.predict.smooth.spline(sf.x, median.x)$y
    sf.y <- smooth.spline(basevar.y[, 1], basevar.y[, 2], df = df)
    var.y <- fixbounds.predict.smooth.spline(sf.y, median.y)$y

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
fdr.adjust <- function(lpe.result,adjp="resamp",target.fdr=c(10^-3 ,seq(0.01,0.10,0.01), 0.15, 0.20, 0.50), iterations=5){
    
if(adjp=="resamp"){
    x.location <- grep("^x",names(lpe.result))
    y.location <- grep("^y",names(lpe.result))
    x <- lpe.result[,x.location]
    y <- lpe.result[,y.location]
    z.null.iterations <- resamp.adj(x,y, q=0.01, iterations=iterations)

    num.genes <- nrow(x)
    # Computing FDR for all z-values
    z.real <- sort(abs(lpe.result$z.stats), decreasing = TRUE)
    num.all <-  length(z.null.iterations)
    FDR <- c(NA, num.genes)
    cat("Computing FDR...\n")
    for (i in 1:num.genes) {
      false.pos <- (1:num.all)[z.null.iterations >= z.real[i]]
      n.false.pos <- length(false.pos)
      avg.n.false.pos <- n.false.pos/ncol(z.null.iterations)
      true.pos <- (1:length(z.real))[z.real >= z.real[i]]
      n.true.pos <- length(true.pos)
      FDR[i] <- avg.n.false.pos/n.true.pos
      FDR[i] <- min(FDR[i], 1)
    }
   
    # Computing FDR for only desired values (target.fdr)
    z.critical <- c(NA, length(target.fdr))
    for (j in 1:length(target.fdr)) {
      temp <- (target.fdr[j] >= FDR)
      genes.signif <- (1:num.genes)[temp]
      num.genes.signif <- length(genes.signif)
      z.critical[j] <- z.real[num.genes.signif]
    }

    data.out <- cbind(target.fdr,z.critical)
     
    } else if(adjp!="resamp"){
      x.location <- grep("^x",names(lpe.result))
      y.location <- grep("^y",names(lpe.result))
      x <- lpe.result[,x.location]
      y <- lpe.result[,y.location]

      pnorm.diff <- pnorm(lpe.result$median.diff, mean = 0, sd = lpe.result$pooled.std.dev)
      p.out <- 2 * apply(cbind(pnorm.diff, 1 - pnorm.diff), 1, min)
      p.adj <- mt.rawp2adjp(p.out, proc=adjp)
      data.out <- data.frame(x=x, median.1 = lpe.result$median.1, std.dev.1 = 
		lpe.result$std.dev.1, y=y, median.2 = lpe.result$median.2,
		std.dev.2 = lpe.result$std.dev.2, median.diff = 
		lpe.result$median.diff, pooled.std.dev=
		lpe.result$pooled.std.dev,abs.z.stats=abs(lpe.result$z.stats),
    		p.adj=p.adj)
    }
}

        

######################
resamp.adj <- function(x, y, q=0.01, iterations=5) {
 median.x <- apply(x,1,median)
 rank.x <- rank(median.x)
 median.y <- apply(y,1,median)
 rank.y <- rank(median.y)
 diff.rank <- abs(rank.x-rank.y)/length(rank.x)

 # finding overall median of genes across all conditions
 median.all <- apply(cbind(x,y),1, median) 
 data.median <- cbind(median.all, diff.rank)
 # sorting the median data
 data.median.ord <- data.median[order(median.all),]
 # Finding the quantiles of overall median
 q.all <- quantile(median.all, probs=seq(0,1,q))
 n.gene <- nrow(x)
 
 
 # finding the 50% rank-invariant genes from each quantile
 intervals <- as.integer(1/q)
 invar.genes.all <- c()
 for (i in 1:intervals) {
   temp <- (data.median.ord[,1] >= q.all[i]) & (data.median.ord[,1] < q.all[i+1])
   data.interval <- data.median.ord[temp,]
   len.interval <- round(0.5*nrow(data.interval),digits=0)
   rank.diff.interval <- order(data.interval[, 2])
   invar.genes.interval <- rank.diff.interval[1:len.interval]
   invar.genes.all <- c(invar.genes.all, data.interval[invar.genes.interval,1])
}

 # finding the order of invariant genes
 genes.invar <- (1:nrow(x))[row.names(x)%in%names(invar.genes.all)]
 n.genes.invar <- length(genes.invar)
 
 orig.data <- data.frame(x,y)
 invar.data <- orig.data[genes.invar,]
 median.new.x <- median.x[genes.invar]
 median.new.y <- median.y[genes.invar]
 
 # finding the quantiles of rank-invariant gene medians
 quantile.x <- quantile(median.new.x, probs=seq(0.00,1,q))
 quantile.y <- quantile(median.new.y, probs=seq(0.00,1,q))

 z.stats.null.iterations <- matrix(NA, nrow=nrow(x), ncol=iterations)
 for (m in 1:iterations) {
   cat("iteration number", m, "is in progress\n")
   pooled.data <- paste("pooled.data",1:intervals,sep="")

   # Pooling the quantiles of invariant genes 
   for (i in 1:intervals) {
     temp1 <- (1:n.genes.invar)[((median.new.x > quantile.x[i]) &
                   (median.new.x <= quantile.x[i+1]))]
     	
     temp2 <- (1:n.genes.invar)[((median.new.y > quantile.y[i]) &
     	           (median.new.y <= quantile.y[i+1]))]
     
     pool.ref <- as.vector(as.matrix(invar.data[temp1,1:ncol(x)]))
     pool.targ <- as.vector(as.matrix(invar.data[temp2,((1+ncol(x)):(ncol(x)+ncol(y)))]))
     
     assign(pooled.data[i], c(pool.ref,pool.targ))
   }
   # Above loop assigns pooled quantiles in pooled.data

   # Generating the 'NULL' distribution:
   null.data <- matrix(NA, nrow=nrow(x), ncol=(ncol(x)+ ncol(y)))
   num.row.per.quant <- nrow(x)%/%intervals
   for (i in 1:nrow(null.data)) {
     j <- 1+ (i%/%num.row.per.quant)
     if (j > intervals) {
       j <- sample(1:intervals,1)
     }
     temp <- get(pooled.data[j])
     null.data[i,] <- sample(temp,(ncol(x)+ncol(y)),replace=TRUE)
   }
   
   # Calculating the z-statistics of 'NULL' distribution:
   # Remember, all significant genes in this 'NULL' distribution
   # are "FALSE-POSITIVES"
   basevar.x <- baseOlig.error(null.data[,1:ncol(x)])
   basevar.y <- baseOlig.error(null.data[,(1+ncol(x)):(ncol(x)+ncol(y))])

   lpe.null <- lpe(null.data[,1:ncol(x)],null.data[,(1+ncol(x)):(ncol(x)+ncol(y))], 
  		basevar.x, basevar.y, probe.set.name=as.character(1:nrow(x)))

   z.stats.null <- lpe.null$z.stats
   z.stats.null.iterations[,m] <- abs(z.stats.null)
   cat("iteration number", m, "finished\n")
 }
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
