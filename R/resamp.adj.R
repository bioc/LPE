
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
