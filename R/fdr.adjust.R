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

        

