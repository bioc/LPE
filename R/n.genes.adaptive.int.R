
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

