baseOlig.error <- function(y, stats=median, q = 0.01, min.genes.int =10, div.factor=1) {
 # Calls baseOlig.error.step1 and baseOlig.error.step2 functions
 
 baseline.step1 <- baseOlig.error.step1(y, stats=stats)
 baseline.step2 <- baseOlig.error.step2(y, stats=stats, baseline.step1, min.genes.int=min.genes.int, div.factor=div.factor)
 return(baseline.step2) 
  
}



# The above function evaluates baseline distribution of M at percentile intervals of A.
# y is (log transformed intensity) of replicated Oligo arrays after normalization and
# q = quantile width 
