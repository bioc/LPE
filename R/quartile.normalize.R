
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

