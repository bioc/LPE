preprocess <- function(x, data.type = "MAS5", threshold=1,LOWESS=FALSE) {
  
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

