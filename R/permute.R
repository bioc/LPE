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

