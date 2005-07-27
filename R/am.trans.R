am.trans <- function (y)
  {
    if(ncol(y)>5)
      y <- y[,sample(1:ncol(y),5)]
    n <- ncol(y)

    if (n < 2)
      {
        stop("There are no replicated arrays!")
      }

    A <- c()
    M <- c()
    cc <- permute(1:n)

    for (i in 1:(n - 1))
      {
        A <- c(A, c((y + y[, cc[i, ]])/2), recursive = TRUE)
        M <- c(M, c(y - y[, cc[i, ]]), recursive = TRUE)
      }
    return(cbind(A, M))
  }
