## The following function has been taken from 
# multtest library (by Dudoit et. al.)

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


