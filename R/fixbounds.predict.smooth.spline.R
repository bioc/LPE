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

