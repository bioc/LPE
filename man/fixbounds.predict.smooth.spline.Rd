\name{fixbounds.predict.smooth.spline}
\alias{fixbounds.predict.smooth.spline}

\title{
	Makes the predicted variance non negative
}

\description{
	Makes the predicted variance non negative
}

\usage{
fixbounds.predict.smooth.spline(object,x, deriv=0)
}

\arguments{
  \item{object}{variance from baseOlig.error function}
  \item{x}{vector for which variance needs to be predicted}
  \item{deriv}{derivative of the vetor required, default =0}
}

\value{
  Returns the predicted variance for the given vector based on the
  baseline error distribution. Maximum and minimum predicted values 
  for the vetor are same as those of baseline error distribution 
}

\references{
	J.K. Lee and M.O.Connell(2003). \emph{An S-Plus library for the analysis of differential expression}. In The Analysis of Gene Expression Data: Methods and Software. Edited by G. Parmigiani, ES Garrett, RA Irizarry ad SL Zegar. Springer, NewYork.
	
	Jain et. al. (2003) \emph{Local pooled error test for identifying
      differentially expressed genes with a small number of replicated microarrays}, Bioinformatics, 1945-1951.



}


\examples{
# Loading the library and the data
 library(LPE)
 data(Ley)
 
 dim(Ley)
 # Gives 12488*7 
 # First column is ID.


 # Subsetting the data
 subset.Ley <- Ley[1:1000,]
  
  subset.Ley[,2:7] <- preprocess(subset.Ley[,2:7],data.type="MAS5")
  # preprocess the data
  
 # Finding the baseline distribution of condition 1 and 2.
 var.1 <- baseOlig.error(subset.Ley[,2:4], q=0.01)
  median.x <- apply(subset.Ley[,2:4], 1, median)

 sf.x <- smooth.spline(var.1[, 1], var.1[, 2], df = 10)
  
 var.test <- fixbounds.predict.smooth.spline(sf.x, median.x)$y

}

\keyword{methods} %from KEYWORD.db
