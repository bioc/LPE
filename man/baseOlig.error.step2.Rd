\name{baseOlig.error.step2}
\alias{baseOlig.error.step2}

\title{
 Evaluates LPE variance function of M for quantiles of A within and 
 experimental condition. It is based on the adaptive number of intervals.
}

\description{
 Similar to baseOlig.error.step1 function, except that now the number of bins
 are chosen adaptively instead of fixed 100.
 }

\usage{
  baseOlig.error.step2(y,baseOlig.error.step1.res, df=10, stats=median, min.genes.int=10, div.factor=1)
}

\arguments{
 \item{y}{y is a preprocessed matrix or data frame of expression 
      intensities in which columns are expression intensities for
      a particular experimental condition and rows are genes.}
 \item{baseOlig.error.step1.res}{It is the result obtained from 
 	baseOlig.error.step1 function, in which number of bins are fixed=100}
 \item{df}{df stands for degrees of freedom. It is used in 
 	   smooth.spline function to interpolate the variances
	   of all genes. Default value is 10.}
 \item{stats}{It determines whether mean or median is to be used for the replicates}
  \item{min.genes.int}{Determines the minimum number of genes in a subinterval for selecting the adaptive intervals.}
 \item{div.factor}{Determines the factor by which sigma needs to be divided for
 selecting adaptive intervals.}

}
 
\value{
 Returns object of class baseOlig comprising a data frame with 2 columns: A
  and var M, and rows for each quantile specified. The A column contains
  the median values of A for each quantile/bin and the M columns contains 
  the pooled variance of the replicate chips for genes within each quantile/bin.
}


\references{
	J.K. Lee and M.O.Connell(2003). \emph{An S-Plus library for the analysis of differential expression}. In The Analysis of Gene Expression Data: Methods and Software. Edited by G. Parmigiani, ES Garrett, RA Irizarry ad SL Zegar. Springer, NewYork.

	Jain et. al. (2003) \emph{Local pooled error test for identifying
      differentially expressed genes with a small number of replicated microarrays}, Bioinformatics, 1945-1951.

}


\seealso{
  \code{\link{lpe}}
}

\examples{

  # Loading the library and the data
  library(LPE)
  data(Ley)
  
  dim(Ley)
  # Gives 12488 by 7
  Ley[1:3,]
   # Returns 
  #       ID           c1   c2   c3    t1    t2    t3
#   1  AFFX-MurIL2_at 4.06 3.82 4.28 11.47 11.54 11.34
#   2 AFFX-MurIL10_at 4.56 2.79 4.83  4.25  3.72  2.94
#   3  AFFX-MurIL4_at 5.14 4.10 4.59  4.67  4.71  4.67

  Ley[1:1000,2:7] <- preprocess(Ley[1:1000,2:7],data.type="MAS5")
  # Finding the baseline distribution of subset of the data
  # condition one (3 replicates)
  var.1 <- baseOlig.error.step1(Ley[1:1000,2:4], q=0.01, df=10)
  dim(var.1)
  var.11 <- baseOlig.error.step2(Ley[1:1000,2:4], var.1, df=10)
  # Returns a matrix of 1000 by 2 (A,M) format
}
\keyword{methods} % from KEYWORDS.db
