\name{baseOlig.error}
\alias{baseOlig.error}

\title{
 Evaluates LPE variance function of M for quantiles of A within and 
 experimental condition and then interpolates it for all genes.
}

\description{
 Calls baseOlig.error.step1 and baseOlig.error.step2 functions in order to
 calculate the baseline distribution.
}

\usage{
  baseOlig.error(y, stats=median, q=0.01, min.genes.int=10,div.factor=1)
}

\arguments{
 \item{y}{y is a preprocessed matrix or data frame of expression intensities in which columns are expression intensities for a particular experimental
          condition and rows are genes.}
 \item{stats}{It determines whether mean or median is to be used for the replicates}
 \item{q}{q is the quantile width; q=0.01 corresponds to 100 quantiles
 	  i.e. percentiles. Bins/quantiles have equal number of genes and
 	  are split according to the average intensity A.}
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

  Ley[,2:7] <- preprocess(Ley[,2:7],data.type="MAS5")
  
  subset <- 1:1000
  Ley.subset <- Ley[subset,]
  
  # Finding the baseline distribution of subset of the data
  # condition one (3 replicates)
  var.1 <- baseOlig.error(Ley.subset[,2:4], q=0.01)
  dim(var.1)
  # Returns a matrix of 1000 by 2 (A,M) format
}
\keyword{methods} % from KEYWORDS.db
