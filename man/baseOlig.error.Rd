\name{baseOlig.error}
\alias{baseOlig.error}

\title{
 Evaluates LPE variance function of M for quantiles of A within and experimental condition
}

\description{
  Genes are placed in bins/quantiles according to their average expression intensity. The function baseOlig.error calculates a pooled variance of M for genes within these bins/quantiles of A for the replicates of the experimental condition contained in y.
}

\usage{
  baseOlig.error(y, q=0.01)
}

\arguments{
 \item{y}{y is a preprocessed matrix or data frame of expression intensities in which columns are expression intensities for a particular experimental
          condition and rows are genes.}
 \item{q}{q is the quantile width; q=0.01 corresponds to 100 quantiles
 	  i.e. percentiles. Bins/quantiles have equal number of genes and
 	  are split according to the average intensity A.}
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
  # Finding the baseline distribution of subset of the data
  # condition one (3 replicates)
  var.1 <- baseOlig.error(Ley[1:1000,2:4], q=0.01)
  dim(var.1)
  # Returns a matrix of 100 by 2 (A,M) format
}
\keyword{methods} % from KEYWORDS.db
