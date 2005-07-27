\name{n.genes.adaptive.int}
\alias{n.genes.adaptive.int}

\title{
 Calcuates the number of genes in various intervals adaptively.
}

\description{
  Instead of dividing the genes equally in 100 intervals, this
  function divides them adaptively based on three rules:
  a) min. number of genes (default =10), b) max. number of genes = total/100;
  c) based on Median + fraction(SD) from the starting gene of each
  interval
}

\usage{
  n.genes.adaptive.int(baseOlig.error.step1.res,
  		min.genes.int=10, div.factor=1)
}

\arguments{
 \item{baseOlig.error.step1.res}{It is the result 
 		from baseOlig.error.step1 function.}
 \item{min.genes.int}{It is the minimum number of genes in the interval,
 	default=10.}
 \item{div.factor}{(1/div.factor) is the fraction of Standard Deviation
 		which we wish to include in each interval to calculate
		number of genes in each interval}
}
 
\value{
  Returns a vector respresenting the number of genes in each interval.
}
\author{ Nitin Jain\email{nitin.jain@pfizer.com} }


\references{
	J.K. Lee and M.O.Connell(2003). \emph{An S-Plus library for the analysis of differential expression}. In The Analysis of Gene Expression Data: Methods and Software. Edited by G. Parmigiani, ES Garrett, RA Irizarry ad SL Zegar. Springer, NewYork.

	Jain et. al. (2003) \emph{Local pooled error test for identifying
      differentially expressed genes with a small number of replicated microarrays}, Bioinformatics, 1945-1951.

Jain et. al. (2005) \emph{Rank-invariant resampling based estimation of false discovery rate for analysis of small sample microarray data}, BMC Bioinformatics, Vol 6, 187.

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
  var.1 <- baseOlig.error.step1(Ley[1:1000,2:4], q=0.01)
  dim(var.1)
  # Returns a matrix of 1000 by 2 (A,M) format
  n.genes.subint <- n.genes.adaptive.int(var.1, min.genes.int=10, div.factor=1)
}
\keyword{methods} % from KEYWORDS.db
