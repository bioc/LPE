\name{quan.norm}
\alias{quan.norm}

\title{
  Finding quartile range
}

\description{
  Finds quartile range of the data (default is IQR = 75th percentile - 25th percentile).
}

\usage{
  quan.norm(x, percent=50)
}

\arguments{
  \item{x}{x is a vector for which quartile range has to be found.}
  \item{percent}{Percentage for which quartile range is needed}
}

\value{
  Returns a numeric value representing the difference of 75th percentile 
  and 25th percentile of the vector. It is used for normalization across 
  the chips - basic assumption is that net differential expression of 
  the middle half of the genes in microarray experiment is zero, which 
  is conservative assumption as typically only 5-10% genes show 
  differential expression.
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
  library(LPE)
  # Loading the LPE library
 
  quan.norm(1:5) 
  % returns 2 [= 75th percentile(4)- 25th percentile(2)]
}

\keyword{methods} % from KEYWORDS.db
