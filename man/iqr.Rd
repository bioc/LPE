\name{iqr}
\alias{iqr}

\title{
  Inter-quartile range
}

\description{
  Finds inter-quartile range of the data = {75th percentile - 25th percentile}.
}

\usage{
  iqr(x)
}

\arguments{
  \item{x}{x is a vector for which IQR has to be found.}
}

\value{
  Returns a numeric value representing the difference of 75th percentile 
  and 25th percentile of the vector. It is used for normalization across 
  the chips - basic assumption is that net differential expression of 
  the middle half of the genes in microarray experiment is zero, which 
  is conservative assumption as typically only 5-10% genes show 
  differential expression.
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
  library(LPE)
  # Loading the LPE library
 
  iqr(1:5) 
  % returns 2 [= 75th percentile(4)- 25th percentile(2)]
}

\keyword{methods} % from KEYWORDS.db
