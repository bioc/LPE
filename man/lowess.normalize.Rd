\name{lowess.normalize}
\alias{lowess.normalize}

\title{
  lowess normalization of the data (based on M vs A graph)
}

\description{
 All the chips are normalized w.r.t. 1st chip 
}

\usage{
  lowess.normalize(x,y)
}

\arguments{
  \item{x}{x is the chip data w.r.t. which other chips would be normalized}
  \item{y}{y is the chip data which would be normalized}
}

\value{
  Returns the lowess normalized chip intensity.
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
 
  data(Ley)
  # Loading the data set
  dim(Ley) #gives 12488 * 7
  Ley[1:3,]

  Ley[1:1000,2:7] <- preprocess(Ley[1:1000,2:7],data.type="MAS5")
  Ley[1:3,]
 
}

\keyword{methods} % from KEYWORDS.db
