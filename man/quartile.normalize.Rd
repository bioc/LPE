\name{quartile.normalize}
\alias{quartile.normalize}

\title{
  Normalization based on quartile range
}

\description{
Does Normalization based on quartile range 
}

\usage{
  quartile.normalize(x, percent=50)
}

\arguments{
  \item{x}{x is a matrix or dataframe on which normalization has to be performed.}
  \item{percent}{Percentage for which normalization is needed}
}

\value{
Returns the normalized data based on quartile normalization 
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

  data(Ley) 
 
 dim(Ley)
 # Gives 12488*7
 # First column is ID.

  subset <- 1:1000
  Ley[subset,2:7] <- quartile.normalize(Ley[subset,2:7],percent=50)


}

\keyword{methods} % from KEYWORDS.db
