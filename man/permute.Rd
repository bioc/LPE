\name{permute}
\alias{permute}

\title{
	Calculating all possible permutations of a vector
}
\description{ 
	Given a vector, all possible combinations of vector are obtained 
}
\usage{
	permute(a)
}
\arguments{
 	\item{a}{a is any numeric vector.}
}

\value{
	A vector containing the possible combinations.
}

\details{
 Used in am.trans. Does all permutations for columns within an experimental 
 condition so that A and M can be calculated for all permutations of chips within a treatment.
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
  # Loading LPE library
  library(LPE)
  
  # A test vector
  permute(1:3)  
  
 # Returns a 2 by 3 matrix 
 #       [,1] [,2] [,3]
 #  [1,]   2    3    1
 #  [2,]   3    1    2

}
\keyword{methods} %from KEYWORD.db
