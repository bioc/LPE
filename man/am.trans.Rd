\name{am.trans}
\alias{am.trans}

\title{
  Transform replicated arrays into (A,M) format
}

\description{
  Transforms expression intensity measurements for replicated arrays of a single experimental condition into (A,M) format: A = (xi+xj)/2, M = (xi-xj) where x1, x2,..., xn are individual chips. This function is used in the estimation of within-bin variances in the LPE function, and not typically by the user.
}

\usage{
  am.trans(y)
}

\arguments{
  \item{y}{y is an ngene by n matrix of expression intensity measurements for
     replicated arrays under a single experimental condition.}
}

\value{
  Returns matrix with 2 columns cbind(A.M) and rows comprising all permutations of individual chip columns of the input matrix y. Note that for each pair of chips M is calculated twice, once for xi-xj and once for xj-xi. The resulting marix thus has twice the number of rows as the input matrix y.  
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
 
  # Test data with duplicates
  temp1 <- matrix(c(1,20,1.5,23),nrow=2)
  am.trans(temp1)
  # It gives a matrix of (4*2) as only two permutaions
  # are possible for each row with duplicates (X1-X2, and X2-X1)
  
  
  # Another test data with three replicates
  temp2 <- matrix(c(1,20,1.5,23,0.8,19),nrow=2)
  am.trans(temp2)
  # Now it returns matrix of (12*2) as there are
  # 6 possible permutations for each row with triplicates
  # (X1-X2, X1-X3, X2-X3, X2-X1, X3-X1 and X3-X2)
  
}

\keyword{methods} % from KEYWORDS.db
