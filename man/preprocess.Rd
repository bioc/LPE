\name{preprocess}
\alias{preprocess}

\title{
  Preprocessing the data (IQR normalization, thresholding, and log 
  transformation)
}

\description{
  Finds inter-quartile range of the data = (75th percentile - 25th percentile),
  thresholds low intensity MAS4 and dChip data to 1, MAS5 data to 0.1 and does log based 2 transformation
  
}

\usage{
  preprocess(x, data.type="MAS5")
}

\arguments{
  \item{x}{x is the data-set which needs preprocessing.}
  \item{data.type}{Three types of data accepted in the current version :
          MAS4 (Microarray suite software) , MAS5 and dChip}
}

\value{
  Returns a data-set of same dimensions as that of the input data. It has 
  IQR normalization for MAS4 and MAS5 data. Low intensities of MAS4 and dChip data are thresholded to 1, and those of MAS5 to 0.1. Finally dat is transformed to base 2.
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
  # Loading the data set
  dim(Ley) #gives 12488 * 7
  Ley[1:3,]

  Ley[,2:7] <- preprocess(Ley[,2:7],data.type="MAS5")
  Ley[1:3,]
 
}

\keyword{methods} % from KEYWORDS.db
