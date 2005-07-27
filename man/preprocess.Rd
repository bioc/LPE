\name{preprocess}
\alias{preprocess}

\title{
  Preprocessing the data (IQR normalization, thresholding, log-
  transformation, and lowess normalization)
}

\description{
  Finds inter-quartile range of the data = (75th percentile - 25th percentile),
  thresholds low intensity MAS4, MAS5 and dChip data to 1, then log transforms
  the data (base 2), and dows lowess normalization  
}

\usage{
  preprocess(x, data.type="MAS5",threshold=1,LOWESS=FALSE)
}

\arguments{
  \item{x}{x is the data-set which needs preprocessing.}
  \item{data.type}{Three types of data accepted in the current version :
          MAS4 (Microarray suite software) , MAS5 and dChip}
  \item{threshold}{threshold is the 'thresholding value' below which
                  all data would be thresholded (default = 1).}
  \item{LOWESS}{LOWESS is a logical variable which determines if lowess 
  	normalization needs to be performed.}

}

\value{
  Returns a data-set of same dimensions as that of the input data. It has 
  IQR normalization for MAS4 and MAS5 data. Low intensities of MAS4, MAS5
  and dChip data are thresholded to 1. Then data is transformed to base 2. If
LOWESS normalization parameter is set as TRUE, then lowess  
normalization is performed.
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

  Ley[1:1000,2:7] <- preprocess(Ley[1:1000,2:7],data.type="MAS5",
  	threshold=1, LOWESS=TRUE)
  Ley[1:3,]
 
}

\keyword{methods} % from KEYWORDS.db
