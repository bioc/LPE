\name{fdr.adjust}
\alias{fdr.adjust}

\title{
  FDR adjustment procedures
}

\description{
  Based on the type of adjustment, eg: resampling, BH, BY, etc, 
  calls appropriate functions for fdr adjustment
}

\usage{
 fdr.adjust(lpe.result,adjp="resamp",target.fdr=c(10^-3 ,seq(0.01,0.10,0.01), 0.15, 0.20, 0.50),iterations=5,ALL=FALSE )
}

\arguments{
 \item{lpe.result}{Data frame obtained from calling lpe function}
 \item{adjp}{Type of adjustment procedure}
 \item{target.fdr}{Desired FDR level (used only for resampling
 	based adjustment)}
 \item{iterations}{Number of iterations for stable z-critical.}
 \item{ALL}{If TRUE, the FDR corresponding to all the z-statistics, i.e.
 	for every gene intensity is given.}
}


\details{
 Returns the output similar to lpe function, including adjusted fdr
 }
\author{ Nitin Jain\email{nitin.jain@pfizer.com} }

\references{
	J.K. Lee and M.O.Connell(2003). \emph{An S-Plus library for the analysis of differential expression}. In The Analysis of Gene Expression Data: Methods and Software. Edited by G. Parmigiani, ES Garrett, RA Irizarry ad SL Zegar. Springer, NewYork.
	

	Jain et. al. (2003) \emph{Local pooled error test for identifying
      differentially expressed genes with a small number of replicated microarrays}, Bioinformatics, 1945-1951.

Jain et. al. (2005) \emph{Rank-invariant resampling based estimation of false discovery rate for analysis of small sample microarray data}, BMC Bioinformatics, Vol 6, 187.

	}

\examples{

 # Loading the library and the data
 library(LPE)
 data(Ley)
 
 dim(Ley)
 # Gives 12488*7 
 # First column is ID.

 Ley[,2:7] <- preprocess(Ley[,2:7],data.type="MAS5")

 # Subsetting the data
 subset.Ley <- Ley[1:1000,]
  
   
 # Finding the baseline distribution of condition 1 and 2.
 var.1 <- baseOlig.error(subset.Ley[,2:4], q=0.01)
 var.2 <- baseOlig.error(subset.Ley[,5:7], q=0.01)

 # Applying LPE
 lpe.result <- lpe(subset.Ley[,2:4],subset.Ley[,5:7], var.1, var.2,
                probe.set.name=subset.Ley[,1])


 final.result <- fdr.adjust(lpe.result, adjp="resamp", target.fdr=c(0.01,0.05), iterations=1)
 final.result
  
}

\keyword{methods} % from KEYWORDS.db



 


  
