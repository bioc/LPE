\name{resamp.adj}
\alias{resamp.adj}

\title{
   Resampling based fdr adjustment
}

\description{
  Adjusts the fdr based on rank invariant genes
}

\usage{
 resamp.adj(x,y, q=0.01, iterations=5, min.genes.int=10) 
}

\arguments{
 \item{x}{Replicated data from first experimental condition (as matrix 
			or data-frame)}.
 \item{y}{Replicated data from second experimental condition (as matrix 
			or data-frame)}.
 \item{q}{q is the quantile width; q=0.01 corresponds to 100 quantiles}.
 \item{iterations}{Number of iterations to be performed to obtain critical z-statistics}.
 \item{min.genes.int}{Determines the minimum number of genes in a subinterval for selecting the adaptive intervals.}
}

\details{
 Returns the z-statistics for the null distribution, obtained from resampling
 the rank invariant genes within each quantile. These z-statistic values are
 compared with z-statiscs from the original data, and fdr is calculated.
}


\references{
	J.K. Lee and M.O.Connell(2003). \emph{An S-Plus library for the analysis of differential expression}. In The Analysis of Gene Expression Data: Methods and Software. Edited by G. Parmigiani, ES Garrett, RA Irizarry ad SL Zegar. Springer, NewYork.
	
	Jain et. al. (2003) \emph{Local pooled error test for identifying
      differentially expressed genes with a small number of replicated microarrays}, Bioinformatics, 1945-1951.
      
      }


\examples{
 
  # Loading the library and the data
 library(LPE)
 data(Ley)
 
 dim(Ley)
 # Gives 12488*7 
 # First column is ID.

 # Subsetting the data
 subset.Ley <- Ley[1:1000,]
  
  subset.Ley[,2:7] <- preprocess(subset.Ley[,2:7],data.type="MAS5")
   
 # Finding the baseline distribution of condition 1 and 2.
 var.1 <- baseOlig.error(subset.Ley[,2:4], q=0.01)
 var.2 <- baseOlig.error(subset.Ley[,5:7], q=0.01)
 
 # Applying LPE
 lpe.result <- lpe(subset.Ley[,2:4],subset.Ley[,5:7], var.1, var.2,
		probe.set.name=subset.Ley[,1])
  

 
 z.stats.null <- resamp.adj(subset.Ley[,2:4], subset.Ley[,5:7], q=0.01, iterations=2,min.genes.int=10 )

}

\keyword{methods} % from KEYWORDS.db


