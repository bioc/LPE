\name{Ley}
\alias{Ley}

\title{Gene Expression Data from Mouse Immune response study, (2002)}

\usage{
  data(Ley)
}

\description{
  Affymetrix GeneChip (12488 genes and 3 different conditions, each with 3 replicates) experiment was conducted by Dr. Ley (2002) to understand mouse immune response study. For demonstration purposes, we show data from 2 conditions, i.e. 
6 chips, only.
}

\format{
 \item{Ley}{matrix of 12488 genes by 7 columns. First column is the GeneID, 
	  next three columns are replicates in condition one, and 
last three columns are replicates in condition 2.
 For details, contact. Dr. Klaus Ley (kfl3f@virginia.edu) 
 http://hsc.virginia.edu/medicine/basic-sci/biomed/ley/}
}
\author{ Nitin Jain\email{nitin.jain@pfizer.com} }


\references{
	J.K. Lee and M.O.Connell(2003). \emph{An S-Plus library for the analysis of differential expression}. In The Analysis of Gene Expression Data: Methods and Software. Edited by G. Parmigiani, ES Garrett, RA Irizarry ad SL Zegar. Springer, NewYork.
	
	Jain et. al. (2003) \emph{Local pooled error test for identifying
      differentially expressed genes with a small number of replicated microarrays}. (To appear in Bioinformatics)
	See \url{http://hesweb1.med.virginia.edu/bioinformatics/research/expression/index.html}.
}

\keyword{datasets}
