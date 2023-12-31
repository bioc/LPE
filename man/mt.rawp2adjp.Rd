\name{mt.rawp2adjp.LPE}
\alias{mt.rawp2adjp.LPE}
\title{Adjusted p-values for simple multiple testing procedures}
\description{
  This function computes adjusted \eqn{p}-values for simple
  multiple testing procedures from a vector of raw (unadjusted)
  \eqn{p}-values. The procedures include the Bonferroni, Holm (1979),
  Hochberg (1988), and Sidak procedures for strong control of the
  family-wise Type I error rate (FWER), and the Benjamini & Hochberg
  (1995) and Benjamini & Yekutieli (2001) procedures for (strong)
  control of the false discovery rate (FDR). 
}
\usage{
mt.rawp2adjp.LPE(rawp, proc=c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY"))
}
\arguments{
  \item{rawp}{A vector of raw (unadjusted) \eqn{p}-values for each
    hypothesis under consideration. These could be nominal
    \eqn{p}-values, for example, from t-tables, or permutation
    \eqn{p}-values as given in \code{mt.maxT} and \code{mt.minP}. If the
    \code{mt.maxT} or \code{mt.minP} functions are used, raw
    \eqn{p}-values should be given in the original data order,
    \code{rawp[order(index)]}.} 
  \item{proc}{A vector of character strings containing the names of the
    multiple testing procedures for which adjusted \eqn{p}-values are to
    be computed. This vector should include any of the following:
    \code{"Bonferroni"}, \code{"Holm"}, \code{"Hochberg"},
    \code{"SidakSS"}, \code{"SidakSD"}, \code{"BH"}, \code{"BY"}. 
}
}
\details{
  Adjusted \eqn{p}-values are computed for simple FWER and FDR
  controlling procedures based on a vector of raw (unadjusted)
  \eqn{p}-values.
\describe{
   \item{Bonferroni}{Bonferroni single-step adjusted \eqn{p}-values
   for strong control of the FWER.}
   \item{Holm}{Holm (1979) step-down adjusted \eqn{p}-values for
   strong control of the FWER.}
   \item{Hochberg}{Hochberg (1988) step-up adjusted \eqn{p}-values
   for 
   strong control of the FWER (for raw (unadjusted) \eqn{p}-values
   satisfying the Simes inequality).}
   \item{SidakSS}{Sidak single-step adjusted \eqn{p}-values for
   strong control of the FWER (for positive orthant dependent test
   statistics).}
   \item{SidakSD}{Sidak step-down adjusted \eqn{p}-values for
   strong control of the FWER (for positive orthant dependent test
   statistics).}
   \item{BH}{adjusted \eqn{p}-values for the Benjamini & Hochberg
   (1995) step-up FDR controlling procedure (independent and positive
   regression dependent test statistics).}
   \item{BY}{adjusted \eqn{p}-values for the Benjamini & Yekutieli
   (2001) step-up FDR controlling procedure (general dependency
   structures).} 
 }
}

\value{
A list with components
  \item{adjp}{A matrix of adjusted \eqn{p}-values, with rows
    corresponding to hypotheses and columns to multiple testing
    procedures. Hypotheses are sorted in increasing order of their raw
    (unadjusted) \eqn{p}-values.} 
  \item{index}{A vector of row indices, between 1 and
    \code{length(rawp)}, where rows are sorted according to 
    their raw (unadjusted) \eqn{p}-values. To obtain the adjusted
    \eqn{p}-values in the original data order, use
    \code{adjp[order(index),]}. 
}
}
\references{
Y. Benjamini and Y. Hochberg (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple
testing. \emph{J. R. Statist. Soc. B}. Vol. 57: 289-300.\cr 

Y. Benjamini and D. Yekutieli (2001). The control of the false discovery
rate in multiple hypothesis testing under dependency. \emph{Annals of
  Statistics}. Accepted.\cr 

S. Dudoit, J. P. Shaffer, and J. C. Boldrick (Submitted). Multiple
hypothesis testing in microarray experiments.\cr 

Y. Ge, S. Dudoit, and T. P. Speed (In preparation). Fast algorithm for 
resampling-based \eqn{p}-value adjustment in multiple testing. \cr

Y. Hochberg (1988). A sharper Bonferroni procedure for multiple tests of
significance, \emph{Biometrika}. Vol. 75: 800-802.\cr 

S. Holm (1979). A simple sequentially rejective multiple test
procedure. \emph{Scand. J. Statist.}. Vol. 6: 65-70. 
}
 
\author{
  Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine},\cr
  Yongchao Ge, \email{gyc@stat.berkeley.edu}. 
}

\seealso{\code{\link{lpe}}}

%\seealso{\code{\link{mt.maxT}}, \code{\link{mt.minP}},
%  \code{\link{mt.plot}}, \code{\link{mt.reject}}, \code{\link{golub}}.}

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
 fdr.BH <- fdr.adjust(lpe.result, adjp="BH")
}

\keyword{htest}





