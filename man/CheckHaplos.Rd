\name{CheckHaplos}
\alias{CheckHaplos}
\title{Estimation of Haplotype Frequencies}
\description{
This function estimates haplotype frequencies using the EM algorithm
applied to the multilocus genotype data output by RecodeHaplos. 
Haplotypes with frequencies below a user-specified tolerance (zero.tol)
are assumed to not exist and are removed from further consideration.
(Pseudo-individuals having haplotypes with zero frequency are deleted and 
the column corresponding to that haplotype is deleted.)
From the remaining haplotypes, those with frequency below a user-defined
pooling tolerance are pooled into a single category called \"pooled\". 
}
\usage{CheckHaplos(haplos.list,numSNPs,pooling.tol = 0.05, zero.tol = 1/(2 * sum(haplos.list$wt) * 10)}
\arguments{
\item{haplo.list}{list output by the function RecodeHaplos}
\item{numSNPs}{number of SNPs per haplotype}
\item{pooling.tol}{pooling tolerance -- set to 0.05 by default}
\item{zero.tol}{tolerance for haplotype frequencies below which haplotypes
are assumed not to exist -- set to \eqn{\frac{1}{2*N*10}}{1/(2*N*10)} 
where N is the number of subjects by default}
}
\value{
\item{haplotest}{T/F, True if some pooling of the haplotypes was done}
\item{initGamma}{initial estimates of haplotype frequencies}
\item{zeroFreqHaplos}{list of haplos assumed not to exist}
\item{pooledHaplos}{list of haplos pooled into a single category}
\item{nonHaploDM}{non-genetic portion of the AUGMENTED data frame}
\item{haploDM}{a data frame with \eqn{2^{numSNPs}}{2^numSNPs} columns scoring number of copies of each haplotype for each pseudo-individual}
\item{haploMat}{matrix with 2 columns giving haplotypes for each pseudo-individual}
\item{wt}{vector giving initial weights for each pseudo-individual for 
the EM algorithm}
\item{ID}{index for each individual in the original data frame. Note that all pseudo-individuals have the same ID value}
\item{unknown}{vector indicating whether the haplotype information was missing for each row in the augmented data}
}
\examples{
data(hypoDat)
example.haplos<-RecodeHaplos(hypoDat,2)
example.newhaplos<-CheckHaplos(example.haplos, 2)

# To get the haplotype counts

example.newhaplos$initGamma

# Result:
#       n00         n01         n10         n11 
#0.009465535 0.246534465 0.437534465 0.306465535 
# The '10' haplotype is the most frequent
}
\seealso{
\code{\link{RecodeHaplos}},\code{\link{EM}},\code{\link{summary.EM}}.
}
\keyword{methods}
