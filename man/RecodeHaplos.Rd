\name{RecodeHaplos}
\alias{RecodeHaplos}
\title{Recode SNP genotype data into possible haplotypes}
\description{
This function takes as an argument the original data frame with 
non-SNP and SNP data and converts the genotype data at single SNPs
(the single-locus genotypes) into haplotype data. 
The rows of the original data frame should correspond to
subjects and each SNP should have two columns, one for each
allele of the single-locus genotype. The SNP data should comprise the 
last 2*numSNPs columns.  If the haplotypes for a subject cannot be 
inferred from his or her genotype data, "pseudo-individuals" 
representing all possible haplotype combinations consistent with 
the single-locus genotypes are considered.
Missing single-locus genotypes, up to a maximum of maxMissingGenos (see
below), are allowed, but subjects with missing data in more than 
maxMissingGenos, or with missing non-SNP data, are removed.
}
\synopsis{RecodeHaplos(dat, numSNPs, maxMissingGenos = 1, logriskmodel="additive")}
\usage{RecodeHaplos(dat, numSNPs, maxMissingGenos = 1)}
\arguments{
\item{dat}{the non-SNP and SNP data as a data frame}
\item{numSNPs}{the number of SNPs in a haplotype}
\item{maxMissingGenos}{maximum number of single-locus genotypes with missing data to allow for each subject. (Subjects with more missing data, or with missing non-SNP data are removed.) The default is 1.}
}
\value{
\item{nonHaploDM}{non-genetic portion of an AUGMENTED data frame which includes
any possible pseudo-individuals}
\item{haploDM}{a data frame with \eqn{2^{numSNPs}}{2^numSNPs} columns scoring number of copies of each haplotype for each pseudo-individual}
\item{haploMat}{matrix with 2 columns giving the haplotypes as a string of 
zeros and ones for each pseudo-individual}
\item{wt}{vector of crude weights for each pseudo-individual to 
be updated by CheckHaplos}
\item{ID}{index for each individual in the original data frame. Note that 
all pseudo-individuals for a subject have the same ID value}
}
\examples{
# To convert SNP data from a simulated diabetes example into two-locus haplotypes:
data(hypoDat)
example.haplos<-RecodeHaplos(hypoDat, 2)
}
\seealso{
\code{\link{CheckHaplos}},\code{\link{EM}},\code{\link{summary.EM}}.
}
\keyword{methods}
