\name{pre.hapassoc}
\alias{pre.hapassoc}
\title{Pre-process the data before fitting it with hapassoc} 
\description{
This function takes as an argument the original dataframe with
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
Initial estimates of haplotype frequencies are then obtained using the 
EM algorithm applied to the multilocus genotype data. 
Haplotypes with frequencies below a user-specified tolerance (zero.tol)
are assumed not to exist and are removed from further consideration.
(Pseudo-individuals having haplotypes of negligible frequency are deleted and 
the column in the design matrix corresponding to that haplotype is deleted.)
For the remaining haplotypes, those with non-negligible frequency below a 
user-defined pooling tolerance (pooling.tol) are pooled into a single 
category called \"pooled\" in the design matrix for the risk model. 
However, the frequencies of each of these pooled haplotypes are 
still calculated separately. 
}

\usage{pre.hapassoc(dat,numSNPs,maxMissingGenos=1,pooling.tol = 0.05, zero.tol = 1/(2 * nrow(dat) * 10)}
\arguments{
\item{dat}{the non-SNP and SNP data as a data frame. The SNP data should comprise the last 2*numSNPs columns.}
\item{numSNPs}{number of SNPs per haplotype}
\item{maxMissingGenos}{maximum number of single-locus genotypes with missing data to allow for each subject. (Subjects with more missing data, or with missing non-SNP data are removed.) The default is 1.}
\item{pooling.tol}{pooling tolerance -- by default set to 0.05}
\item{zero.tol}{tolerance for haplotype frequencies below which haplotypes
are assumed not to exist -- by default set to 
\eqn{\frac{1}{2*N*10}}{1/(2*N*10)} where N is the number of subjects}
}
\value{
\item{haplotest}{T/F, true if some haplotypes were pooled in the risk model}
\item{initGamma}{initial estimates of haplotype frequencies}
\item{zeroFreqHaplos}{list of haplos assumed not to exist}
\item{pooledHaplos}{list of haplos pooled into a single category in the design matrix}
\item{nonHaploDM}{non-haplotype portion of the AUGMENTED data frame}
\item{haploDM}{data frame with \eqn{2^{numSNPs}}{2^numSNPs} columns scoring number of copies of each haplotype for each pseudo-individual}
\item{haploMat}{matrix with 2 columns giving haplotypes for each pseudo-individual}
\item{wt}{vector giving initial weights for each pseudo-individual for 
the EM algorithm}
\item{ID}{index for each individual in the original data frame. Note that all pseudo-individuals have the same ID value}
\item{unknown}{vector indicating whether the haplotype information was missing for each row in the augmented data}
}
\examples{
data(hypoDat)
example.pre.hapassoc<-pre.hapassoc(hypoDat, numSNPs=3)

# To get the initial haplotype frequencies:
example.pre.hapassoc$initGamma
#      h000       h001       h010       h011       h100       h101       h110 
#0.25179111 0.26050418 0.23606001 0.09164470 0.10133627 0.02636844 0.01081260 
#      h111 
#0.02148268 
# The '001' haplotype is estimated to be the most frequent

example.pre.hapassoc$pooledHaplos
# "h101" "h110" "h111"
# These haplotypes are to be pooled in the design matrix for the risk model

names(example.pre.hapassoc$haploDM)
# "h000"   "h001"   "h010"   "h011"   "h100"   "pooled"
}
\references{Burkett K, McNeney B, Graham J (2004). 
A note on inference of trait associations with SNP
haplotypes and other attributes in generalized linear models.
Human Heredity, In press}}
\seealso{
\code{\link{hapassoc}},\code{\link{summary.hapassoc}}.
}
\keyword{methods}
