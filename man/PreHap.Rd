\name{pre.hapassoc}
\alias{pre.hapassoc}
\title{Pre-process the data before fitting it with hapassoc} 
\description{
This function takes as an argument a dataframe with
non-SNP and SNP data and converts the genotype data at single SNPs
(the single-locus genotypes) into haplotype data.
The rows of the input data frame should correspond to
subjects. Single-locus SNP genotypes  may be specified in 
one of two ways: (i) as pairs of columns, with one column for
each allele of the single-locus genotypes (\dQuote{allelic format}),
or (ii) as columns of two-character genotypes (\dQuote{genotypic format}).
The SNP data should comprise the
last 2*numSNPs columns (allelic format) or the last numSNPs columns (genotypic
format) of the data frame. 

If the haplotypes for a subject cannot be
inferred from his or her genotype data, ``pseudo-individuals''
representing all possible haplotype combinations consistent with
the single-locus genotypes are considered.
Missing single-locus genotypes, up to a maximum of maxMissingGenos (see
below), are allowed, but subjects with missing data in more than
maxMissingGenos, or with missing non-SNP data, are removed.
Initial estimates of haplotype frequencies are then obtained using the 
EM algorithm applied to the genotype data. 
Haplotypes with frequencies below a user-specified tolerance (zero.tol)
are assumed not to exist and are removed from further consideration.
(Pseudo-individuals having haplotypes of negligible frequency are deleted and 
the column in the design matrix corresponding to that haplotype is deleted.)
For the remaining haplotypes, those with non-negligible frequency below a 
user-defined pooling tolerance (pooling.tol) are pooled into a single 
category called \dQuote{pooled} in the design matrix for the risk model. 
However, the frequencies of each of these pooled haplotypes are 
still calculated separately. 
}

\usage{pre.hapassoc(dat,numSNPs,maxMissingGenos=1,pooling.tol = 0.05, 
       zero.tol = 1/(2 * nrow(dat) * 10), allelic=TRUE, verbose=TRUE)}
\arguments{
\item{dat}{the non-SNP and SNP data as a data frame. The SNP data should comprise the last 2*numSNPs columns (allelic format) or last numSNPs columns (genotypic format).
Missing allelic data should be coded as \code{NA} or \code{""} and missing genotypic data should be coded as, e.g., \code{"A"} if one allele is missing and \code{""} if both alleles are missing.  }
\item{numSNPs}{number of SNPs per haplotype}
\item{maxMissingGenos}{maximum number of single-locus genotypes with missing data to allow for each subject. (Subjects with more missing data, or with missing non-SNP data are removed.) The default is 1.}
\item{pooling.tol}{pooling tolerance -- by default set to 0.05}
\item{zero.tol}{tolerance for haplotype frequencies below which haplotypes
are assumed not to exist -- by default set to 
\eqn{\frac{1}{2*N*10}}{1/(2*N*10)} where N is the number of subjects}
\item{allelic}{TRUE if single-locus SNP genotypes are in allelic format and FALSE if in genotypic format; default is TRUE.}
\item{verbose}{indicates whether or not a list of the genotype
variables used to form haplotypes and a list of other non-genetic
variables should be printed; default is TRUE.}
}
\value{
\item{haplotest}{logical, TRUE if some haplotypes were pooled in the risk model}
\item{initFreq}{initial estimates of haplotype frequencies}
\item{zeroFreqHaplos}{list of haplotypes assumed not to exist}
\item{pooledHaplos}{list of haplotypes pooled into a single category in the design matrix}
\item{haploDM}{Haplotype portion of the data frame \strong{augmented} with pseudo-individuals. Has \eqn{2^{numSNPs}}{2^numSNPs} columns scoring number of copies of each haplotype for each pseudo-individual}
\item{nonHaploDM}{non-haplotype portion of the data frame \strong{augmented} with pseudo-individuals}
\item{haploMat}{matrix with 2 columns listing haplotype labels for each pseudo-individual}
\item{wt}{vector giving initial weights for each pseudo-individual for 
the EM algorithm}
\item{ID}{index for each individual in the original data frame. Note that all pseudo-individuals have the same ID value}
}
\examples{
#First example data set has single-locus genotypes in "allelic format"
data(hypoDat)
example.pre.hapassoc<-pre.hapassoc(hypoDat, numSNPs=3)

# To get the initial haplotype frequencies:
example.pre.hapassoc$initFreq
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

####
#Second example data set has single-locus genotypes in "genotypic format"
data(hypoDatGeno)
example2.pre.hapassoc<-pre.hapassoc(hypoDatGeno, numSNPs=3, allelic=FALSE)

# To get the initial haplotype frequencies:
example2.pre.hapassoc$initFreq
#      hAAA       hAAC       hACA       hACC       hCAA       hCAC
#0.25179111 0.26050418 0.23606001 0.09164470 0.10133627 0.02636844
#      hCCA           hCCC 
#0.01081260	0.02148268 
# The 'hAAC' haplotype is estimated to be the most frequent

example2.pre.hapassoc$pooledHaplos
#  "hCAC" "hCCA" "hCCC"
# These haplotypes are to be pooled in the design matrix for the risk model

names(example2.pre.hapassoc$haploDM)
#  "hAAA"   "hAAC"   "hACA"   "hACC"   "hCAA"   "pooled"
}
\references{Burkett K, McNeney B, Graham J (2004). 
A note on inference of trait associations with SNP
haplotypes and other attributes in generalized linear models.
Human Heredity, \bold{57}:200-206

Burkett K, Graham J and McNeney B (2006). hapassoc: Software for 
Likelihood Inference of Trait Associations with SNP Haplotypes and 
Other Attributes. Journal of Statistical Software, \bold{16(2)}:1-19
}
\seealso{
\code{\link{hapassoc}},\code{\link{summary.hapassoc}}.
}
\keyword{methods}
