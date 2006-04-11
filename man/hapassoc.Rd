\name{hapassoc}
\alias{hapassoc}
\title{EM algorithm to fit maximum likelihood estimates of trait associations with SNP haplotypes}
\description{
This function takes a dataset of haplotypes in which rows for
individuals of uncertain phase have been augmented by \dQuote{pseudo-individuals}
who carry the possible multilocus genotypes consistent with
the single-locus phenotypes.  The EM algorithm is used to find
MLE's for trait associations with covariates in generalized linear models.
}
\usage{hapassoc(form,haplos.list,baseline = "missing" ,family = binomial(),
freq = NULL, maxit = 50, tol = 0.001, start = NULL, verbose=FALSE)}
\arguments{
\item{form}{model equation in usual R format}
\item{haplos.list}{list of haplotype data from \code{\link{pre.hapassoc}}}
\item{baseline}{optional, haplotype to be used for baseline coding. Default
is the most frequent haplotype according to the initial haplotype frequency
estimates returned by \code{\link{pre.hapassoc}}.}
\item{family}{binomial, poisson, gaussian or freq are supported, 
default=binomial}
\item{freq}{initial estimates of haplotype frequencies, default values are 
calculated in \code{\link{pre.hapassoc}} using standard haplotype-counting 
(i.e. EM algorithm without adjustment for non-haplotype covariates)}
\item{maxit}{maximum number of iterations of the EM algorithm; default=50}
\item{tol}{convergence tolerance in terms of either the maximum difference in 
parameter estimates between interations or the maximum relative
difference in parameter estimates between iterations, which ever is larger.}
\item{start}{starting values for parameter estimates in the risk model}
\item{verbose}{should the iteration number and value of the covergence criterion be printed
at each iteration of the EM algorithm? Default=FALSE}
}
\value{
\item{it}{number of iterations of the EM algorithm}
\item{beta}{estimated regression coefficients}
\item{freq}{estimated haplotype frequencies}
\item{fits}{fitted values of the trait}
\item{wts}{final weights calculated in last iteration of the EM algorithm. These are estimates of the
conditional probabilities of each multilocus genotype given the observed 
single-locus genotypes.}
\item{var}{joint variance-covariance matrix of the estimated regression
coefficients and the estimated haplotype frequencies}
\item{dispersionML}{maximum likelihood estimate of dispersion parameter 
(to get the moment estimate, use \code{\link{summary.hapassoc}})}
\item{family}{family of the generalized linear model (e.g. binomial, 
gaussian, etc.)}
\item{response}{trait value}
\item{converged}{TRUE/FALSE indicator of convergence. If the algorithm 
fails to converge, only the \code{converged} indicator is returned.}
\item{model}{model equation}
\item{loglik}{The log-likelihood evaluated at the maximum likelihood estimates 
of all parameters}
}
\examples{
data(hypoDat)
example.pre.hapassoc<-pre.hapassoc(hypoDat, 3)

example.pre.hapassoc$initFreq # look at initial haplotype frequencies
#      h000       h001       h010       h011       h100       h101       h110 
#0.25179111 0.26050418 0.23606001 0.09164470 0.10133627 0.02636844 0.01081260 
#      h111 
#0.02148268 


names(example.pre.hapassoc$haploDM)
# "h000"   "h001"   "h010"   "h011"   "h100"   "pooled"

# Columns of the matrix haploDM score the number of copies of each haplotype 
# for each pseudo-individual.

# Logistic regression for a multiplicative odds model having as the baseline 
# group homozygotes '001/001' for the most common haplotype

example.regr <- hapassoc(affected ~ attr + h000+ h010 + h011 + h100 + pooled,
                  example.pre.hapassoc, family=binomial())

# Logistic regression with separate effects for 000 homozygotes, 001 homozygotes 
# and 000/001 heterozygotes

example2.regr <- hapassoc(affected ~ attr + I(h000==2) + I(h001==2) +
                   I(h000==1 & h001==1), example.pre.hapassoc, family=binomial())

}
\references{Burkett K, McNeney B, Graham J (2004).
A note on inference of trait associations with SNP
haplotypes and other attributes in generalized linear models.
Human Heredity, \bold{57}:200-206

Burkett K, Graham J and McNeney B. hapassoc: Software for likelihood
inference of trait associations with SNP haplotypes and other
attributes. Under revision for Journal of Statistical Software.
}
\seealso{
\code{\link{pre.hapassoc}},\code{\link{summary.hapassoc}},\code{\link{glm}},\code{\link{family}}.
}
\keyword{methods}
