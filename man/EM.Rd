\name{EM}
\alias{EM}
\title{EM algorithm to fit maximum likelihood estimates of trait associations with SNP haplotypes}
\description{
This function takes a dataset of haplotypes in which rows for
individuals of uncertain phase have been augmented by "pseudo-individuals"
who carry the possible multilocus genotypes consistent with
the single-locus phenotypes.  The EM algorithm is used to find
MLE's for trait associations with covariates in generalized linear models.
}
\usage{EM(form,haplos.list,baseline = "missing" ,family = binomial(),
gamma = FALSE, maxit = 50, tol = 0.001, \ldots)}
\arguments{
\item{form}{model equation in usual R format}
\item{haplos.list}{list of haplotype data from \code{\link{PreEM}}}
\item{baseline}{optional, haplotype to be used for baseline coding. Default
is the most frequent haplotype.}
\item{family}{binomial, poisson, gaussian or gamma are supported, default=binomial}
\item{gamma}{initial estimates of haplotype frequencies, default values are calculated in \code{\link{PreEM}} using standard haplotype-counting 
(i.e. EM algorithm without adjustment for non-haplotype covariates)}
\item{maxit}{maximum iterations of the EM loop, default=50}
\item{tol}{convergence tolerance in terms of the maximum difference in 
parameter estimates between interations; default=0.001}
\item{\ldots}{additional arguments to be passed to the glm function such 
as starting values for parameter estimates in the risk model}
}
\value{
\item{it}{number of iterations of the EM algorithm}
\item{beta}{estimated regression coefficients}
\item{gamma}{estimated haplotype frequencies}
\item{fits}{fitted values of the trait}
\item{wts}{final weights calculated in last iteration of the EM loop. These are estimates of the
conditional probabilities of each multilocus genotype given the observed 
single-locus genotypes.}
\item{var}{joint variance-covariance matrix of the estimated regression
coefficients and the estimated haplotype frequencies}
\item{dispersionML}{maximum likelihood estimate of dispersion parameter 
(to get the moment estimate, use \code{\link{summary.EM}})}
\item{family}{family of the generalized linear model (e.g. binomial, 
gaussian, etc.)}
\item{response}{trait value}
\item{converged}{TRUE/FALSE indicator of convergence. If the algorithm 
fails to converge, only the converged indicator is returned.}
}
\examples{
data(hypoDat)
example.preEM<-PreEM(hypoDat, 3)

names(example.preEM$haploDM)
# "h000"   "h001"   "h010"   "h011"   "h100"   "pooled"

# Logistic regression, baseline group: '001/001'

example.regr <- EM(affected ~ attr + h000+ h010 + h011 + h100 + pooled,
                     example.preEM, family=binomial())
}
\references{Burkett K, McNeney B, Graham J (2004).
A note on inference of trait associations with SNP
haplotypes and other attributes in generalized linear models.
Human Heredity, In press}}
\seealso{
\code{\link{PreEM}},\code{\link{summary.EM}},\code{\link{glm}},\code{\link{family}}.
}
\keyword{methods}
