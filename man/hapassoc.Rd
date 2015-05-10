\name{hapassoc}
\alias{hapassoc}
\title{EM algorithm to fit maximum likelihood estimates of trait associations with SNP haplotypes}
\description{
This function takes a dataset of haplotypes in which rows for
individuals of uncertain phase have been augmented by \dQuote{pseudo-individuals}
who carry the possible multilocus genotypes consistent with
the single-locus phenotypes. For cohort or cross-sectional
data, the EM algorithm is used to find
MLE's for trait associations with covariates in generalized linear models.
For case-control data, the algorithm solves a set of unbiased
estimating equations (see \bold{Details}).
}
\usage{hapassoc(form,haplos.list,baseline = "missing" ,family = binomial(),
design="cohort",disease.prob=NULL,freq = NULL, maxit = 50, tol = 0.001, 
start = NULL, verbose=FALSE)}
\arguments{
\item{form}{model equation in usual R format}
\item{haplos.list}{list of haplotype data from \code{\link{pre.hapassoc}}}
\item{baseline}{optional, haplotype to be used for baseline coding if the model formula
either includes all haplotypes or is of the form "y~." for example. Default
is the most frequent haplotype according to the initial haplotype frequency
estimates returned by \code{\link{pre.hapassoc}}.}
\item{family}{binomial, poisson, gaussian or gamma are supported, 
default=binomial}
\item{design}{study design. Default is \dQuote{cohort} for cohort or 
cross-sectional sampling. Users may optionally specify \dQuote{cc} for 
case-control or retrospective sampling of exposures (i.e. genotypes and 
non-genetic attributes) conditional on disease status. 
When \code{design="cc"}, \code{family=binomial()} is assumed
and the robust MPSE estimator of the regression parameters 
(Spinka et al., 2005) is returned; see \bold{Details} for more information.
}

\item{disease.prob}{marginal disease probability [P(D=1)]  to use in 
the MPSE estimator, if \code{design="cc"}. If \code{disease.prob=NULL} 
(the default value), a rare disease is assumed.
This argument is ignored if \code{design="cohort"}.}
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
\item{dispersion}{maximum likelihood estimate of dispersion parameter 
(to get the moment estimate, use \code{\link{summary.hapassoc}})
if applicable, otherwise 1}
\item{family}{family of the generalized linear model (e.g. binomial, 
gaussian, etc.)}
\item{response}{trait value}
\item{converged}{TRUE/FALSE indicator of convergence. If the algorithm 
fails to converge, only the \code{converged} indicator is returned.}
\item{model}{model equation}
\item{loglik}{the log-likelihood evaluated at the maximum likelihood estimates 
of all parameters if \code{design="cohort"}, or \code{NA} if 
\code{design="cc"}}
\item{call}{the function call}
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

Burkett K, Graham J and McNeney B (2006). hapassoc: Software for
Likelihood Inference of Trait Associations with SNP Haplotypes and Other
Attributes. Journal of Statistical Software, \bold{16(2)}:1-19

Chen, Z. (2006): Approximate likelihood inference for haplotype risks
in case-control studies of a rare disease, Masters thesis, Statistics
and Actuarial Science, Simon Fraser University, available at
\url{https://www.stat.sfu.ca/content/dam/sfu/stat/alumnitheses/MiscellaniousTheses/Chen-2006.pdf}.

Spinka, C., Carroll, R. J. & Chatterjee, N. (2005). Analysis of 
case-control studies of genetic and environmental factors with missing 
genetic information and haplotype-phase ambiguity. 
Genetic Epidemiology, \bold{29}, 108-127. 
}
\seealso{
\code{\link{pre.hapassoc}},\code{\link{summary.hapassoc}},\code{\link{glm}},\code{\link{family}}.
}
\details{
See the hapassoc vignette, of the same name as the package, for details.

When the study design is case-control, i.e. genotypes and
non-genetic attributes have been sampled retrospectively given disease
status, naive application of prospective maximum likelihood methods
can yield biased inference (Spinka et al., 2005, Chen, 2006). 
Therefore, when \code{design="cc"}, 
the algorithm solves the modified prospective score equations or MPSE
(Spinka et al. 2005) for regression
and haplotype frequency parameters. The implementation in
\pkg{hapassoc} is due to Chen (2006). 
In general, the MPSE approach requires that
the marginal probability of disease,  P(D=1), be known. 
An exception is when the disease is rare; hence, when 
\code{disease.prob=NULL} (the default) a rare disease is assumed.
The variance-covariance matrix of the regression parameter and
haplotype frequency estimators is approximated as described
in Chen (2006). Limited simulations indicate that the resulting standard errors
for regression parameters perform well, but not the
standard errors for haplotype frequencies, which should be ignored.
For case-control data, we hope to
implement the variance-covariance estimator of Spinka et al. (2005)
in a future version of \pkg{hapassoc}.
}
\keyword{methods}
