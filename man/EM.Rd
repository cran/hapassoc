\name{EM}
\alias{EM}
\title{Use the EM algorithm to fit maximum likelihood estimates of trait associations with SNP haplotypes}
\description{
This function takes a data set of haplotypes in which rows for
individuals of uncertain phase have been augmented by "pseudo-individuals"
carrying consistent multilocus genotypes.
The function uses the EM algorithm to find
MLE's for trait associations with covariates in generalized linear models.
}
\usage{EM(form,haplos.list,baseline = "missing" ,family = binomial(),
gamma = FALSE, maxit = 50, tol = 0.001, \ldots)}
\arguments{
\item{form}{model equation in usual R format}
\item{haplos.list}{list of haplotype data from \code{\link{CheckHaplos}}}
\item{baseline}{optional, haplotype to be used for baseline coding. Default
is the most frequent haplotype.}
\item{family}{binomial, poisson, gaussian or gamma are supported, default=binomial}
\item{gamma}{initial values for haplotype frequencies, default values are calculated in \code{\link{CheckHaplos}} using the standard haplotype-counting EM algorithm without adjustment for non-haplotype covariates}
\item{maxit}{maximum iterations of the EM loop, default=50}
\item{tol}{convergence tolerance in terms of the maximum difference in parameter estimates between interations; default=0.001}
\item{\ldots}{additional arguments to be passed to the glm function such as starting values}
}
\value{
\item{it}{number of iterations}
\item{beta}{estimated regression coefficients}
\item{gamma}{estimated haplotype frequencies}
\item{fits}{fitted values of the trait}
\item{wts}{final weights calculated in last iteration of the EM loop. These are estimates of the
conditional probabilities of each multilocus genotype given the observed 
single-locus genotypes.}
\item{var}{joint variance-covariance matrix of the estimated regression
coefficients and the estimated haplotype frequencies}
\item{dispersionML}{maximum likelihood estimate of dispersion parameter}
\item{family}{family of the generalized linear model fit (e.g. binomial, gaussian, etc.)}
\item{response}{trait value}
}
\examples{
data(hypoDat)
example.haplos<-RecodeHaplos(hypoDat,2)
example.newhaplos<-CheckHaplos(example.haplos, 2)

# Logistic regression, baseline group: '10/10'

example.regr <- EM(affected ~ gender + age + BMI + n00+ n01 + n11,
                     example.newhaplos, family=binomial())
}
\seealso{
\code{\link{RecodeHaplos}},\code{\link{CheckHaplos}},\code{\link{summary.EM}},\code{\link{glm}},\code{\link{family}}.
}
\keyword{methods}
