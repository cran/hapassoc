\name{logLik.hapassoc}
\alias{logLik.hapassoc}
\title{Return log-likelihood}
\description{
This function is used to return the log-likelihood at the maximum 
likelihood estimates computed by \code{hapassoc} and to return
the number of parameters fit by \code{hapassoc} 
(i.e. the degrees of freedom in \code{R}) for cohort or
cross-sectional data. 
}
\usage{\method{logLik}{hapassoc}(object, \ldots)}
\arguments{
\item{object}{a list of class \code{hapassoc} output by the 
\code{\link{hapassoc}} function}
\item{\ldots}{additional arguments to the summary 
function (currently unused)}
}
\value{
\item{logLik}{log-likelihood computed at the maximum likelihood estimates
if \code{design="cohort"}, or \code{NA} if \code{design="cc"}}
\item{df}{number of parameters in the model (i.e. regression coefficients,
any dispersion parameters and haplotype 
frequencies). This is not the residual
degrees of freedom, which is the number of subjects minus the number
of parameters estimated.}
}
\examples{
data(hypoDatGeno)
example2.pre.hapassoc<-pre.hapassoc(hypoDatGeno, numSNPs=3, allelic=FALSE)
example.regr <- hapassoc(affected ~ attr + hAAA+ hACA + hACC + hCAA + 
pooled, example2.pre.hapassoc, family=binomial())
logLik(example.regr)

# Returns:
# Log Lik: -322.1558 (df=14)
}
\references{Burkett K, McNeney B, Graham J (2004).
A note on inference of trait associations with SNP
haplotypes and other attributes in generalized linear models.
Human Heredity, \bold{57}:200-206

Burkett K, Graham J and McNeney B (2006). hapassoc: Software for
Likelihood Inference of Trait Associations with SNP Haplotypes and Other
Attributes. Journal of Statistical Software, \bold{16(2)}:1-19
}
\details{
See the hapassoc vignette, of the same name as the package, for details.
}
\seealso{
\code{\link{pre.hapassoc}},\code{\link{hapassoc}},
\code{\link{summary.hapassoc}}.
}
\keyword{methods}
