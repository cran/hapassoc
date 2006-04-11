\name{anova.hapassoc}
\alias{anova.hapassoc}
\title{Return likelihood ratio test of haplotype effect}
\description{
This function returns the likelihood ratio test statistic comparing
two nested models fit with \code{hapassoc}
}
\usage{\method{anova}{hapassoc}(object, redfit, display=TRUE, \ldots)}
\arguments{
\item{object}{a list of class \code{hapassoc} output by the 
\code{\link{hapassoc}} function}.
\item{redfit}{A \code{hapassoc} object resulting from fitting a reduced 
model}
\item{display}{An indicator to suppress output displayed on screen} 
\item{\ldots}{additional arguments to the summary function currently unused}
}
\value{
\item{LRTstat}{The likelihood ratio statistic comparing the two models}
\item{df}{Degrees of freedom of the likelihood ratio statistic}
\item{pvalue}{The p-value of the test}
}
\examples{
data(hypoDatGeno)
example2.pre.hapassoc<-pre.hapassoc(hypoDatGeno, numSNPs=3, allelic=FALSE)
example2.regr <- hapassoc(affected ~ attr + hAAA+ hACA + hACC + hCAA + 
pooled, example2.pre.hapassoc, family=binomial())
example2.regr2 <- hapassoc(affected ~ attr + hAAA, example2.pre.hapassoc, 
family=binomial())
anova(example2.regr,example2.regr2)

# Returns:

#	hapassoc: likelihood ratio test

#Full model: affected ~ attr + hAAA + hACA + hACC + hCAA + pooled 
#Reduced model: affected ~ attr + hAAA 

#LR statistic = 1.5433 , df = 4 , p-value =  0.8189 
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
\code{\link{pre.hapassoc}},\code{\link{hapassoc}}, 
\code{\link{summary.hapassoc}}.
}
\keyword{methods}
