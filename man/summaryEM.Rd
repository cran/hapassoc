\name{summary.EM}
\alias{summary.EM}
\title{Summarize results of the EM function}
\description{
Summary function for reporting the results of the EM function in a similar  
style to the lm and glm summaries.
}
\usage{summary.EM(object, \ldots)}
\arguments{
\item{object}{a list of class EM}
\item{\ldots}{additional arguments to the summary function (currently unused)}
}
\value{
\item{coefficients}{Table of estimated coefficients, standard errors and Wald tests for each variable}
\item{frequencies}{Table of estimated haplotype frequencies and standard errors}
\item{dispersion}{Estimate of dispersion parameter (Moment estimator for gamma model)}
}
\examples{
data(hypoDat)
example.haplos<-RecodeHaplos(hypoDat,2)
example.newhaplos<-CheckHaplos(example.haplos, 2)
example.regr <- EM(affected ~ gender + age + BMI + n00+ n01 + n11,
                     example.newhaplos, family=binomial())
# Summarize the results of the function EM
summary.EM(example.regr) # or just summary(example.regr)

# Results:
#$coefficients
#               Estimate Std. Error     zscore     Pr(>|z|)
#(Intercept) -7.49563443 1.48540169 -5.0462003 4.506830e-07
#gender      -0.21123481 0.20849812 -1.0131257 3.110001e-01
#age          0.03774116 0.01803805  2.0923083 3.641095e-02
#BMI          0.16280168 0.03890977  4.1840828 2.863198e-05
#n00         -0.24241680 1.24418218 -0.1948403 8.455180e-01
#n01         -0.26860899 0.19322030 -1.3901696 1.644774e-01
#n11          0.42651099 0.17602854  2.4229650 1.539441e-02
#
#$frequencies
#         Estimate  Std. Error
#f.n00 0.009349114 0.004653122
#f.n01 0.246650886 0.013743345
#f.n10 0.437650886 0.015693892
#f.n11 0.306349114 0.014638656
#
#$dispersion
#[1] 1
}
\seealso{
\code{\link{RecodeHaplos}},\code{\link{CheckHaplos}},\code{\link{EM}}.
}
\keyword{methods}
