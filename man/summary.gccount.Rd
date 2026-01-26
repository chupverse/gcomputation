\name{summary.gccount}
\alias{summary.gccount}
\title{
S3 Method for Summarizing an 'gccount' Object
}
\description{
Summarize an object returned by the function \code{gc_count}.
}
\usage{
\method{summary}{gccount}(object, digits=4, ci.type=NULL, ci.level=0.95,
                         unadjusted=TRUE, ...)
}
\arguments{
\item{object}{An object returned by the function \code{gc_count}.}
\item{digits}{An optional integer for the number of digits to print when summarizing numeric values.}
\item{ci.type}{The type of confidence intervals. Two types are possible: "norm" (assumed the Normal distribution) and "perc" (non-parametric estimation by percentiles).}
\item{ci.level}{The confidence level required. Default is 95\%.}
\item{unadjusted}{A logicical value to print the unadjusted results. The default is \code{TRUE}}
\item{...}{For future methods.}
}
\value{
\item{GC}{The data frame of the G-computation summary results.}
\item{unadjusted}{The data frame of the unadjusted summary results.}
}
\examples{
data("dataPROPHYVAP")

.f <- formula(VAP ~ GROUP * (AGE + SEX + ALCOHOL + BMI + DIABETES))

### In practice use larger values of boot.number (e.g., 500)
### We set boot.number and cv at 10 for speed in CRAN checks
gc1 <- gc_count(formula=.f, model="all", data=dataPROPHYVAP, group="GROUP",
                    cv=10, boot.type="boot", boot.number=10, boot.tune=TRUE,
                    effect="ATE")

summary(gc1,vci.type="norm", ci.level=0.99)
}

\keyword{Summarizing}
\keyword{Continuous outcome}