\name{summary.gcsurv}
\alias{summary.gcsurv}
\title{

S3 Method for Summarizing an 'gcsurv' Object

}
\description{

Summarize an object returned by the function \code{gc_survival}.

}
\usage{

\method{summary}{gcsurv}(object, digits=4, ci.type=NULL, ci.level=0.95, ...)

}
\arguments{

\item{object}{An object returned by the function \code{gc_survival}.}

\item{digits}{An optional integer for the number of digits to print when summarizing numeric values.}

\item{ci.type}{The type of confidence intervals to calculate. Two types are possible: "norm" (asymptotic) and "perc" (non parametric).}

\item{ci.level}{The confidence level required.}

\item{...}{For future methods.}

}
\value{
\item{GC}{The data.frame of the G-computation summary results.}
\item{unadjusted}{The data.frame of the unadjusted summary results.}
}
\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$DEATH_num <- ifelse(dataPROPHYVAP$DEATH == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

formula <- formula(Surv(TIME_DEATH, DEATH_num) ~ GROUP_num * AGE +
              SEX + BMI + DIABETES)
              
gc1 <- gc_survival(formula=formula, method="lasso", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=500, effect="ATE", progress=TRUE , pro.time=10,
              boot.tune = TRUE)
summary(gc1)

}
\keyword{Summarizing}