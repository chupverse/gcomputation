\name{summary.gcbinary}
\alias{summary.gcbinary}
\title{

S3 Method for Summarizing an 'gcbinary' Object

}
\description{

Summarize an object returned by the function \code{gc_binary}.

}
\usage{

\method{summary}{gcbinary}(object, digits=4, ci.type=NULL, ci.level=0.95,
                         unadjusted=TRUE, ...)

}
\arguments{

\item{object}{An object returned by the function \code{gc_binary}.}

\item{digits}{An optional integer for the number of digits to print when summarizing numeric values.}

\item{ci.type}{The type of confidence intervals to calculate. Two types are possible: "norm" (asymptotic) and "perc" (non parametric).}

\item{ci.level}{The confidence level required.}

\item{unadjusted}{A logicical value to print the unadjusted results. The default is \code{TRUE}}

\item{...}{For future methods.}

}
\value{
\item{GC}{The data.frame of the G-computation summary results.}
\item{unadjusted}{The data.frame of the unadjusted summary results.}
}
\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

.f <- formula(VAP_num ~ GROUP_num * AGE + SEX + BMI + DIABETES)

gc1 <- gc_binary(formula=.f, model="ridge", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", progress=TRUE, boot.tune=TRUE)
summary(gc1)

summary(gc1,ci.type="norm",ci.level=0.99)
}
\keyword{Summarizing}