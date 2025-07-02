\name{summary.gclogi}
\alias{summary.gclogi}
\title{

S3 Method for Summarizing an 'gclogi' Object

}
\description{

Summarize an object returned by the function \code{gc_logistic}.

}
\usage{

\method{summary}{gclogi}(object, digits=4, ...)

}
\arguments{

\item{object}{An object returned by the function \code{gc_logistic}.}

\item{digits}{An optional integer for the number of digits to print when summarizing numeric values.}

\item{...}{For future methods.}

}
\value{No return value for this S3 method.}
\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$VAP_EARLY_num <- ifelse(dataPROPHYVAP$VAP_EARLY == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

formula <- formula(VAP_EARLY_num ~ GROUP_num * AGE + SEXE + BMI + DIABETES)

gc1 <- gc_logistic(formula=formula, method="ridge", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=500, effect="ATE", progress=TRUE, boot.tune = TRUE)
summary(gc1)

}
\keyword{Summarizing}