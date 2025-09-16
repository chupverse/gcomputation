\name{print.gclogi}
\alias{print.gclogi}
\title{

S3 Method for Printing an 'gclogi' Object

}
\description{

Print a summary of the \code{gclogi} object returned from gc_logistic function.

}
\usage{

\method{print}{gclogi}(x, digits=4, ...)

}
\arguments{

\item{x}{An object returned by the function \code{gc_logistic}.}

\item{digits}{An optional integer for the number of digits to print when printing numeric values.}

\item{...}{For future methods.}

}
\value{No return value for this S3 method.}
\examples{

data(dataPROPHYVAP)

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

.f <- formula(VAP_num ~ GROUP_num * AGE + SEX + BMI + DIABETES)

gc1 <- gc_logistic(formula=.f, method="ridge", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=500, effect="ATE", progress=TRUE, boot.tune = TRUE)
print(gc1)

}
\keyword{Printing}