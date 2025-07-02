\name{print.gcsurv}
\alias{print.gcsurv}

\title{
S3 Method for Printing an 'gcsurv' Object
}
\description{
Print a summary of the \code{gcsurv} object returned from gc_survival function
}

\usage{
\method{print}{gcsurv}(x,  digits=4, ...)
}

\arguments{
  \item{x}{An object returned by the function \code{gc_survival}.}
  \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
  \item{...}{For future methods.}
}

\value{No return value for this S3 method.}

\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$DEATH_J60_num <- ifelse(dataPROPHYVAP$DEATH_J60 == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

formula <- formula(Surv(FOLLOWUP_J60, DEATH_J60_num) ~ GROUP_num * AGE + SEXE +
              BMI + DIABETES)

set.seed(12251)
gc1 <- gc_survival(formula=formula, method="lasso", data=dataPROPHYVAP,
group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10, boot.number=500,
              effect="ATE", progress=TRUE , max.time=10, boot.tune = FALSE)

print(gc1)
}

\keyword{Printing}
