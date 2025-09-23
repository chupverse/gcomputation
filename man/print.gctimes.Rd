\name{print.gctimes}
\alias{print.gctimes}

\title{
S3 Method for Printing an 'gctimes' Object
}
\description{
Print a summary of the \code{gctimes} object returned from gc_times function
}

\usage{
\method{print}{gctimes}(x,  digits=4, ...)
}

\arguments{
  \item{x}{An object returned by the function \code{gc_times}.}
  \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
  \item{...}{For future methods.}
}

\value{No return value for this S3 method.}

\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$DEATH_num <- ifelse(dataPROPHYVAP$DEATH == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

formula <- formula(Surv(TIME_DEATH, DEATH_num) ~ GROUP_num * AGE + SEX +
              BMI + DIABETES)

set.seed(12251)
gc1 <- gc_times(formula=formula, model="lasso", data=dataPROPHYVAP,
group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10, boot.number=500,
              effect="ATE", progress=TRUE , pro.time=10, boot.tune=FALSE)

print(gc1)
}

\keyword{Printing}
