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

.f <- formula(Surv(TIME_DEATH, DEATH) ~ GROUP * AGE + SEX +
              BMI + DIABETES)

### In practice use larger values for cv and boot.number
### We set boot.number and cv at 10 for speed in CRAN checks
gc1 <- gc_times(formula=.f, model="lasso", data=dataPROPHYVAP,
group="GROUP", param.tune=NULL, boot.type="bcv", cv=10, boot.number=10,
              effect="ATE", progress=TRUE , pro.time=10, boot.tune=FALSE)

print(gc1)
}

\keyword{Printing}
\keyword{Time-to-event outcome}