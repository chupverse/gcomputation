\name{print.gcbinary}
\alias{print.gcbinary}
\title{
S3 Method for Printing an 'gcbinary' Object
}
\description{
Print a summary of an object returned from \code{gc_binary} function.
}
\usage{
\method{print}{gcbinary}(x, digits=4, ...)
}
\arguments{
\item{x}{An object returned by the function \code{gc_binary}.}
\item{digits}{An optional integer for the number of digits to print when printing numeric values.}
\item{...}{For future methods.}
}
\value{No return value for this S3 method.}
\examples{
data(dataPROPHYVAP)

.f <- formula(VAP ~ GROUP * AGE + SEX + BMI + DIABETES)

### In practice use larger values for cv and boot.number
### We set boot.number and cv at 10 for speed in CRAN checks
gc1 <- gc_binary(formula=.f, model="ridge", data=dataPROPHYVAP,
              group="GROUP", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", progress=TRUE, boot.tune=TRUE)

print(gc1)
}

\keyword{Printing}
\keyword{Binary outcome}