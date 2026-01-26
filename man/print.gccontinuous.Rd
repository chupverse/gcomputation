\name{print.gccontinuous}
\alias{print.gccontinuous}
\title{
S3 Method for Printing an 'gccontinuous' Object
}
\description{
Print a summary of and object returned from \code{gc_continuous} function.
}
\usage{
\method{print}{gccontinuous}(x, digits=4, ...)
}
\arguments{
\item{x}{An object returned by the function \code{gc_continuous}.}
\item{digits}{An optional integer for the number of digits to print when printing numeric values.}
\item{...}{For future methods.}
}
\value{No return value for this S3 method.}
\examples{
data("dataPROPHYVAP")

.f <- formula(VAP ~ GROUP * (AGE + SEX + ALCOHOL + BMI + DIABETES))

### In practice use larger values of boot.number (e.g., 500)
### We set boot.number and cv at 10 for speed in CRAN checks
gc1 <- gc_continuous(formula=.f, model="all", data=dataPROPHYVAP, group="GROUP",
                    cv=10, boot.type="boot", boot.number=10, boot.tune=TRUE,
                    effect="ATE", progress=TRUE, seed=5192)

print(gc1)
}

\keyword{Printing}