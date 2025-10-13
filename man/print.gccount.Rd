\name{print.gccount}
\alias{print.gccount}
\title{
S3 Method for Printing an 'gccontinuous' Object
}
\description{
Print a summary of and object returned from \code{gc_count} function.
}
\usage{
\method{print}{gccount}(x, digits=4, ...)
}
\arguments{
\item{x}{An object returned by the function \code{gc_count}.}
\item{digits}{An optional integer for the number of digits to print when printing numeric values.}
\item{...}{For future methods.}
}
\value{No return value for this S3 method.}
\examples{
data("dataPROPHYVAP")

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

.f <- formula(VAP_num ~ GROUP_num * (AGE + SEX + ALCOOL + BMI + DIABETES))

### In practice use larger values of boot.number (e.g., 500)
### We set boot.number and cv at 10 for speed in CRAN checks
gc1 <- gc_count(formula=.f, model="all", data=dataPROPHYVAP, group="GROUP_num",
                    cv=10, boot.type="boot", boot.number=10, boot.tune=TRUE,
                    effect="ATE", progress=TRUE, seed=5192)

print(gc1)
}

\keyword{Printing}