\name{print.gcbinary}
\alias{print.gcbinary}
\title{

S3 Method for Printing an 'gcbinary' Object

}
\description{

Print a summary of the \code{gcbinary} object returned from gc_binary function.

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

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

.f <- formula(VAP_num ~ GROUP_num * AGE + SEX + BMI + DIABETES)

gc1 <- gc_binary(formula=.f, model="ridge", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", progress=TRUE, boot.tune=TRUE)
print(gc1)

}
\keyword{Printing}