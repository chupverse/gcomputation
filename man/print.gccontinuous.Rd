\name{print.gccontinuous}
\alias{print.gccontinuous}
\title{

S3 Method for Printing an 'gccontinuous' Object

}
\description{

Print a summary of the \code{gccontinuous} object returned from gc_binary function.

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

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

beta_ttt <- 10
beta_alcool <- 2
beta_bmi <- 8
beta_age <- 5

cont_outcome <- rnorm(n = nrow(dataPROPHYVAP), 
                      mean = beta_ttt*dataPROPHYVAP$GROUP_num +
                        beta_alcool*match(as.factor(dataPROPHYVAP$ALCOOL),
                            levels(as.factor(dataPROPHYVAP$ALCOOL))) + 
                        beta_bmi*dataPROPHYVAP$BMI + 
                        beta_age*dataPROPHYVAP$AGE, 
                      sd = 1)

dataPROPHYVAP$y_cont <- cont_outcome

.f <- formula(y_cont ~ GROUP_num * (AGE + SEX + ALCOOL + BMI + DIABETES))

### In practice use larger values of boot.number (e.g., 500)
### We set boot.number = 10 for speed in CRAN checks
gc1 <- gc_continuous(formula=.f, model="all", data=dataPROPHYVAP, group="GROUP_num",
                    cv=10, boot.type="boot", boot.number=10, boot.tune=TRUE,
                    effect="ATE", progress=TRUE, seed=5192)

print(gc1)
}
\keyword{Printing}