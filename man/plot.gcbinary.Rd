\name{plot.gcbinary}
\alias{plot.gcbinary}
\title{
S3 Method for Plotting a 'gcbinary' Object
}
\description{
Provides graphical summaries for an object returned by the function \code{gc_binary}. This S3 method offers two types of plots for calibration or proportions' representation.
}
\usage{
\method{plot}{gcbinary}(x, method="calibration", n.groups=5, smooth=FALSE, ...)
}
\arguments{
\item{x}{An object returned by the function \code{gc_binary}.}
\item{method}{A character string specifying the type of plot. Possible values are "calibration" (default) or "proportion".}
\item{n.groups}{An integer for the number of groups when \code{method="calibration"} is used. The default is 5.}
\item{smooth}{A logical value by default equals to \code{FALSE}. If \code{boot.mi=TRUE} in the object \code{x}, enables plotting a single smooth calibration curve averaged over all \code{m} values instead of separate curves.}
\item{...}{Additional graphical parameters that can be passed to the underlying plot function.}
}
\details{
If \code{n.groups} is too high, some groups may have too few events to reliably estimate proportions. An error can also occur if the Q-model is singular.
}
\value{No return value for this S3 method.}
\examples{
data(dataPROPHYVAP)

.f <- formula(VAP ~ GROUP * AGE + SEX + BMI + DIABETES)

### In practice use larger values for cv and boot.number
### We set boot.number and cv at 10 for speed in CRAN checks
### cv=30 and boot.number=1000 are more appropriate values
gc1 <- gc_binary(formula=.f, model="ridge", data=dataPROPHYVAP,
              group="GROUP", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", progress=TRUE , boot.tune=TRUE)

#Plot the calibration curve
plot(gc1, n.groups=3, col="red3")
}

\keyword{Plotting}
\keyword{Binary outcome}