\name{plot.gctimes}
\alias{plot.gctimes}
\title{
S3 Method for Plotting a 'gctimes' Object
}
\description{
Provides graphical summaries for an object returned by the function \code{gc_times}. This S3 method offers two types of calibration plots for calibration or survivals' representation.
}
\usage{
\method{plot}{gctimes}(x, method="calibration", n.groups=5,
                      pro.time=NULL, smooth=FALSE, ...)
}
\arguments{
\item{x}{An object returned by the function \code{gc_times}.}
\item{method}{A character string specifying the type of plot. Possible values are "calibration" (default) or "survival". See details.}
\item{n.groups}{An integer for the number of groups when \code{method="calibration"} is used. The default is 5.}
\item{pro.time}{If \code{method="calibration"}, the prognostic time at which the calibration plot is performed. If \code{NULL} (default), the median is used.}
\item{smooth}{A logical value by default equals to \code{FALSE}. If \code{boot.mi=TRUE} in the object \code{x}, it enables plotting a single smooth calibration curve averaged over all \code{m} values instead of separate curves.}
\item{...}{Additional graphical parameters that can be passed to the underlying plot function.}
}
\details{
This method supports two distinct plots:
\itemize{
\item \code{method="survival"}: This plot compares the Kaplan-Meier estimations of the overall survival against the mean of the survival predictions derived from the Q-model.
\item \code{method="calibration"}: This plot compares the mean of the predicted survivals at \code{pro.time} against the Kaplan-Meier estimations for \code{n.groups} defined from the related percentiles. If \code{n.groups} is set too high, some groups may have too few events to reliably estimate survivals.
}
}
\value{No return value for this S3 method.}
\examples{
data(dataPROPHYVAP)

.f <- formula(Surv(TIME_DEATH, DEATH) ~ GROUP + AGE + SEX +
              BMI + DIABETES)
              
### In practice use larger values for cv and boot.number
### We set boot.number and cv at 10 for speed in CRAN checks
gc1 <- gc_times(formula=.f, model="lasso", data=dataPROPHYVAP,
              group="GROUP", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", progress=TRUE , pro.time=10,
              boot.tune=TRUE)

#Plot the calibration curve
plot(gc1, method="calibration", n.groups=5)

#Plot the survival curves
plot(gc1, method="survival")
}

\keyword{Plotting}
\keyword{Time-to-event outcome}