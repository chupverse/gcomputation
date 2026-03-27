\name{plot.gccount}
\alias{plot.gccount}
\title{S3 Method for Plotting a 'gccount' Object }

\description{
Provides a calibration plot for an object returned by the function \code{gc_count}. This method assesses how well the Q-model's predicted values align with the observed counting outcomes.
}

\usage{
\method{plot}{gccount}(x, n.groups=5, smooth=FALSE, ...)
}

\arguments{
  \item{x}{An object returned by the function \code{gc_count}.}
  \item{n.groups}{An integer for the number of groups to divide the predicted means into. The default is 5. Note: If \code{n.groups} is set too high, some groups may have too few observations. An error can also occur if the Q-model is singular.}
  \item{smooth}{A logical value, by default \code{FALSE}. If \code{boot.mi} is TRUE, enables plotting a single smooth calibration curve averaged over all \code{m} values instead of separate curves.}
  \item{...}{Additional graphical parameters that can be passed to the underlying plot function.}
}

\details{
The function visualizes the calibration of the Q-model by dividing predicted probabilities into \code{n.groups} quantiles. For each group, the average predicted probability is plotted against the observed proportion of events, including 95\% confidence intervals for the observed values. An identity line is provided for reference and perfect calibration is indicated when the observed points are directly along this line.
}

\value{No return value for this S3 method.}

\examples{
data("dataPROPHYVAP")

.f <- formula(VAP ~ GROUP * (AGE + SEX + ALCOHOL + BMI + DIABETES))

### In practice use larger values for cv and boot.number
### We set boot.number and cv at 10 for speed in CRAN checks
### cv=30 and boot.number=1000 are more appropriate values
gc1 <- gc_count(formula=.f, model="all", data=dataPROPHYVAP, group="GROUP", 
                cv=10, boot.type="boot", boot.number=10, boot.tune=TRUE, 
                effect="ATE", progress=TRUE, seed=5192)
                     
#Plot the calibration curve
plot(gc1, n.groups=5, col="red3")
}

\keyword{Plotting}
\keyword{Counting outcome}