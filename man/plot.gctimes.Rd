\name{plot.gctimes}
\alias{plot.gctimes}
\title{S3 Method for Plotting a 'gctimes' Object}

\description{
Provides calibration and survival plots for an object returned by the function \code{gc_times}. This method assesses how well the Q-model's predicted survival aligns with observed time-to-event outcomes.
}

\usage{
\method{plot}{gctimes}(x, method="calibration", n.groups=5, pro.time=NULL, smooth=FALSE, ...)
}

\arguments{
  \item{x}{An object returned by the function \code{gc_times}.}
  \item{method}{A character string specifying the type of plot. Options are:
                 \code{"calibration"} (default) for calibration at \code{pro.time}, 
                 \code{"calibration2"} for Kaplan-Meier survival vs. mean predicted survival,
                 \code{"survival"} for predicted survival curves by treatment/exposure group.}
  \item{n.groups}{An integer for the number of quantiles to divide the predicted probabilities into. The default is 5. Note: If \code{n.groups} is set too high, some groups may have too few events to reliably estimate survivals. An error can also occur if the Q-model is singular.}
  \item{pro.time}{A numeric value specifying the time at which calibration or RMST is evaluated. Defaults to the \code{pro.time} used in \code{gc_times}.}
  \item{smooth}{A logical value, by default \code{FALSE}. If \code{TRUE}, applies smoothing to the calibration curves.}
  \item{...}{Additional graphical parameters that can be passed to the underlying plot functions.}
}

\details{
\strong{Methods:}  
\itemize{
  \item \code{"calibration"}: Visualizes calibration at \code{pro.time} by dividing predicted probabilities into \code{n.groups} quantiles. For each group, the average predicted probability is plotted against the observed proportion of events at \code{pro.time}, including 95\% confidence intervals. An identity line is provided for reference, and perfect calibration is indicated when the observed points lie directly along this line.
  
  \item \code{"calibration2"}: Plots the Kaplan-Meier estimations of the overall survival against the mean of the survival predictions derived from the Q-model. Note that this method cannot be used when \code{boot.mi=TRUE}.
  
  \item \code{"survival"}: Plots the predicted survival curves derived from the Q-model for each treatment/exposure group.
}
 
Survival predictions are computed using the linear predictor from the fitted model and the baseline cumulative hazard estimated via Breslow's method. Bootstrap iterations are used in \code{gc_times} to compute adjusted and unadjusted survival curves, RMST, and differences between groups.
}

\value{No return value for this S3 method.}

\examples{
data("dataPROPHYVAP")

.f <- formula(Surv(TIME_DEATH, DEATH) ~ GROUP * (AGE + SEX + ALCOHOL + BMI + DIABETES))

### In practice use larger values for cv and boot.number
### We set boot.number and cv at 10 for speed in CRAN checks
### cv=30 and boot.number=1000 are more appropriate values
gc1 <- gc_times(formula=.f, model="lasso", data=dataPROPHYVAP, group="GROUP",
                cv=10, boot.type="boot", boot.number=10, boot.tune=TRUE,
                effect="ATE", progress=TRUE, pro.time=20, seed=5192)

#Calibration plot at 20 days
plot(gc1, method="calibration", n.groups=5, pro.time=10)

#Kaplan-Meier vs mean predicted survival
plot(gc1, method="calibration2")

#Predicted survival curves for treatment groups
plot(gc1, method="survival", col=c("red3","blue3"))
legend("bottomleft", c("Placebo", "Ceftriaxone"), col=c("red3","blue3"), lty=1)
}

\keyword{Plotting}
\keyword{Time-to-event outcome}
\keyword{Survival analysis}
\keyword{Calibration}