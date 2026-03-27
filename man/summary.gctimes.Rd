\name{summary.gctimes}
\alias{summary.gctimes}
\alias{print.summary.gctimes}
\title{
S3 Method for Summarizing an 'gctimes' Object
}
\description{
Summarize an object returned by the function \code{gc_times}.
}
\usage{
\method{summary}{gctimes}(object, digits=4, ci.type=NULL, ci.level=0.95,
                         unadjusted=TRUE, ...)
}
\arguments{
\item{object}{An object returned by the function \code{gc_times}.}
\item{digits}{An optional integer for the number of digits to print when summarizing numeric values.}
\item{ci.type}{The type of confidence intervals. Two types are possible: "norm" (assumed the Normal distribution) and "perc" (non-parametric estimation by percentiles).}
\item{ci.level}{The confidence level required. Default is 95\%.}
\item{unadjusted}{A logical value to print the unadjusted results. The default is \code{TRUE}}
\item{...}{For future methods.}
}

\value{
  \item{adjusted}{The data frame of the G-computation summary results.}
  \item{unadjusted}{The data frame of the unadjusted summary results.}
  \item{model}{The method used for the Q-model.}
  \item{formula}{The formula provided by the user.}
  \item{tuning.parameters}{The estimated tuning parameters for the Q-model. For "aic" and "bic" methods, this represents the final model.}
  \item{n}{The sample size of the dataset after missing data removal (if \code{boot.mi = FALSE}).}
  \item{nevent}{The total number of events in the dataset.}
  \item{nimput}{The number of imputations to perform if boot.mi is \code{TRUE}.}
  \item{missing}{The number of observations that were removed from the original dataset due to missing values.}
  \item{m}{The number of multiple imputations performed (only present if \code{mi.boot = TRUE}).}
  \item{pro.time}{The time point of interest for the evaluation.}
  \item{digits}{The number of digits to print.}
  \item{unadjusted.flag}{A logical value to indicate if the unadjusted results are provided.}
}

\examples{
data(dataPROPHYVAP)

formula <- formula(Surv(TIME_DEATH, DEATH) ~ GROUP * AGE +
              SEX + BMI + DIABETES)

### In practice use larger values for cv and boot.number
### We set boot.number and cv at 10 for speed in CRAN checks
### cv=30 and boot.number=1000 are more appropriate values
gc1 <- gc_times(formula=formula, model="lasso", data=dataPROPHYVAP,
              group="GROUP", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", pro.time=20, boot.tune=TRUE)
              
summary(gc1, ci.type="norm", ci.level=0.99)
}

\keyword{Summarizing}
\keyword{Time-to-event outcome}