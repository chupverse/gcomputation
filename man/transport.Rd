\name{transport}
\alias{transport}
\title{
Transport the Marginal for Another Populations
}
\description{
Applies an already fitted G-computation model (an object of class to estimate marginal effects for that new population.
}
\usage{
transport(object, newdata, effect="ATE", boot.number=500, seed=NULL)
}
\arguments{
  \item{object}{An object returned by the function \code{\link{gc_binary}},  \code{\link{gc_continuous}}, or \code{\link{gc_times}}.}
  \item{newdata}{The new data frame with all the covariates in the \code{formula} of the \code{object}.}
  \item{effect}{The type of the marginal effect to be estimated on the \code{newdata}. Three types are possible: "ATE" (by default), "ATT", and "ATU". These correspond to Average Treatment Effect on the Entire population, on the Treated, and on the Untreated, respectively.}
  \item{boot.number}{The number of bootstrap resamples to perform for estimating standard errors and confidence intervals of the marginal effects on the \code{newdata}. The default value is 500.}
  \item{seed}{A random seed to ensure reproducibility during the bootstrap process. If \code{NULL} (the default), a seed is randomly assigned.}
}
\value{
The function returns an object of the same class with the same arguments as the input \code{object}, but with updated results for the \code{newdata}.
}
\references{
Pearl J, Bareinboim E. External validity: From do-calculus to transportability across populations. Statistical Science. 2012;29(4):579-595.
}

\examples{
### For a binary outcome
data("dataPROPHYVAP")

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

.f <- formula(VAP_num ~ GROUP_num * AGE + SEX + BMI + DIABETES)

### In practice use larger values of boot.number (e.g., 500)
### We set boot.number and cv at 10 for speed in CRAN checks
gc1 <- gc_binary(formula=.f, model="ridge", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", progress=TRUE, boot.tune=FALSE)

summary(gc1)

### New dataset (for example a subset of younger patients)
newdata <- subset(dataPROPHYVAP, AGE <= 50)

### transport the gc_binary model to the new dataset
gc2 <- transport(object=gc1, newdata=newdata, effect="ATE", boot.number=10)

summary(gc2)
}

\keyword{Transportation}