\name{transport}
\alias{transport}
\title{Transport the Marginal Effect for Another Populations}

\description{
  Applies an already fitted Q-model to estimate marginal effects for a new population (i.e., with differentially distibuted baseline characteristics used as predictors in the Q-Model).
}

\usage{
transport(object, newdata, n.sim=500, seed=NULL)
}

\arguments{
  \item{object}{An object returned by the function \code{\link{gc_binary}},  \code{\link{gc_continuous}},  \code{\link{gc_count}} or \code{\link{gc_times}}.}
  \item{newdata}{The new data frame with all the covariates in the \code{formula} of the \code{object}.}
  \item{n.sim}{The number of Monte Carlo simulations to perform for estimating standard errors and confidence intervals of the marginal effects on the \code{newdata}. This is only applicable for unpenalized models. The default value is 500.}
  \item{seed}{A random seed to ensure reproducibility during the simulation process. If \code{NULL} (the default), a seed is randomly assigned.}
}

\details{
The function employs two distinct estimation strategies:
\itemize{
  \item \bold{Monte Carlo Simulation:} For unpenalized models (likelihood-based models like "all", "aic", or "bic"), the function uses the variance-covariance matrix of the Q-model regression coefficients to simulate Q-models assuming a multivariate normal distribution. This provides a distribution of adjusted results and the related confidence intervals.
  \item \bold{Single Point Estimate:} For penalized models ("lasso", "ridge", "elasticnet") and survival models (\code{gctimes}), the function does not estimate the distribution of adjusted results and the related confidence intervals.
}
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

.f <- formula(VAP ~ GROUP * AGE + SEX + BMI + DIABETES)

### In practice use larger values of boot.number (e.g., 500)
### We set boot.number and cv at 10 for speed in CRAN checks
gc1 <- gc_binary(formula=.f, model="all", data=dataPROPHYVAP,
              group="GROUP", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", progress=TRUE, boot.tune=FALSE)

summary(gc1, ci.type="norm", ci.level=0.95)

### New dataset (for example a subset of younger patients)
newdata <- subset(dataPROPHYVAP, AGE <= 50)

### transport the gc_binary model to the new dataset
gc2 <- transport(object=gc1, newdata=newdata, n.sim=10)

summary(gc2, ci.type="perc", ci.level=0.95)
}

\keyword{Transportation}