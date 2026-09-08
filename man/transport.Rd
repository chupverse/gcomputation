\name{transport}
\alias{transport}
\title{Transport the Marginal Effect for Another Populations}

\description{
  Applies an already fitted Q-model to estimate marginal effects for a new population (i.e., with differentially distributed baseline characteristics used as predictors in the Q-Model).
}

\usage{
transport(object, newdata, estim_var, nboot=1000, n.sim=500, seed=NULL)
}

\arguments{
  \item{object}{An object returned by the function \code{\link{gc_binary}},  \code{\link{gc_continuous}},  \code{\link{gc_count}} or \code{\link{gc_times}}.}
  \item{newdata}{The new data frame with all the covariates in the \code{formula} of the \code{object}. Rows containing missing values are
removed.}
  \item{estim_var}{Method used to estimate the transported effect and its  uncertainty. Four options are available: \code{"monte_carlo"} for Monte Carlo simulation, \code{"m_estimation"} for M-estimation, \code{"bootstrap"} for bootstrap, and \code{"point_estimate"} for obtaining the point estimate without estimating its variance.}
  \item{nboot}{The number of bootstrap resamples. The default is 1000. This is
only applicable when \code{estim_var = "bootstrap"}.}
  \item{n.sim}{The number of Monte Carlo simulations to perform for estimating standard errors and confidence intervals of the marginal effects on the \code{newdata}. This is only applicable for unpenalized (likelihood-based) models. The default value is 500.}
  \item{seed}{A random seed to ensure reproducibility during the simulation process. If \code{NULL} (the default), a seed is randomly assigned.}
}

\details{
The function employs three distinct estimation strategies, depending on
the value of the \code{estim_var} argument:
\itemize{
\item \bold{\code{estim_var = "monte_carlo"} (Monte Carlo simulation):}
For unpenalized models (likelihood-based models such as \code{"all"},
\code{"aic"}, or \code{"bic"}), the function uses the variance-covariance
matrix of the Q-model regression coefficients to simulate Q-models
assuming a multivariate normal distribution. This provides a distribution
of adjusted results and the related confidence intervals.

\item \bold{\code{estim_var = "point_estimate"} (Point estimate):}
For penalized models (\code{"lasso"}, \code{"ridge"},
\code{"elasticnet"}) and survival models (\code{gctimes}),
the function only computes the transported point estimate.
Standard errors and confidence intervals are not estimated.

\item \bold{\code{estim_var = "m_estimation"} (M-estimation):}
For unpenalized models (likelihood-based) models such as \code{"all"},
\code{"aic"}, or \code{"bic"}), the function uses estimating equations,
followed by root solving and an empirical sandwich estimator to estimate
standard errors. Root solving and sandwich estimation are performed
using the \pkg{geex} package.

\item \bold{\code{estim_var = "bootstrap"} (Bootstrap):}
Bootstrap resampling is applied on both populations. 
At each bootstrap iteration, the original dataset is resampled with replacement
and a new outcome model is fitted to this bootstrap sample using the initial
parameters. The new population is also resampled with replacement, and the
fitted model is used to estimate the marginal effects in this bootstrap new
population. The 2.5th and 97.5th percentiles of the bootstrap estimates are
then used for the 95\% confidence interval. Compared to the previous two
methods, this one adds sampling variability in both the original and
new populations.
}
Transportation cannot be applied to an object that has already been
transported, and it is also not implemented when multiple imputation has been
used to create the original object.
}

\value{
The function returns an object of the same class with the same arguments as the input \code{object}, but with updated results for the \code{newdata}.
}
\references{
Pearl J, Bareinboim E. External validity: From do-calculus to transportability across populations. Statistical Science. 2012;29(4):579-595.

Ross RK, Zivich PN, Stringer JSA, Cole SR. M-estimation for common epidemiological measures: introduction and applied examples. Int J Epidemiol. 2024 Feb 14;53(2):dyae030. doi: 10.1093/ije/dyae030. PMID: 38423105; PMCID: PMC10904145.

Saul BC, Hudgens MG . The Calculus of M-Estimation in R with geex. Journal of Statistical Software. 2020;92(2), 1–15. doi:10.18637/jss.v092.i02.
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
newdata <- na.omit(newdata)

### transport the gc_binary model to the new dataset

gc2 <- transport(object=gc1, newdata=newdata, estim_var = "monte_carlo", n.sim=10)
summary(gc2, ci.type="perc", ci.level=0.95)


gc3 <- transport(object=gc1, newdata=newdata, estim_var = "m_estimation")
summary(gc3, ci.type="norm", ci.level=0.95)

}

\keyword{Transportation}