\name{transpose}
\alias{transpose}
\title{
Transposition for G-computation Logistic and Survival Models
}
\description{
Applies an already fitted G-computation model (an object of class \code{gclogi} or \code{gcsurv}) to a new dataset (\code{newdata}) to estimate marginal effects for that new population.
}
\usage{
transpose(object, newdata, effect="ATE", boot.number=50000, seed=NULL)
}
\arguments{
  \item{object}{
An object returned by the function \code{\link{gc_logistic}} or \code{\link{gc_survival}}.
}
  \item{newdata}{
A data frame for the new dataset to which the G-computation model will be applied. This dataset should contain all covariates present in the \code{formula} used to fit the original \code{object}.
}
  \item{effect}{
The type of the marginal effect to be estimated on the \code{newdata}. Three types are possible: "ATE" (by default), "ATT", and "ATU". These correspond to Average Treatment Effect on the Entire population, on the Treated, and on the Untreated, respectively. See details in \code{\link{gc_logistic}} or \code{\link{gc_survival}} for more information.
}
  \item{boot.number}{
The number of bootstrap resamples to perform for estimating standard errors and confidence intervals of the marginal effects on the \code{newdata}. The default value is 500.
}
  \item{seed}{
A random seed to ensure reproducibility during the bootstrap process. If \code{NULL} (the default), a seed is randomly assigned.
}
}
\details{
The \code{transpose} function allows for applying a previously fitted G-computation model to a different population specified by \code{newdata}. The Q-model (logistic or survival) from the original \code{object}, calculated on the previous entire population, is used to predict counterfactual outcomes for each individual in \code{newdata}, as if they were exposed to both treatment (group = 1) and control (group = 0) conditions. These predictions are then used to estimate marginal effects in \code{newdata}.

The type of bootstrap ("bcv" or "boot") used for the variance estimation on \code{newdata} is inherited from the original \code{object} (i.e., \code{object$boot.type}).
}
\value{
The function returns an object of the same class as the input \code{object} (\code{gclogi} or \code{gcsurv}), but with results updated for the \code{newdata}. The returned object is a list containing the same arguments as in \code{\link{gc_logistic}} or \code{\link{gc_survival}}.
}
\references{
Joe de Keizer et al. G-computation for increasing performances of clinical trials with individual randomization and binary response. ArXiv 2024-11.
Chatton et al. G-computation and doubly robust standardisation for continuous-time data: A comparison with inverse probability weighting. Stat Methods Med Res. 31(4):706-718. 2022.
}

\examples{
### For a binary outcome
data("dataPROPHYVAP")

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)


.f <- formula(VAP_num ~ GROUP_num * (AGE + SEX + BMI + DIABETES))

gc1 <- gc_logistic(formula=.f, method="ridge", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=500,  effect="ATE", progress=TRUE, boot.tune = TRUE)

print(gc1)
summary(gc1)
plot(gc1)

### New dataset (for example a subset of younger patients)
newdata_logistic <- subset(dataPROPHYVAP, AGE <= 50)

### Transpose the gc_logistic model to the new dataset
gclogi_transposed <- transpose(object = gc1,
                                newdata = newdata_logistic,
                                effect = "ATE",
                                boot.number = 5)

print(gclogi_transposed)
summary(gclogi_transposed, ci.type = "norm")


### For continuous-time data
data("dataPROPHYVAP")

dataPROPHYVAP$DEATH_num <- ifelse(dataPROPHYVAP$DEATH == "Yes", 1, 0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo", 0, 1)

.f <- formula(Surv(TIME_DEATH, DEATH_num) ~ GROUP_num * (AGE + SEX + BMI + DIABETES))

gcsurv1 <- gc_survival(formula = .f, method = "ridge",
                       data = dataPROPHYVAP, group = "GROUP_num",
                       pro.time = 60, param.tune = NULL, boot.type = "bcv", cv = 10,
                       boot.number = 5, effect = "ATE", progress = TRUE, boot.tune = TRUE)

print(gcsurv1)
summary(gcsurv1)
plot(gcsurv1)

### New dataset (for example a subset of younger patients)
newdata_survival <- subset(dataPROPHYVAP, AGE <= 50)

gcsurv_transposed <- transpose(object = gcsurv1,
                                newdata = newdata_survival,
                                effect = "ATE",
                                boot.number = 5,
                                seed = 4132)

print(gcsurv_transposed)
summary(gcsurv_transposed, ci.type = "perc")
}
\keyword{Causal inference}
\keyword{G-computation}
\keyword{transposition}
\keyword{prediction}