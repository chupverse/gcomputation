\name{transport}
\alias{transport}
\title{
Transportability for G-computation binary and times Models
}
\description{
Applies an already fitted G-computation model (an object of class \code{gcbinary} or \code{gctimes}) to a new dataset (\code{newdata}) to estimate marginal effects for that new population.
}
\usage{
transport(object, newdata, effect="ATE", boot.number=500, seed=NULL)
}
\arguments{
  \item{object}{
An object returned by the function \code{\link{gc_binary}} or \code{\link{gc_times}}.
}
  \item{newdata}{
A data frame for the new dataset to which the G-computation model will be applied. This dataset should contain all covariates present in the \code{formula} used to fit the original \code{object}.
}
  \item{effect}{
The type of the marginal effect to be estimated on the \code{newdata}. Three types are possible: "ATE" (by default), "ATT", and "ATU". These correspond to Average Treatment Effect on the Entire population, on the Treated, and on the Untreated, respectively. See details in \code{\link{gc_binary}} or \code{\link{gc_times}} for more information.
}
  \item{boot.number}{
The number of bootstrap resamples to perform for estimating standard errors and confidence intervals of the marginal effects on the \code{newdata}. The default value is 500.
}
  \item{seed}{
A random seed to ensure reproducibility during the bootstrap process. If \code{NULL} (the default), a seed is randomly assigned.
}
}
\details{
The \code{transport} function allows for applying a previously fitted G-computation model to a different population specified by \code{newdata}. The Q-model (binary or times) from the original \code{object}, calculated on the previous entire population, is used to predict counterfactual outcomes for each individual in \code{newdata}, as if they were exposed to both treatment (group = 1) and control (group = 0) conditions. These predictions are then used to estimate marginal effects in \code{newdata}.

The type of bootstrap ("bcv" or "boot") used for the variance estimation on \code{newdata} is inherited from the original \code{object} (i.e., \code{object$boot.type}).
}
\value{
The function returns an object of the same class as the input \code{object} (\code{gcbinary} or \code{gctimes}), but with results updated for the \code{newdata}. The returned object is a list containing the same arguments as in \code{\link{gc_binary}} or \code{\link{gc_times}}.
}
\references{
Joe de Keizer et al. G-computation for increasing performances of clinical trials with individual randomization and binary response. ArXiv 2024-11.

Chatton et al. G-computation and doubly robust standardisation for continuous-time data: A comparison with inverse probability weighting. Stat Methods Med Res. 31(4):706-718. 2022.

Pearl J, Bareinboim E. External validity: From do-calculus to transportability across populations. Statistical Science. 2012;29(4):579-595.

Colnet B, Mayer I, Chen G, Dieng A, Li R, Varoquaux G, Vert J-P, Josse J, Yang S. Causal Inference Methods for Combining Randomized Trials and Observational Studies: A Review. Statistical Science. 2024;39(1):79-104. doi:10.1214/23-STS889.
}

\examples{
### For a binary outcome
data("dataPROPHYVAP")

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

.f <- formula(VAP_num ~ GROUP_num * (AGE + SEX + BMI + DIABETES))

### In practice use larger values of boot.number (e.g., 500).
### Here we set boot.number = 10 for speed in CRAN checks.
gc1 <- gc_binary(formula=.f, model="ridge", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=10, effect="ATE", progress=TRUE, boot.tune=FALSE,
              seed=7718)

print(gc1)
summary(gc1)
plot(gc1)

### New dataset (for example a subset of younger patients)
newdata_binary <- subset(dataPROPHYVAP, AGE <= 50)

### transport the gc_binary model to the new dataset
gc1_tr <- transport(object=gc1, newdata=newdata_binary,
                                effect="ATE",boot.number=10)

print(gc1_tr)
summary(gc1_tr, ci.type = "norm")


### For time-to-event data
data("dataPROPHYVAP")

dataPROPHYVAP$DEATH_num <- ifelse(dataPROPHYVAP$DEATH == "Yes", 1, 0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo", 0, 1)

.f <- formula(Surv(TIME_DEATH, DEATH_num) ~ GROUP_num * (AGE + SEX + BMI + DIABETES))

### boot.number is set to 10 here only for speed.
gc2 <- gc_times(formula=.f, model="ridge",
                       data=dataPROPHYVAP, group="GROUP_num",
                       pro.time=60, param.tune=NULL, boot.type="bcv", cv=10,
                       boot.number=10, effect="ATE", progress=TRUE, boot.tune=FALSE,
                       seed=9341)

print(gc2)
summary(gc2)
plot(gc2)

### New dataset (for example a subset of younger patients)
newdata_times <- subset(dataPROPHYVAP, AGE <= 50)

gc2_tr <- transport(object=gc2, newdata=newdata_times,
                                effect="ATE", boot.number=10, seed=4132)

print(gc2_tr)
summary(gc2_tr, ci.type="perc")
}
\keyword{Causal inference}
\keyword{G-computation}
\keyword{transposition}
\keyword{prediction}