\name{gc_logistic}
\alias{gc_logistic}
\title{
G-Computation to Estimate a Marginal Effect
}
\description{
This function computes G-computation (GC) with different working models or algorithms (Q-models) to estimate the marginal logistic functions regarding two exposures.
}
\usage{
gc_logistic(formula, data, group, effect="ATE", method,
            param.tune=NULL, cv=10, boot.type="bcv", boot.number=500,
            boot.tune=FALSE, progress=TRUE)
}

\arguments{
  \item{formula}{A logistic regression formula related to the Q-model with the the variable \code{group} among the predictors.}
  \item{data}{A data frame in which to look for the variables related to the outcome, the studied (\code{group}) and the covariables included in the previous model.}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones.}
   \item{effect}{The type of the marginal effect to be estimated. Three types are possible : "ATE" (by default), "ATT" and "ATU". See details.}   
   \item{method}{The method used to create the Q-model. Current implemented methods are : "all", "lasso", "ridge", "elasticnet" and "aic". See details.} 
    \item{param.tune}{An optional argument to specify the tuning parameters for the Q-model.
If \code{NULL} (the default), the tuning parameters of each algorithm are estimated by \code{cv}-fold cross-validation.
Otherwise, the user can propose a tuning grid or specific values to use for each method. See details.}
    \item{cv}{The number of splits for cross-validation. The default value is 10.}
    \item{boot.type}{The type of bootstrap to use. Two types are possible: "bcv" (by default) and "boot". See details.}
     \item{boot.number}{The number of bootstrap resamples. The default value is 500.}
     \item{boot.tune}{A logical value to determine whether the tuning parameter should be estimated inside of each bootstrap iteration. See details.}


   \item{progress}{A logical value to print a progress bar in the R console. The default is \code{TRUE}}
 }

\details{
The ATE (Average Treatment effect on the Entire population) corresponds to the marginal effect if the entire sample were treated versus entirely untreated. The ATT (Average Treatment effect on the Treated) is the marginal effect if the treated patients (group = 1) would have been untreated. The ATU (Average Treatment effect on the Untreated) represents the marginal effect if the untreated patients (group = 0) would have been treated.

For the group variable, it requires two modalities encoded as 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients). Similarly, the outcome variable must also have two modalities encoded as 0 (for censored patients/non-events) and 1 (for events).

Several methods are implemented to create the Q-model:

\tabular{lll}{
\code{"all"} \tab Usual logistic/Cox model, all given variables in the formula are used \cr
\code{"lasso"} \tab L1 regularization for logistic/Cox regression, allows predictor selection, considers interactions (A*X), B-splines for continuous covariates, uses the glmnet package, L1 norm penalization is the tuning parameter \cr
\code{"ridge"} \tab L2 regularization for logistic/Cox regression, uses L2 norm for penalization, equivalent to Elastic Net with an alpha value of 0, uses the glmnet package \cr
\code{"elasticnet"} \tab Combines both L1 and L2 regularizations for logistic/Cox regression, uses the glmnet package, alpha parameter mixes L1 and L2 (typically 0-1), default alpha sequence is seq(0,1,.1) if param.tune is NULL \cr
\code{"aic"} \tab Forward selection for logistic/Cox regression based on the Akaike Information Criterion (AIC), uses stepAIC with a penalty (k) of 2 \cr
\code{"bic"} \tab Forward selection for logistic/Cox regression based on the Bayesian Information Criterion (BIC), uses stepAIC with a penalty (k) of log(nrow(data)) \cr
}

The  \code{param.tune} argument allows users to specify tuning parameters for the Q-model.
If  \code{NULL} (the default), the tuning parameters of each algorithm are estimated by cv-fold cross-validation.
Otherwise, the user can propose a tuning grid for each method:
For "lasso" and "ridge" methods, it should be a scalar or a vector representing the lambda penalization parameter.
For "elasticnet", it should be a list or a vector of length 2, containing lambda (penalization parameter) and alpha (mixing parameter between L1 and L2 regularizations) values. The alpha value typically ranges from 0 to 1.
Each object in the list declared in param.tune must have the same name as the corresponding parameter of the method.

Bootstrap Tuning \code{boot.tune}:
This logical value determines whether the tuning parameter should be estimated inside of each bootstrap iteration.
If \code{FALSE} (the default), the tuning parameter will be estimated once on the complete dataset.
If \code{TRUE} and only a single \code{lambda} value is provided for "lasso" or "ridge" methods, \code{boot.tune} will be automatically set to \code{FALSE} with a warning.
}

\value{
\item{calibration}{A list of predictions and the fit obtained on the whole data set.}
\item{tuning.parameters}{The estimated tunning parameters.}
\item{data}{The data frame used in the function.}
\item{formula}{The given formula.}
\item{method}{The given method.}
\item{cv}{The number of splits for cross-validation.}
\item{missing}{The number of observations with missing data excluded.}
\item{boot.number}{Number of bootstrap resamples.}
\item{outcome}{A character vector with the name of the outcome variable.}
\item{group}{A character vector with the name of the group variable.}
\item{n}{Sample size.}
\item{nevent}{Number of events.}
\item{coefficients}{All results of the proportions and their difference.}
\item{mOR}{All results of the marginal Odds Ratio.}
\item{call}{The complete call.}
}

\references{
Joe de Keizer et al. G-computation for increasing performances of clinical trials with individual randomization and binary response. ArXiv 2024-11. <doi: 10.48550/arXiv.2411.10089>.

}

\examples{
data("dataPROPHYVAP")

dataPROPHYVAP$VAP_EARLY_num <- ifelse(dataPROPHYVAP$VAP_EARLY == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)


formula <- formula(VAP_EARLY_num ~ GROUP_num * (AGE + SEXE + BMI + DIABETES))

gc1 <- gc_logistic(formula=formula, method="ridge", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=500,  effect="ATE", progress=TRUE, boot.tune = TRUE)

print(gc1)
summary(gc1)
plot(gc1)


### To specify specific tuning parameters:
gc2 <- gc_logistic(formula=formula, method="elasticnet", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=list(lambda=0.05, alpha=0.5),
              boot.type="bcv", cv=10, boot.number=500,  effect="ATE",
              progress=TRUE, boot.tune = FALSE)

}

\keyword{G-computation}
\keyword{logistic}