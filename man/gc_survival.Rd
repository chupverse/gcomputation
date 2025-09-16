\name{gc_survival}
\alias{gc_survival}

\title{
G-Computation to Estimate a Marginal Effect
}
\description{
This function computes G-computation (GC) with different working models or algorithms (Q-models) to estimate the marginal survival functions regarding two exposures.
}
\usage{
gc_survival(formula, data, group, pro.time, effect="ATE", method,
            param.tune=NULL, cv=10, boot.type="bcv", boot.number=50000,
            boot.tune=FALSE, progress=TRUE, seed=NULL)
}

\arguments{
  \item{formula}{A survival formula related to the Q-model with the the variable \code{group} among the predictors.}
  \item{data}{A data frame in which to look for the variables related to the status of the follow-up time, the event , the studied (\code{group}) and the covariables included in the previous model.}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones.}
   \item{pro.time}{An optional value for censoring the follow-up times and obtaining survival curves and related restricted mean survival times up to \code{pro.time}. The default value is median of the follow-up times.}
   \item{effect}{The type of the marginal effect to be estimated. Three types are possible : "ATE" (by default), "ATT" and "ATU". See details.}   
   \item{method}{The method used to create the Q-model. Current implemented methods are : "all", "lasso", "ridge", "elasticnet", "aic" and "bic". See details.} 
    \item{param.tune}{An optional argument to specify the tuning parameters for the Q-model.
If \code{NULL} (the default), the tuning parameters of each algorithm are estimated by \code{cv}-fold cross-validation.
Otherwise, the user can propose a tuning grid or specific values to use for each method. See details.}
    \item{cv}{The number of splits for cross-validation. The default value is 10.}
    \item{boot.type}{The type of bootstrap to use. Two types are possible: "bcv" (by default) and "boot". See details.}
     \item{boot.number}{The number of bootstrap resamples. The default value is 500.}
     \item{boot.tune}{A logical value to estimate the tuning parameter inside of each bootstrap iteration, if \code{FALSE} it will estimate it once on the complete data set. The default is \code{FALSE}}  
   \item{progress}{A logical value to print a progress bar in the R console. The default is \code{TRUE}}
    \item{seed}{A random seed to ensure reproducibility during the cv process. If \code{NULL}, a seed is randomly assigned.}
 }

\details{
The ATE corresponds to Average Treatment effect on the Entire population, i.e. the marginal effect if all the sample is treated versus all the sample is untreated. The ATT corresponds to Average Treatment effect on the Treated, i.e. the marginal effect if the treated patients (\code{group = 1}) would have been untreated. The ATU corresponds to Average Treatment effect on the Untreated , i.e. the marginal effect if the untreated patients (\code{group = 0}) would have been treated.

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
\item{calibration}{A list containing the final fitted Q-model (\code{fit}), the predicted time points (\code{time}), cumulative hazards (\code{cumhaz}), survival probabilities (\code{surv}), baseline cumulative hazard (\code{H0.multi}), and linear predictors (\code{lp}) obtained on the complete dataset.}
\item{tuning.parameters}{The estimated tuning parameters for the Q-model. For "aic" and "bic" methods, this is equal to \code{NULL}. Else it is a list containing the optimal lambda and alpha values.}
\item{data}{The data frame used in the function, after any rows with missing data have been removed.}
\item{formula}{The original formula provided to the function.}
\item{method}{The method used for the Q-model (e.g., "lasso", "ridge", "aic").}
\item{cv}{The number of splits used for cross-validation.}
\item{missing}{The number of observations that were removed from the original dataset due to missing values.}
\item{pro.time}{The specific time point up to which RMST and survival probabilities are estimated.}
\item{boot.number}{The total number of bootstrap resamples performed.}
\item{boot.type}{The type of bootstrap resamples performed.}
\item{outcome}{A list of character strings, specifying the names of the time (\code{times}) and outcome  variable.}
\item{group}{A character string specifying the name of the grouping (exposure/treatment) variable.}
\item{n}{The sample size of the dataset used in calculations, after missing data removal.}
\item{nevent}{The total number of events (outcome coded as 1) in the dataset used in calculations.}
\item{AHR}{A numeric vector containing the raw bootstrap samples of the estimated Average Hazard Ratio.}
\item{RMST0}{A numeric vector containing the raw bootstrap samples of the estimated Restricted Mean Survival Time for the control group, up to \code{pro.time}.}
\item{RMST1}{A numeric vector containing the raw bootstrap samples of the estimated Restricted Mean Survival Time for the experimental group, up to \code{pro.time}.}
\item{deltaRMST}{A numeric vector containing the raw bootstrap samples of the estimated difference in Restricted Mean Survival Times (\code{RMST1} - \code{RMST0}).}
\item{surv0}{A numeric vector containing the raw bootstrap samples of the estimated survival probability for the control group, at \code{pro.time}.}
\item{surv1}{A numeric vector containing the raw bootstrap samples of the estimated survival probability for the experimental group, at \code{pro.time}.}
\item{deltasurv}{A numeric vector containing the raw bootstrap samples of the estimated difference in survival probabilities (\code{surv1} - \code{surv0}), at \code{pro.time}.}
\item{AHR.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted Average Hazard Ratio.}
\item{RMST0.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted Restricted Mean Survival Time for the control group, up to \code{pro.time}.}
\item{RMST1.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted Restricted Mean Survival Time for the experimental group, up to \code{pro.time}.}
\item{deltaRMST.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted difference in Restricted Mean Survival Times (\code{RMST1.unadj} - \code{RMST0.unadj}).}
\item{surv0.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted survival probability for the control group, at \code{pro.time}.}
\item{surv1.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted survival probability for the experimental group, at \code{pro.time}.}
\item{deltasurv.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted difference in survival probabilities (\code{surv1.unadj} - \code{surv0.unadj}), at \code{pro.time}.}
\item{call}{The complete function call that generated the \code{gcsurv} object.}
}

\references{
Chatton et al. G-computation and doubly robust standardisation for continuous-time data: A comparison with inverse probability weighting. Stat Methods Med Res. 31(4):706-718. 2022. <doi: 10.1177/09622802211047345>. 
}

\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$DEATH_num <- ifelse(dataPROPHYVAP$DEATH == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

formula <- formula(Surv(TIME_DEATH, DEATH_num) ~ GROUP_num * (AGE +
              SEX + BMI + DIABETES))

gc1 <- gc_survival(formula=formula, method="lasso", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=500,  effect="ATE", progress=TRUE , pro.time=10,
              boot.tune = TRUE)

print(gc1)
summary(gc1)
plot(gc1)


### To specify specific tuning parameters:
gc2 <- gc_survival(formula=formula, method="elasticnet", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=list(lambda=c(0.1,0.25,0.5),
              alpha=0.5), boot.type="bcv", cv=10, boot.number=500, effect="ATE",
              progress=TRUE , pro.time=10, boot.tune = FALSE)
              
gc2$tuning.parameters
}

\keyword{G-computation}
\keyword{Survival}