\name{gc_times}
\alias{gc_times}

\title{
G-Computation to Estimate a Marginal Effect for Time-to-Event Outcome
}
\description{
This function computes G-computation (GC) with different working models or algorithms (Q-models) for a time-to-event outcome and a 2-class exposure/treatment.
}
\usage{
gc_times(formula, data, group, pro.time=NULL, effect="ATE", model,
            param.tune=NULL, cv=30, boot.type="bcv", boot.number=500,
            boot.tune=FALSE, progress=TRUE, seed=NULL, boot.mi=FALSE, m=5, ...)
}

\arguments{
  \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame in which to look for the variables related to the outcome, the studied (\code{group}) and the predictors included in the previous model.}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones.}
   \item{pro.time}{An optional value for censoring the follow-up times and obtaining survival curves and related restricted mean survival times up to \code{pro.time}. If not specified, it is automatically set to the time point at which 10\% of patients remain at risk.}
   \item{effect}{The type of the marginal effect to be estimated. Three types are possible : "ATE" (by default), "ATT" and "ATU". See details.}   
   \item{model}{The modelling method used to create the Q-model. Current implemented methods are : "all", "lasso", "ridge", "elasticnet" and "aic". See details.} 
   \item{param.tune}{An optional argument to specify the tuning parameters for the Q-model. If \code{NULL} (the default), the tuning parameter(s) is(are) estimated by \code{cv}-fold cross-validation. Otherwise, the user can propose a tuning grid. See details.}
    \item{cv}{The number of splits for cross-validation. The default value is 30.}
    \item{boot.type}{The type of bootstrap to use. Two types are possible: "bcv" (by default) and "boot". See details.}
    \item{boot.number}{The number of bootstrap resamples. The default value is 500.}
    \item{boot.tune}{A logical value to determine whether the tuning parameter should be estimated inside of each bootstrap iteration. See details.}
    \item{progress}{A logical value to print a progress bar in the R console. The default is \code{TRUE}}
    \item{seed}{A random seed to ensure reproducibility during the cv process. If \code{NULL}, a seed is randomly assigned.}
    \item{boot.mi}{A logical value to apply multiple imputation using the \code{\link[mice]{mice}} package. See details.}
    \item{m}{Number of multiple imputations to perform if boot.mi is \code{TRUE}.}
  \item{...}{Additional arguments to be passed directly to the \code{\link[mice]{mice}} function for customizing the multiple imputation process (e.g., \code{method}, \code{maxit}, \code{diagnostics}).}
 }

\details{
The option \code{effect="ATE"} corresponds to the Average Treatment effect on the Entire population, i.e. the marginal effect if the entire sample were treated versus untreated. The "ATT" modality allows the estimation of the Average Treatment effect on the Treated, i.e. the marginal effect if the treated patients (\code{group = 1}) would have been untreated. The "ATU" modality modality allows the estimation of the Average Treatment effect on the Untreated, i.e. the marginal effect if the untreated patients (\code{group = 0}) would have been treated.

Several modelling methods can be used for the Q-model estimation:
\tabular{lll}{
\code{"all"} \tab A proportional hazard regression, the baseline hazard function being estimated by using the Breslow estimator and all the covariates in the formula being used.\cr
\code{"lasso"} \tab The same PH regression with L1 regularization allowing predictors' selection.\cr
\code{"ridge"} \tab The same PH regression with L2 regularization.\cr
\code{"elasticnet"} \tab The same PH regression which combines both L1 and L2 regularizations.\cr
\code{"aic"} \tab The same PH regression with a AIC-based forward selection.\cr
\code{"bic"} \tab The same PH regression with a BIC-based forward selection.
}

The \code{param.tune} argument allows users to specify tuning parameters penalized regression. If \code{NULL} (the default), the tuning parameters of each algorithm are estimated by cv-fold cross-validation and the default grids on the \code{glmnet} package. Otherwise, the user can propose a tuning grid for each modelling method. For "lasso" and "ridge" penalizations, it should be a vector representing the lambda penalization parameter. For "elasticnet", it should be a list or a vector of length 2, containing lambda (penalization parameter) and alpha (mixing parameter between L1 and L2 regularizations) values. The alpha value typically ranges from 0 to 1. Each object in the list declared in param.tune must have the same name as the corresponding parameter of the method. The user may propose a single value of each tuning parameter if she/he aims to define her/his own penality.

The \code{boot.tune} argument is logical value which determines whether the tuning parameter should be estimated inside of each bootstrap iteration. If \code{FALSE} (the default), the tuning parameter will be estimated once on the complete dataset.

If \code{boot.mi = TRUE}, multiple imputation is considered by using the MI-BOOT approach proposed by Schomaker & Heumann (2018). The input \code{data} is first imputed with the \code{\link[mice]{mice}} function to create \code{m} complete datasets. GC is then performed on each imputed dataset with the specified Q-model, and bootstrap estimates from all runs are concatenated. Note that all columns in the provided dataset are used as predictors during the imputation process.
}

\value{
\item{calibration}{A list containing the fitted Q-model (\code{fit}); the time points (\code{time}) with the related cumulative hazards (\code{cumhaz}), survival probabilities (\code{surv}), and baseline cumulative hazard (\code{H0.multi}); the linear predictors (\code{lp}) obtained on the complete dataset; survival probabilities for group 0 (\code{surv0}); and survival probabilities for group 1 (\code{surv1})}
\item{tuning.parameters}{The estimated tuning parameters for the Q-model. For "aic" and "bic" methods, this represents the final model.}
\item{data}{The data frame with individual with no missing data in the formula parameters.}
\item{formula}{The formula provided to the user.}
\item{model}{The method used for the Q-model (e.g., "lasso", "ridge", "aic").}
\item{cv}{The number of splits used for cross-validation.}
\item{missing}{The number of observations that were removed from the original dataset due to missing values.}
\item{pro.time}{The time point up to which RMST and survival probabilities are estimated.}
\item{boot.number}{The total number of bootstrap resamples performed.}
\item{boot.type}{The type of bootstrap resamples performed.}
\item{group}{A character string specifying the name of the grouping variable.}
\item{n}{The sample size of the dataset used in calculations, after missing data removal.}
\item{nevent}{The total number of events in the dataset used in calculations.}
\item{adjusted.results}{A data frame containing the adjusted results for each bootstrap sample, including: \code{AHR} (average hazard ratio), \code{RMST0} and \code{RMST1} (restricted mean survival times), \code{deltaRMST} (difference in RMST), \code{s0} and \code{s1} (survival probabilities at pro.time), and \code{delta} (difference in survival probabilities).}
\item{unadjusted.results}{A data frame containing the unadjusted results for each bootstrap sample, including: \code{AHR} (average hazard ratio), \code{RMST0} and \code{RMST1} (restricted mean survival times), \code{deltaRMST} (difference in RMST), \code{s0} and \code{s1} (survival probabilities at pro.time), and \code{delta} (difference in survival probabilities).}
\item{call}{The complete function call that generated the \code{gctimes} object.}
}

\references{
Chatton et al. G-computation and doubly robust standardisation for continuous-time data: A comparison with inverse probability weighting. Stat Methods Med Res. 31(4):706-718. 2022. <doi:10.1177/09622802211047345>. 
}

\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$DEATH_num <- ifelse(dataPROPHYVAP$DEATH == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

.f <- formula(Surv(TIME_DEATH, DEATH_num) ~ GROUP_num + AGE +
              SEX + BMI + DIABETES)

### In practice use larger values of boot.number (e.g., 500)
### We set boot.number at 10 for speed in CRAN checks
gc1 <- gc_times(formula=.f, model="lasso", data=dataPROPHYVAP,
              group="GROUP_num", boot.type="bcv",
              boot.number=10,  effect="ATE", progress=TRUE , pro.time=10,
              boot.tune=FALSE)

summary(gc1)
}

\keyword{G-computation}
\keyword{Time-to-event outcome}