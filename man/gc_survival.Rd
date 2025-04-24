\name{gc_survival}
\alias{gc_survival}

\title{
G-Computation to Estimate a Marginal Effect
}
\description{
This function computes G-Computation (GC) with different working models or algorithms (Q-models)to estimate the marginal survival functions reagrding two exposures.
}
\usage{
gc_survival(formula, data, group, max.time, effect="ATE", method, param.tune=NULL, cv=10, boot.type=NULL,
            boot.number=500,  boot.tune=FALSE, progress=TRUE)
}

\arguments{
  \item{formula}{A survival formula related to the Q-model with the the variable \code{group} among the predictors.}
  \item{data}{A data frame in which to look for the variables related to the status of the follow-up time, the event , the studied (\code{group}) and the covariables included in the previous model.}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones.}
   \item{max.time}{An optional value for censoring the follow-up times and obtaining survival curves and related restricted mean survival times up to \code{max.time}. The default value is median of the follow-up times.}
   \item{effect}{To be completed BY JOE.}   
   \item{method}{To be completed BY JOE.} 
    \item{param.tune}{To be completed BY JOE. QQ PHRASES SYNTHETIQUES. See details.}
    \item{cv}{The number of splits for cross-validation. The default value is 10.}
    \item{boot.type}{To be completed BY JOE.}
     \item{boot.number}{The number of bootstrap resamples. The default value is 500.}
     \item{boot.tune}{To be completed BY JOE.}  
   \item{progress}{A logical value to print a progress bar in the R console. The default is \code{TRUE}}
 }

\details{
TO BE COMPLETED BY JOE. Each object of the list declared in \code{param.tune} must have the same name than the names of the \code{methods} included in the SL. If \code{param.tune} = \code{NULL}, the tunning parameters of each algorithm are estimated by \code{cv}-fold cross-validation. Otherwise, the user can propose a tunning grid for each method, as explained in the following table. The following metrics can be used: "ci" for the concordance index at the prognostic time \code{pro.time}, "bs" for the Brier score at the prognostic time \code{pro.time}, "loglik" for the log-likelihood, "ibs" for the integrated Brier score up to the last observed time of event, "ibll" for the Integrated binomial log-likelihood up to the last observed time of event, "bll" for the binomial log-likelihood, "ribs" for the restricted integrated Brier score up to the prognostic time \code{pro.time}, "ribll" for the restricted integrated binomial log-likelihood up to the last observed time of event, "bll" for the binomial log-likelihood, and "auc" for the area under the time-dependent ROC curve up to the prognostic time \code{pro.time}.


The following learners are available:
  \tabular{llll}{
  Names \tab Description \tab Package  \cr
  \code{"LIB_AFTgamma"} \tab Gamma-distributed AFT model \tab flexsurv  \cr
  \code{"LIB_AFTggamma"} \tab  Generalized Gamma-distributed AFT model \tab flexsurv  \cr
  \code{"LIB_AFTweibull"} \tab  Weibull-distributed AFT model \tab flexsurv  \cr
  \code{"LIB_PHexponential"} \tab  Exponential-distributed PH model \tab flexsurv  \cr
  \code{"LIB_PHgompertz"} \tab  Gompertz-distributed PH model \tab flexsurv  \cr
  \code{"LIB_PHspline"} \tab  Spline-based PH model \tab flexsurv  \cr
  \code{"LIB_COXall"} \tab  Usual Cox model \tab survival  \cr
  \code{"LIB_COXaic"} \tab  Cox model with AIC-based forward selection\tab MASS  \cr
  \code{"LIB_COXen"} \tab  Elastic Net Cox model \tab glmnet  \cr
  \code{"LIB_COXlasso"} \tab  Lasso Cox model \tab glmnet  \cr
  \code{"LIB_COXridge"} \tab   Ridge Cox model \tab glmnet  \cr
  \code{"LIB_RSF"} \tab Survival Random Forest \tab
  randomForestSRC  \cr
  \code{"LIB_SNN"} \tab (Python-based) Survival Neural Network \tab survivalmodels \cr
  \code{"LIB_PLANN"} \tab (Python-based) Survival Neural Network \tab survivalPLANN \cr}

The following loss functions for the estimation of the super learner weigths are available (\code{metric}):
\itemize{
  \item Area under the ROC curve (\code{"auc"})
  \item Concordance index (\code{"ci"})
  \item Brier score (\code{"bs"})
  \item Binomial log-likelihood (\code{"bll"})
  \item Integrated Brier score (\code{"ibs"})
  \item Integrated binomial log-likelihood (\code{"ibll"})
  \item Restricted integrated Brier score (\code{"ribs"})
  \item Restricted integrated binomial log-Likelihood (\code{"ribll"})
}
}

\value{
\item{calibration}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A list of matrices with the predictions of survivals of each subject (lines) for each observed time (columns). Each matrix corresponds to the included \code{methods} and the resulted SL (the last item entitled "sl"). If \code{keep.predictions=TRUE}, it corresponds to a matrix with predictions related to the SL.}
\item{data}{The data frame used for learning. The first column is entitled \code{times} and corresponds to the observed follow-up times. The second column is entitled \code{failures} and corresponds to the event indicators. The other columns correspond to the predictors.}
\item{predictors}{A list with the predictors involved in \code{group}, \code{cov.quanti} and \code{cov.quali}.}
\item{ROC.precision}{The percentiles (between 0 and 1) of the prognostic variable used for computing each point of the time dependent ROC curve.}
\item{cv}{The number of splits for cross-validation.}
\item{pro.time}{The maximum delay for which the capacity of the variable is evaluated.}
\item{models}{A list with the estimated models/algorithms included in the SL.}
\item{weights}{A list composed by two vectors: the regressions \code{coefficients} of the logistic multinomial regression and the resulting weights' \code{values}}
\item{metric}{A list composed by two vectors: the loss function used to estimate the weights of the algorithms in the SL and its value.}
\item{param.tune}{The estimated tunning parameters.}
}

\references{
Chatton et al. G-computation and doubly robust standardisation for continuous-time data: A comparison with inverse probability weighting. Stat Methods Med Res. 31(4):706-718. 2022. <doi: 10.1177/09622802211047345>. 
}

\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$DEATH_J60_num <- ifelse(dataPROPHYVAP$DEATH_J60 == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

formula <- formula(Surv(FOLLOWUP_J60, DEATH_J60_num) ~ GROUP_num * AGE + SEXE + BMI + DIABETES)

gc1 <- gc_survival(formula=formula, method="lasso", data=dataPROPHYVAP, group="GROUP_num",
              param.tune=NULL, boot.type="bcv", cv=10, boot.number=50,  effect="ATE",
              progress=TRUE , max.time=10, boot.tune = TRUE)

print(gc1)
}

\keyword{G-Computation}
\keyword{Survival}