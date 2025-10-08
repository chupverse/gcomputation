\name{gc_continuous}
\alias{gc_continuous}
\title{
G-Computation to Estimate a Marginal Effect
}
\description{
This function computes G-computation (GC) with different working models or algorithms (Q-models) to estimate the marginal functions for continuous outcomes regarding two exposures.
}
\usage{
gc_continuous(formula, data, group, effect="ATE", model,
            param.tune=NULL, cv=10, boot.type="bcv", boot.number=500,
            boot.tune=FALSE, boot.mi=FALSE, progress=TRUE, seed=NULL, m=5, ...)
}

\arguments{
  \item{formula}{A regression formula related to the Q-model with the the variable \code{group} among the predictors.}
  \item{data}{A data frame in which to look for the variables related to the outcome, the studied (\code{group}) and the predictors included in the previous model.}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones.}
   \item{effect}{The type of the marginal effect to be estimated. Three types are possible : "ATE" (by default), "ATT" and "ATU". See details.}   
   \item{model}{The modelling method used to create the Q-model. Current implemented methods are : "all", "lasso", "ridge", "elasticnet" and "aic". See details.} 
   \item{param.tune}{An optional argument to specify the tuning parameters for the Q-model. If \code{NULL} (the default), the tuning parameters of each algorithm are estimated by \code{cv}-fold cross-validation. Otherwise, the user can propose a tuning grid or specific values to use for each method. See details.}
    \item{cv}{The number of splits for cross-validation. The default value is 10.}
    \item{boot.type}{The type of bootstrap to use. Two types are possible: "bcv" (by default) and "boot". See details.}
    \item{boot.number}{The number of bootstrap resamples. The default value is 500.}
    \item{boot.tune}{A logical value to determine whether the tuning parameter should be estimated inside of each bootstrap iteration. See details.}
    \item{boot.mi}{A logical value to apply multiple imputation using the \code{\link[mice]{mice}} package and then applying GC on each imputed dataset. See details.}
    \item{progress}{A logical value to print a progress bar in the R console. The default is \code{TRUE}}
    \item{seed}{A random seed to ensure reproducibility during the cv process. If \code{NULL}, a seed is randomly assigned.}
    \item{m}{Number of multiple imputations to perform if boot.mi is \code{TRUE}.}
  \item{...}{Additional arguments to be passed directly to the \code{\link[mice]{mice}} function for customizing the multiple imputation process (e.g., \code{method}, \code{maxit}, \code{diagnostics}).}
 }

\details{
The ATE (Average Treatment effect on the Entire population) corresponds to the marginal effect if the entire sample were treated versus entirely untreated. The ATT (Average Treatment effect on the Treated) is the marginal effect if the treated patients (group = 1) would have been untreated. The ATU (Average Treatment effect on the Untreated) represents the marginal effect if the untreated patients (group = 0) would have been treated.

For the group variable, it requires two modalities encoded as 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients). Similarly, the outcome variable must also have two modalities encoded as 0 (for censored patients/non-events) and 1 (for events).

Several modelling methods are implemented to create the Q-model:
\tabular{lll}{
\code{"all"} \tab Usual model \code{glm}, all given variables in the formula are used \cr
\code{"lasso"} \tab L1 regularization, allows predictor selection, considers interactions (A*X), B-splines for continuous covariates, uses the glmnet package, L1 norm penalization is the tuning parameter \cr
\code{"ridge"} \tab L2 regularization, uses L2 norm for penalization, equivalent to Elastic Net with an alpha value of 0, uses the glmnet package \cr
\code{"elasticnet"} \tab Combines both L1 and L2 regularizations, uses the glmnet package, alpha parameter mixes L1 and L2 (typically 0-1), default alpha sequence is seq(0,1,.1) if param.tune is NULL \cr
\code{"aic"} \tab Forward selection based on the Akaike Information Criterion (AIC), uses stepAIC with a penalty (k) of 2 \cr
\code{"bic"} \tab Forward selection based on the Bayesian Information Criterion (BIC), uses stepAIC with a penalty (k) of log(nrow(data)) \cr
}

The \code{param.tune} argument allows users to specify tuning parameters for the Q-model. If  \code{NULL} (the default), the tuning parameters of each algorithm are estimated by cv-fold cross-validation. Otherwise, the user can propose a tuning grid for each modelling method:
For "lasso" and "ridge" methods, it should be a scalar or a vector representing the lambda penalization parameter.
For "elasticnet", it should be a list or a vector of length 2, containing lambda (penalization parameter) and alpha (mixing parameter between L1 and L2 regularizations) values. The alpha value typically ranges from 0 to 1. Each object in the list declared in param.tune must have the same name as the corresponding parameter of the method.

Bootstrap Tuning \code{boot.tune}:
This logical value determines whether the tuning parameter should be estimated inside of each bootstrap iteration. If \code{FALSE} (the default), the tuning parameter will be estimated once on the complete dataset. If \code{TRUE} and only a single \code{lambda} value is provided for "lasso" or "ridge" methods, \code{boot.tune} will be automatically set to \code{FALSE} with a warning.

If \code{boot.mi = TRUE}, \code{gc_continuous} integrates multiple imputation with GC using the MI-BOOT approach (Schomaker & Heumann 2018). The input \code{data} is first imputed with the \code{\link[mice]{mice}} function to create \code{m} complete datasets. GC is then performed on each imputed dataset with the specified Q-model, and  bootstrap estimates from all runs are concatenated. 

The returned object remains of class \code{"gccontinuous"} and works with the usual S3 methods, but several fields differ 
compared to a standard call without multiple imputation:
\itemize{
  \item \code{data}: list of the \code{m} imputed datasets.
  \item \code{initial.data}: the original dataset before imputation.
  \item \code{calibration}, \code{tuning.parameters}: lists of results from each imputed dataset.
  \item \code{boot.number}: total number of bootstrap samples across imputations (\code{m} * \code{boot.number}).
  \item \code{m}: number of imputations used.
}
All bootstrap-based outcome measures (\code{m0}, \code{m1}, \code{delta}, \code{ratio}, and their unadjusted 
counterparts) are pooled by concatenating samples across imputations.
}

\value{
\item{calibration}{A list containing the final fitted Q-model (\code{fit}) and the outcome predictions (\code{predict}) for each patient on the complete dataset.}
\item{tuning.parameters}{The estimated tuning parameters for the Q-model. For "aic" and "bic" methods, this represents the final formula chosen; otherwise, it is a list containing the optimal lambda and alpha values.}
\item{data}{The data frame used in the function, after any rows with missing data have been removed.}
\item{formula}{The original formula provided to the function.}
\item{model}{The method used for the Q-model (e.g., "lasso", "ridge", "aic").}
\item{cv}{The number of splits used for cross-validation.}
\item{missing}{The number of observations that were removed from the original dataset due to missing values.}
\item{boot.number}{The total number of bootstrap resamples performed.}
\item{boot.type}{The type of bootstrap resamples performed.}
\item{outcome}{A character string specifying the name of the outcome variable.}
\item{group}{A character string specifying the name of the grouping (exposure/treatment) variable.}
\item{n}{The sample size of the dataset used in calculations, after missing data removal.}
\item{mean_outcome}{The mean of the outcome in the whole dataset used in calculations.}
\item{m0}{A numeric vector containing the raw bootstrap samples of the estimated outcome's means in the control group.}
\item{m1}{A numeric vector containing the raw bootstrap samples of the estimated outcome's means in the experimental group.}
\item{delta}{A numeric vector containing the raw bootstrap samples of the estimated difference in outcome's means (\code{m1} - \code{m0}).}
\item{ratio}{A numeric vector containing the raw bootstrap samples of the estimated ratio of outcome's means (\code{m1} / \code{m0}).}
\item{m0.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted outcome's means in the control group.}
\item{m1.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted outcome's means in the experimental group.}
\item{delta.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted difference in outcome's means (\code{p1.unadj} - \code{p0.unadj}).}
\item{ratio.unadj}{A numeric vector containing the raw bootstrap samples of the unadjusted ratio of outcome's means (\code{p1.unadj} / \code{p0.unadj}).}
\item{call}{The function call that generated the \code{gccontinuous} object.}
}

\references{
Joe de Keizer et al. G-computation for increasing performances of clinical trials with individual randomization and binary response. ArXiv 2024-11. <doi: 10.48550/arXiv.2411.10089>.

Schomaker M, Heumann C. Bootstrap Inference When Using Multiple Imputation. Statistics in medicine. 2018;37(14):2252-2266. doi: 10.1002/sim.7654.

Van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67.

}

\examples{
data("dataPROPHYVAP")

dataPROPHYVAP$VAP_num <- ifelse(dataPROPHYVAP$VAP == "Yes",1,0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo",0,1)

beta_ttt <- 10
beta_alcool <- 2
beta_bmi <- 8
beta_age <- 5

cont_outcome <- rnorm(n = nrow(dataPROPHYVAP), 
                      mean = beta_ttt*dataPROPHYVAP$GROUP_num +
                        beta_alcool*match(as.factor(dataPROPHYVAP$ALCOOL),levels(as.factor(dataPROPHYVAP$ALCOOL))) + 
                        beta_bmi*dataPROPHYVAP$BMI + 
                        beta_age*dataPROPHYVAP$AGE, 
                      sd = 1)

#mean(cont_outcome)

dataPROPHYVAP$y_cont <- cont_outcome

.f <- formula(y_cont ~ GROUP_num * (AGE + SEX + ALCOOL + BMI + DIABETES))

### In practice use larger values of boot.number (e.g., 500).
### Here we set boot.number = 10 for speed in CRAN checks.
gc1 <- gc_continuous(formula=.f, model="all", data=dataPROPHYVAP, group="GROUP_num",
                    cv=10, boot.type="boot", boot.number=10, boot.tune=TRUE,
                    effect="ATE", progress=TRUE, seed=5192,boot.mi = F)

print(gc1)
summary(gc1)
plot(gc1)


### To specify specific tuning parameters:
### boot.number is set to 10 here only for speed.
gc2 <- gc_continuous(formula=.f, model="elasticnet", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=list(lambda=0.1, alpha=0.2),
              boot.type="bcv", cv=10, boot.number=10, effect="ATE",
              progress=TRUE, seed=5321)
              
             
### To apply multiple imputation to handle missing data:
.f <- formula(VAP_num ~ GROUP_num + AGE + SEX + ALCOOL + BMI + DIABETES + GLASGOW + INJURY)

### boot.number is set to 10 here only for speed.
gc3 <- gc_continuous(formula=.f, model="all", data=dataPROPHYVAP,
                 group="GROUP_num", boot.type="boot", boot.number=10, effect="ATE",
                 progress=TRUE, boot.tune=FALSE, boot.mi=TRUE, m=5, seed=8812)
                 
plot(gc3,smooth=TRUE)
}

\keyword{G-computation}
\keyword{continuous}