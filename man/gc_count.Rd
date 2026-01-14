\name{gc_count}
\alias{gc_count}
\title{
G-Computation to Estimate a Marginal Effect with a Counting Outcome
}
\description{
This function computes G-computation (GC) with different working models or algorithms (Q-models) for a counting outcome and a 2-class exposure/treatment.
}
\usage{
gc_count(formula, data, group, effect="ATE", model,
            param.tune=NULL, cv=30, boot.type="bcv", boot.number=500,
            boot.tune=FALSE, progress=TRUE, seed=NULL, boot.mi=FALSE, m=5, ...)
}

\arguments{
  \item{formula}{A regression formula related to the Q-model with the the variable \code{group} among the predictors.}
  \item{data}{A data frame in which to look for the variables related to the outcome, the studied (\code{group}) and the predictors included in the previous model.}
  \item{group}{The name of the variable related to the exposure/treatment. This variable shall have only two modalities encoded 0 for the untreated/unexposed patients and 1 for the treated/exposed ones.}
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
The option \code{methods="ATE"} corresponds to the Average Treatment effect on the Entire population, i.e. the marginal effect if the entire sample were treated versus untreated. The "ATT" modality allows the estimation of the Average Treatment effect on the Treated, i.e. the marginal effect if the treated patients (\code{group = 1}) would have been untreated. The "ATU" modality modality allows the estimation of the Average Treatment effect on the Untreated, i.e. the marginal effect if the untreated patients (\code{group = 0}) would have been treated.

Several modelling methods can be used for the Q-model estimation:
\tabular{lll}{
\code{"all"} \tab A linear regression,  all the covariates in the formula being used.\cr
\code{"lasso"} \tab L1 regularized linear regression allows predictors' selection.\cr
\code{"ridge"} \tab L2 regularized, equivalent to Elastic Net with an alpha value of 0.\cr
\code{"elasticnet"} \tab A linear regression which combines both L1 and L2 regularizations.\cr
\code{"aic"} \tab A linear regression with a AIC-based forward selection.\cr
\code{"bic"} \tab A linear regression with a BIC-based forward selection.
}

The \code{param.tune} argument allows users to specify tuning parameters penalized regression. If \code{NULL} (the default), the tuning parameters of each algorithm are estimated by cv-fold cross-validation and the default grids on the \code{glmnet} package. Otherwise, the user can propose a tuning grid for each modelling method. For "lasso" and "ridge" penalizations, it should be a vector representing the lambda penalization parameter. For "elasticnet", it should be a list or a vector of length 2, containing lambda (penalization parameter) and alpha (mixing parameter between L1 and L2 regularizations) values. The alpha value typically ranges from 0 to 1. Each object in the list declared in param.tune must have the same name as the corresponding parameter of the method. The user may propose a single value of each tuning parameter if she/he aims to define her/his own penality.

The \code{boot.tune} argument is logical value which determines whether the tuning parameter should be estimated inside of each bootstrap iteration. If \code{FALSE} (the default), the tuning parameter will be estimated once on the complete dataset.

If \code{boot.mi = TRUE}, multiple imputation is considered by using the MI-BOOT approach proposed by Schomaker & Heumann (2018). The input \code{data} is first imputed with the \code{\link[mice]{mice}} function to create \code{m} complete datasets. GC is then performed on each imputed dataset with the specified Q-model, and bootstrap estimates from all runs are concatenated. Note that all columns in the provided dataset are used as predictors during the imputation process.
}

\value{
\item{calibration}{A list containing the fitted Q-model (\code{fit}) and the outcome predictions (\code{predict}).}
\item{tuning.parameters}{The estimated tuning parameters for the Q-model. For "aic" and "bic" methods, this represents the final model.}
\item{data}{The data frame with individual with no missing data in the formula parameters.}
\item{formula}{The formula provided to the user.}
\item{model}{The method used for the Q-model.}
\item{cv}{The number of splits used for cross-validation.}
\item{missing}{The number of observations that were removed from the original dataset due to missing values.}
\item{boot.number}{The total number of bootstrap resamples performed.}
\item{boot.type}{The type of bootstrap resamples performed.}
\item{group}{A character string specifying the name of the grouping variable.}
\item{n}{The sample size of the dataset used in calculations, after missing data removal.}
\item{adjusted.results}{A data frame containing the adjusted results for each bootstrap sample: \code{c1} (mean count experimental), \code{c0} (mean count control), \code{delta} (difference), and \code{ratio} (rate ratio).}
\item{unadjusted.results}{A data frame containing the unadjusted results for each bootstrap sample: \code{c1} (mean count experimental), \code{c0} (mean count control), \code{delta} (difference), and \code{ratio} (rate ratio).}
\item{call}{The function call that generated the \code{gccontinuous} object.}
}

\references{
Joe de Keizer et al. G-computation for increasing performances of clinical trials with individual randomization and binary response. ArXiv 2024-11. <doi:10.48550/arXiv.2411.10089>.

Schomaker M, Heumann C. Bootstrap inference when using multiple imputation. Statistics in medicine. 2018;37(14):2252-2266. <doi:10.1002/sim.7654>.}

\examples{
### Generate a dataset with 1000 individuals and a treatment effect of 2
sim_data <- data_gen_gcomp(n_indiv = 1000, beta_ttt = 2)

### Create formulas
.f<- formula(y_count ~ ttt + x1 + x2 + x3 + x4)

### Apply G-computation
gc1 <- gc_count(formula=.f, model="all", data=sim_data, group="ttt",
                 boot.type="boot", boot.number=10, seed=3210)
print(gc1)

#For count the C1/C0 ratio corresponds to exp(beta_ttt)
}

\keyword{G-computation}
\keyword{Counting outcome}