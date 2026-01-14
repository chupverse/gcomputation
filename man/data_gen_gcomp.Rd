\name{data_gen_gcomp}
\alias{data_gen_gcomp}
\title{
Generate Simulated Data for G-Computation Analysis
}
\description{
Generates a simulated dataset containing continuous and count outcomes, along with several covariates and a treatment variable, to be used for G-Computation analysis.
}
\usage{
data_gen_gcomp(n_indiv, beta_ttt)
}
\arguments{
\item{n_indiv}{An integer specifying the number of individuals (rows) to generate for the dataset.}
\item{beta_ttt}{A numeric value representing the regression coefficient (effect size) for the treatment variable (\code{ttt}).}
}
\details{
The function generates five predictors and two outcome variables:
\itemize{
  \item \code{x1}: Continuous variable from a Normal distribution (mean=1, sd=0.5).
  \item \code{x2}: Continuous variable from a Normal distribution (mean=2, sd=1.5).
  \item \code{x3}: Binary variable (0/1) from a Binomial distribution (prob=0.4).
  \item \code{x4}: Categorical variable (0/1/2) from a Binomial distribution (size=2, prob=0.8).
  \item \code{ttt}: Binary treatment variable (0/1) from a Binomial distribution (prob=0.5).
}

The linear predictor (\code{lp}) is calculated using fixed coefficients for the covariates (\code{beta_1=4}, \code{beta_2=0.5}, \code{beta_3=1.5}, \code{beta_4=-2.5}) and the user-defined \code{beta_ttt}.

Two outcomes are produced:
\itemize{
  \item \code{y_continuous}: A continuous outcome generated from a Normal distribution with a mean equal to the linear predictor.
  \item \code{y_count}: A count outcome generated from a Poisson distribution with a lambda equal to \code{exp(lp)}.
}
}
\value{
A \code{data.frame} with \code{n_indiv} rows and 7 columns: \code{x1}, \code{x2}, \code{x3}, \code{x4}, \code{ttt}, \code{y_continuous}, and \code{y_count}.
}
\examples{
### Generate a dataset with 1000 individuals and a treatment effect of 2
sim_data <- data_gen_gcomp(n_indiv = 1000, beta_ttt = 2)

### Create formulas
.f <- formula(y_continuous ~ ttt + x1 + x2 + x3 + x4)
.f_count <- formula(y_count ~ ttt + x1 + x2 + x3 + x4)

### Apply G-computation (example for continuous and count outcomes)
gc1 <- gc_continuous(formula=.f, model="all", data=sim_data, group="ttt", 
                boot.type="boot", boot.number=100, seed=8512)
print(gc1)

gc2 <- gc_count(formula=.f_count, model="all", data=sim_data, group="ttt",
                 boot.type="boot", boot.number=10, seed=3210)
print(gc2)

#For continuous the M1-M0 difference corresponds to beta_ttt
#For count the C1/C0 ratio corresponds to exp(beta_ttt)
}

\keyword{datagen}
\keyword{G-computation}