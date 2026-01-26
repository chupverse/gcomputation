\name{dataPROPHYVAP}
\alias{dataPROPHYVAP}
\docType{data}
\title{
A Simulated Randomized Clinical Trial from the PROPHYVAP Study
}
\description{
Ventilator-associated pneumonia (VAP) is the first cause of healthcare-associated infections in intensive care units. The PROPHYVAP is a French multicentre, randomized, double-blind, placebo-controlled, clinical trial. The main objective of this study was to determine whether a single dose of ceftriaxone within the 12 hours post-intubation after severe brain injury can decrease the risk of early-onset VAP and mortality at 60 days. All variables in this dataset were simulated according to the real parameters and correlations observed in the original trial. The only exception is \code{HYPOTENSION}, which is a fictional count variable added solely for running exemples in the package.
}
\usage{data(dataPROPHYVAP)}
\format{
  A data frame with 319 observations for the following variables:
  \describe{
    \item{\code{GROUP}}{This character vector represents the randomized treatment (1 for Ceftriaxone and 0 for Placebo)}
    \item{\code{AGE}}{This numeric vector represents the age in years}
    \item{\code{SEX}}{This character vector represents the sex (F=female/M=male)}
    \item{\code{BMI}}{This numeric vector represents the body mass index in kg/m2}  
    \item{\code{DIABETES}}{This character vector represents the diabetes status (Yes/No)}  
     \item{\code{ALCOHOL}}{This character vector represents the alcohol consumption (Yes/No) } 
     \item{\code{SMOKING}}{This character vector represents the tabac status (Yes/No)}  
     \item{\code{INJURY}}{This character vector represents the cause of the injury in 4 classes}  
     \item{\code{GLASGOW}}{This character vector represents the Glasgow scale in 3 classes}  
     \item{\code{PAO2FIO2}}{This character vector represents the PAO2-FIO2 ratio in 2 classes}      
     \item{\code{LEUKO}}{This character vector represents leukocytosis at ICU admission per mm3 in 2 classes}  
     \item{\code{TIME_INTUBATION}}{This numeric vector represents time-to-intubation in hours}  
     \item{\code{TIME_TRT}}{This numeric vector represents time-to-treatment in hours}  
     \item{\code{VAP}}{This character vector represents ventilatory associated pneumonia (1 for event and 0 otherwise)}  
     \item{\code{TIME_DEATH}}{This numeric vector represents time-to-death in days (for a follow-up at 60 days)}  
     \item{\code{DEATH}}{This character vector represents death status at the end of the follow-up (1 for event and 0 otherwise)}
     \item{\code{HYPOTENSION}}{This numeric vector represents the number of episodes of severe hypotension. This is a fictional variable generated via a Poisson distribution specifically for demonstration purposes in this package.}
  }
}

\references{
Dahyot-Fizelier et al. Ceftriaxone to prevent early ventilator-associated pneumonia in patients with acute brain injury: a multicentre, randomised, double-blind, placebo-controlled, assessor-masked superiority trial. Lancet Respir Med, 12:375-385, 2024. <doi:10.1016/S2213-2600(23)00471-X>.
}
\examples{
data(dataPROPHYVAP)

### Kaplan and Meier estimation of the survival at day 60
plot(survfit(Surv(TIME_DEATH, DEATH) ~ GROUP, data = dataPROPHYVAP),
 xlab="Post-randomization time in days", ylab="Patient survival",
 mark.time=FALSE, col=c(2,1), lty=c(2,1))
 
legend("bottomleft", c("Placebo", "Ceftriaxone"), col=1:2, lty=1:2) 
}

\keyword{Dataset}

