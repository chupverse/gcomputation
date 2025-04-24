\name{dataPROPHYVAP}
\alias{dataPROPHYVAP}
\docType{data}
\title{
A Simulated Data Set from the PROPHYVAP Study: a Randomized Clinical Trial in Intensive Care Unit
}
\description{
Ventilator-associated pneumonia (VAP) is the first cause of healthcare-associated infections in intensive
care units. The goal of this study is to determine whether a single dose of ceftriaxone within
the 12 hours postintubation after severe brain injury can decrease the risk of early-onset VAP.
The PROPHYVAP is a French multicentre, randomised, double-blind, placebo-controlled,
clinical trial.
}
\usage{data(dataPROPHYVAP)}
\format{
  A data frame with 319 observations for the following variables:
  \describe{
    \item{\code{GROUP}}{This vector of characters represents }
    \item{\code{AGE}}{This numeric vector represents }
    \item{\code{SEXE}}{This numeric vector represents }
    \item{\code{BMI}}{This numeric vector represents }  
    \item{\code{DIABETES}}{This numeric vector represents }  
     \item{\code{ALCOOL}}{This numeric vector represents } 
     \item{\code{SMOKING}}{This numeric vector represents }  
     \item{\code{INJURY_TYPE}}{This numeric vector represents }  
     \item{\code{GLASGOW}}{This numeric vector represents }  
     \item{\code{PAO2FIO2}}{This numeric vector represents }      
     \item{\code{LEUKO}}{This numeric vector represents }  
     \item{\code{DELAI_INTUBATION}}{This numeric vector represents }  
     \item{\code{DELAI_TRT}}{This numeric vector represents }  
     \item{\code{VAP_EARLY}}{This numeric vector represents }  
     \item{\code{FOLLOWUP_J60}}{This numeric vector represents }  
     \item{\code{DEATH_J60}}{This numeric vector represents }
  }
}

\references{
Dahyot-Fizelier et al. Ceftriaxone to prevent early ventilator-associated pneumonia in patients with acute brain injury: a multicentre, randomised, double-blind, placebo-controlled, assessor-masked superiority trial. Lancet Respir Med, 12:375-385, 2024. <doi:10.1016/S2213-2600(23)00471-X>.
}
\examples{
data(dataPROPHYVAP)

### Kaplan and Meier estimation of the recipient survival
plot(survfit(Surv(FOLLOWUP_J60, DEATH_J60) ~ GROUP, data = dataPROPHYVAP),
 xlab="Post-randomiszation time (in days)", ylab="VAP-free survival",
 mark.time=FALSE, col=c(2,1), lty=c(2,1))
 
legend("bottomleft", c("Placebo",
 "Ceftriaxone"), col=1:2, lty=1:2) # VERIF JOE LEGENDE
}

\keyword{Dataset}

