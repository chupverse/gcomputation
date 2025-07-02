\name{dataPROPHYVAP}
\alias{dataPROPHYVAP}
\docType{data}
\title{
A simulated data set from the PROPHYVAP study: a randomized clinical trial in intensive care unit
}
\description{
Ventilator-associated pneumonia (VAP) is the first cause of healthcare-associated infections in intensive
care units. The goal of this study is to determine whether a single dose of ceftriaxone within
the 12 hours postintubation after severe brain injury can decrease the risk of early-onset VAP and mortality at 60 days.
The PROPHYVAP is a French multicentre, randomised, double-blind, placebo-controlled,
clinical trial.
}
\usage{data(dataPROPHYVAP)}
\format{
  A data frame with 319 observations for the following variables:
  \describe{
    \item{\code{GROUP}}{This character vector represents the randomization group (Ceftriaxone/Placebo)}
    \item{\code{AGE}}{This numeric vector represents the age in years}
    \item{\code{SEXE}}{This character vector represents the sexe (F=female/M=male)}
    \item{\code{BMI}}{This numeric vector represents the body mass index}  
    \item{\code{DIABETES}}{This character vector represents the diabetes status (Yes/No)}  
     \item{\code{ALCOOL}}{This character vector represents the alcohol consumption (Yes/No) } 
     \item{\code{SMOKING}}{This character vector represents the tabac (Yes/No)}  
     \item{\code{INJURY_TYPE}}{This character vector represents the cause of the injury (Ischemic stroke/Subarachnoid hemorrhage/Trauma brain injury/Hemorrhagic stroke)}  
     \item{\code{GLASGOW}}{This character vector represents the Glasgow scale in classes (3/4-8/>=9)}  
     \item{\code{PAO2FIO2}}{This character vector represents the PAO2-FIO2 ratio in classes (<100/100-199/>=200)}      
     \item{\code{LEUKO}}{This character vector represents leukocytosis at ICU admission per mm3 in classes (<20000/>=20000)}  
     \item{\code{DELAI_INTUBATION}}{This numeric vector represents time to intubation in hours}  
     \item{\code{DELAI_TRT}}{This numeric vector represents time to treatment in hours}  
     \item{\code{VAP_EARLY}}{This character vector represents ventilatory associated pneumonia (Yes/No)}  
     \item{\code{FOLLOWUP_J60}}{This numeric vector represents follow-up in days}  
     \item{\code{DEATH_J60}}{This character vector represents status at the end of the follow-up (Yes/No)}
  }
}

\references{
Dahyot-Fizelier et al. Ceftriaxone to prevent early ventilator-associated pneumonia in patients with acute brain injury: a multicentre, randomised, double-blind, placebo-controlled, assessor-masked superiority trial. Lancet Respir Med, 12:375-385, 2024. <doi:10.1016/S2213-2600(23)00471-X>.
}
\examples{
data(dataPROPHYVAP)
dataPROPHYVAP$DEATH_J60_num <- ifelse(dataPROPHYVAP$DEATH_J60 == "Yes",1,0)

### Kaplan and Meier estimation of the survival at day 60
plot(survfit(Surv(FOLLOWUP_J60, DEATH_J60_num) ~ GROUP, data = dataPROPHYVAP),
 xlab="Post-randomiszation time (in days)", ylab="Patient survival",
 mark.time=FALSE, col=c(2,1), lty=c(2,1))
 
legend("bottomleft", c("Placebo",
 "Ceftriaxone"), col=1:2, lty=1:2) 
}

\keyword{Dataset}

