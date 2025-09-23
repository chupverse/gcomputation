\name{dataCOHORT}
\alias{dataCOHORT}
\docType{data}
\title{
A simulated data set from the PROPHYVAP study context: an observational cohort study
}
\description{
This dataset \code{dataCOHORT} is a simulated observational cohort study, with a total of 2000 patients, designed to reflect a real-world scenario similar to the PROPHYVAP clinical trial.The goal of this study is to determine whether a single dose of ceftriaxone within
the 12 hours postintubation after severe brain injury can decrease the risk of early-onset VAP and mortality at 60 days. Unlike the randomized trial \code{\link{dataPROPHYVAP}}, the treatment group (\code{GROUP}) in this simulated cohort is influenced by patient characteristics, reflecting confounding factors.
}
\usage{data(dataCOHORT)}
\format{
A data frame with 2000 observations for the following variables:
\describe{
\item{\code{GROUP}}{This character vector represents the treatment group (Ceftriaxone/Untreated).}
\item{\code{AGE}}{This numeric vector represents the age in years}
\item{\code{SEX}}{This character vector represents the sex (F=female/M=male)}
\item{\code{BMI}}{This numeric vector represents the body mass index}
\item{\code{DIABETES}}{This character vector represents the diabetes status (Yes/No)}
\item{\code{ALCOOL}}{This character vector represents the alcohol consumption (Yes/No) }
\item{\code{SMOKING}}{This character vector represents the tabac (Yes/No)}
\item{\code{INJURY}}{This character vector represents the cause of the injury (Ischemic stroke/Subarachnoid hemorrhage/Trauma brain injury/Hemorrhagic stroke)}
\item{\code{GLASGOW}}{This character vector represents the Glasgow scale in classes (3/4-8/>=9)}
\item{\code{PAO2FIO2}}{This character vector represents the PAO2-FIO2 ratio in classes (<100/100-199/>=200)}
\item{\code{LEUKO}}{This character vector represents leukocytosis at ICU admission per mm3 in classes (<20000/>=20000)}
\item{\code{TIME_INTUBATION}}{This numeric vector represents time to intubation in hours}
\item{\code{VAP}}{This character vector represents ventilatory associated pneumonia (Yes/No)}
\item{\code{TIME_DEATH}}{This numeric vector represents time to death in days (for a follow-up at 60 days)}
\item{\code{DEATH}}{This character vector represents status at the end of the follow-up (Yes/No)}
}
}
\references{
Dahyot-Fizelier et al. Ceftriaxone to prevent early ventilator-associated pneumonia in patients with acute brain injury: a multicentre, randomised, double-blind, placebo-controlled, assessor-masked superiority trial. Lancet Respir Med, 12:375-385, 2024. <doi:10.1016/S2213-2600(23)00471-X>.

(Note: This \code{dataCOHORT} dataset is a simulation for methodological exploration, not the original trial data. The original PROPHYVAP study [2] was a randomized, double-blind, placebo-controlled clinical trial.)
}
\examples{
data(dataCOHORT)

dataCOHORT$DEATH_num <- ifelse(dataCOHORT$DEATH == "Yes",1,0)

### Kaplan and Meier estimation of the survival at day 60
plot(survfit(Surv(TIME_DEATH, DEATH_num) ~ GROUP, data = dataCOHORT),
 xlab="Post-randomiszation time (in days)", ylab="Patient survival",
 mark.time=FALSE, col=c(2,1), lty=c(2,1))
 
legend("bottomleft", c("Untreated",
 "Ceftriaxone"), col=1:2, lty=1:2) 
}
\keyword{Dataset}