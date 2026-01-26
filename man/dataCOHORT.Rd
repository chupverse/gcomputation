\name{dataCOHORT}
\alias{dataCOHORT}
\docType{data}
\title{
A Simulated Cohort from the PROPHYVAP Study
}
\description{
This dataset \code{dataCOHORT} is a simulated observational cohort study designed to reflect a real-world scenario similar to the PROPHYVAP clinical trial. The main objective of this study was to determine whether a single dose of ceftriaxone within the 12 hours post-intubation after severe brain injury can decrease the risk of early-onset VAP and mortality at 60 days. Unlike the randomized trial \code{\link{dataPROPHYVAP}}, the treatment allocation was influenced by patient characteristics.
}
\usage{data(dataCOHORT)}
\format{
A data frame with 2000 observations for the following variables:
\describe{
\item{\code{GROUP}}{This character vector represents the treatment group (1 for Ceftriaxone and 0 for Placebo)}
\item{\code{AGE}}{This numeric vector represents the age in years}
\item{\code{SEX}}{This character vector represents the sex (F=female/M=male)}
\item{\code{BMI}}{This numeric vector represents the body mass index in kg/m2}
\item{\code{DIABETES}}{This character vector represents the diabetes status (Yes/No)}
\item{\code{ALCOHOL}}{This character vector represents the alcohol consumption (Yes/No) }
\item{\code{SMOKING}}{This character vector represents the tabac status (Yes/No)}
\item{\code{INJURY}}{This character vector represents the cause of the injury in 4 classes}
\item{\code{GLASGOW}}{This character vector represents the Glasgow scale in 3 classes }
\item{\code{PAO2FIO2}}{This character vector represents the PAO2-FIO2 ratio in 2 classes}
\item{\code{LEUKO}}{This character vector represents leukocytosis at ICU admission per mm3 in 2 classes}
\item{\code{TIME_INTUBATION}}{This numeric vector represents time-to-intubation in hours}
\item{\code{VAP}}{This character vector represents ventilatory associated pneumonia (1 for event and 0 otherwise)} 
\item{\code{TIME_DEATH}}{This numeric vector represents time-to-death in days (maximum follow-up at 60 days)}
\item{\code{DEATH}}{This character vector represents death status at the end of the follow-up (1 for event and 0 otherwise)} 
}
}

\references{
Dahyot-Fizelier et al. Ceftriaxone to prevent early ventilator-associated pneumonia in patients with acute brain injury: a multicentre, randomised, double-blind, placebo-controlled, assessor-masked superiority trial. Lancet Respir Med, 12:375-385, 2024. <doi:10.1016/S2213-2600(23)00471-X>.
}

\examples{
data(dataCOHORT)

dataCOHORT$DEATH_num <- ifelse(dataCOHORT$DEATH == "Yes",1,0)

### Kaplan and Meier estimation of the survival at day 60
plot(survfit(Surv(TIME_DEATH, DEATH_num) ~ GROUP, data = dataCOHORT),
 xlab="Time in days", ylab="Patient survival",
 mark.time=FALSE, col=c(2,1), lty=c(2,1))
 
legend("bottomleft", c("Untreated","Ceftriaxone"), col=1:2, lty=1:2)
}
\keyword{Dataset}