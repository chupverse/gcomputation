\name{dataCOHORT}
\alias{dataCOHORT}
\docType{data}
\title{
Simulated Real Word Data Similar to the PROPHYVAP Study
}
\description{
This dataset \code{dataCOHORT} is a simulated observational cohort designed to reflect a real-world data similar to the PROPHYVAP clinical trial (see \code{\link{dataPROPHYVAP}}).
}
\usage{data(dataCOHORT)}
\format{
A data frame with 2000 observations for the following variables:
\describe{
\item{\code{GROUP}}{This character vector represents the observed treatment group (1 for Ceftriaxone and 0 otherwise)}
\item{\code{AGE}}{This numeric vector represents the patient age in years}
\item{\code{SEX}}{This character vector represents the gender (F=female/M=male)}
\item{\code{BMI}}{This numeric vector represents the body mass index in kg/m2}
\item{\code{DIABETES}}{This character vector represents the diabetes status (Yes/No)}
\item{\code{ALCOHOL}}{This character vector represents the alcohol consumption (Yes/No) }
\item{\code{SMOKING}}{This character vector represents the tabaco status (Yes/No)}
\item{\code{INJURY}}{This character vector represents the cause of the injury in 4 classes}
\item{\code{GLASGOW}}{This character vector represents the Glasgow scale in 3 classes }
\item{\code{PAO2FIO2}}{This character vector represents the PAO2-FIO2 ratio in 2 classes}
\item{\code{LEUKO}}{This character vector represents leukocytosis at admission per mm3 in 2 classes}
\item{\code{TIME_INTUBATION}}{This numeric vector represents the time to intubation in hours}
\item{\code{VAP}}{This character vector represents ventilatory associated pneumonia (1 for event and 0 otherwise)} 
\item{\code{TIME_DEATH}}{This numeric vector represents the time to death in days (maximum follow-up at 60 days)}
\item{\code{DEATH}}{This character vector represents death status at the end of the follow-up (1 for event and 0 otherwise)} 
}
}

\references{
Dahyot-Fizelier et al. Ceftriaxone to prevent early ventilator-associated pneumonia in patients with acute brain injury: a multicentre, randomised, double-blind, placebo-controlled, assessor-masked superiority trial. Lancet Respir Med, 12:375-385, 2024. <doi:10.1016/S2213-2600(23)00471-X>.
}

\examples{
data(dataCOHORT)

### Kaplan and Meier estimation of the survival at day 60
plot(survfit(Surv(TIME_DEATH, DEATH) ~ 1, data = dataCOHORT),
 xlab="Time in days", ylab="Patient survival", conf.int=TRUE,
 mark.time=FALSE)
 
}
\keyword{Dataset}