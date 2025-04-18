% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.polyreg.r
\name{predict.polyreg}
\alias{predict.polyreg}
\title{Predicting cumulative incidence probabilities at a specific time point from direct polynomial regression}
\usage{
\method{predict}{polyreg}(
  formula,
  exposure,
  data,
  coefficient,
  effect.measure1,
  effect.measure2,
  outcome.type,
  inner.optim.method = "optim",
  prob.bound = 1e-05
)
}
\arguments{
\item{formula}{formula Model formula representing outcome, exposure and covariates}

\item{exposure}{character Column name representing the exposure (1 = exposed, 0 = not exposed).}

\item{data}{data.frame Input dataset containing survival data.}

\item{coefficient}{numeric The coefficient of direct polynomial regression.}

\item{effect.measure1}{character Specifies the effect measure for event (RR, OR, SHR).}

\item{effect.measure2}{character Specifies the effect measure for competing risk (RR, OR, SHR).}

\item{outcome.type}{character Specifies the type of outcome (COMPETINGRISK or SURVIVAL).}

\item{inner.optim.method}{character Specifies the method of optimization (nleqslv, optim, multiroot).}

\item{prob.bound}{numeric Small threshold for clamping probabilities. Defaults to 1e-5.}
}
\value{
A list of predicted cumulative incidence probabilities according to exposure and cause from direct polynomial regression.
}
\description{
Predicting cumulative incidence probabilities at a specific time point from direct polynomial regression
}
\examples{
data(bmt)
result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet',
cens.model = Event(time,cause==0)~+1, data = bmt, effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
msummary(result$out_summary, statistic = c("conf.int"), exponentiate = TRUE)
prediction <- predict.polyreg(formula = Event(time, cause)~age+tcell, exposure = 'platelet', data = bmt, coefficient = result$out_coefficient, effect.measure1='RR', effect.measure2='RR', outcome.type='COMPETINGRISK')
}
