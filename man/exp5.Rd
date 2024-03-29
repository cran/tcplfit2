% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcplObj.R
\name{exp5}
\alias{exp5}
\title{Exponential 5 Model}
\usage{
exp5(ps, x)
}
\arguments{
\item{ps}{Vector of parameters: tp,ga,p,er}

\item{x}{Vector of concentrations (regular units)}
}
\value{
Vector of model responses
}
\description{
\eqn{f(x) = tp*(1-2^{(-(x/ga)^p)})}
}
\examples{
exp5(c(1,2,3),1)

}
