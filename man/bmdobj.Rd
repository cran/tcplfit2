% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmdobj.R
\name{bmdobj}
\alias{bmdobj}
\title{BMD Objective Function}
\usage{
bmdobj(bmd, fname, bmr, conc, resp, ps, mll, onesp, partype = 2)
}
\arguments{
\item{bmd}{Benchmark dose.}

\item{fname}{Function name: "exp2", "exp3", "exp4", "exp5", "hillfn", "gnls",
"poly1", "poly2", or "pow".}

\item{bmr}{Benchmark response.}

\item{conc}{Vector of concentrations NOT in log units.}

\item{resp}{Vector of corresponding responses.}

\item{ps}{Named list of parameters.}

\item{mll}{Maximum log-likelihood of winning model.}

\item{onesp}{One-sided p-value.}

\item{partype}{Number for parameter type. Type 1 is y-scaling: a or tp.
Type 2 is x-scaling: b or ga, when available, a otherwise. Type 3 is
power scaling: p when available, then b or ga, then a if no others.
Since bmd is linked to the x-scale, type 2 should always be used. Other
types can also be vulnerable to underflow/overflow.}
}
\value{
Objective function value to find the zero of.
}
\description{
Utility function for bmdbounds
}
