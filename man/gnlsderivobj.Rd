% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/acy.R
\name{gnlsderivobj}
\alias{gnlsderivobj}
\title{GNLS Derivative Objective Function}
\usage{
gnlsderivobj(x, tp, ga, p, la, q)
}
\arguments{
\item{x}{Concentration.}

\item{tp}{Gain-loss Model Top Parameter.}

\item{ga}{Gain AC50.}

\item{p}{Gain power.}

\item{la}{Loss AC50.}

\item{q}{Loss power.}
}
\value{
Value of gnls derivative at x.
}
\description{
Derivative of the gnls function set to zero for top location solver.
}
