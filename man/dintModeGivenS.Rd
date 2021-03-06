% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dintModeGivenS.R
\name{dintModeGivenS}
\alias{dintModeGivenS}
\title{Mode probabilities}
\usage{
dintModeGivenS(S, mod, pars)
}
\arguments{
\item{S}{number of spcies}

\item{mod}{function specifying the SAD density function (e.g. \code{dfish})}

\item{pars}{vector of parameters for the SAD model}
}
\value{
the probability that the first bin is not the modal bin
}
\description{
Probability of a mode across Preston octaves
}
\details{
For Preston octaves 0 through 12, this function calculates the probability that each
octave could be a modal octave (mode defined as equal to the maximum frequency of species).
More than one octave could meet this criterian, so the returned vector of probabilities will
not neccesarily sum to 1.  The probability of a mode at octave i is
\eqn{\sum_s P(o_i = s) P(\text{all } o_{j \neq i} \leq s)} where \eqn{P(o_i = s)} is simply
\code{dbinom(s, S, p[i])} where \code{S} is the total number of species and \code{p[i]} is the
summed probability of the SAD corresponding to octave i. And where
\eqn{P(\text{all } o_{j \neq i} \leq s)} is given by \code{pmultinom(rep(s, 12), S - s, p[-i])}
}
\author{
Andy Rominger <ajrominger@gmail.com>
}
