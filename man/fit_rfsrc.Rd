% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_surv.R
\name{fit_rfsrc}
\alias{fit_rfsrc}
\title{Wrapper of \code{randomForestSRC::rfsrc}}
\usage{
fit_rfsrc(
  formula,
  data,
  id.var,
  time.var,
  event.var,
  nfold = 1,
  option = list(),
  oob = TRUE,
  tune = TRUE,
  tune.option = list(),
  ...
)
}
\arguments{
\item{formula}{formula used by \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}}

\item{data}{data containing all covariates, follow-up time, event indicator and id}

\item{id.var}{see \code{\link{SDRsurv}}}

\item{time.var}{see \code{\link{SDRsurv}}}

\item{event.var}{see \code{\link{SDRsurv}}}

\item{nfold}{number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting}

\item{option}{a list containing optional arguments passed to \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}. We encourage using a named list. Will be passed to \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} by running a command like \code{do.call(rfsrc, option)}. The user should not specify \code{formula} and \code{data}.}

\item{oob}{whether to use out-of-bag (OOB) fitted values from \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} when sample splitting is not used (\code{nfold=1})}

\item{tune}{whether to tune \code{mtry} and \code{nodesize}.}

\item{tune.option}{a list containing optional arguments passed to \code{\link[randomForestSRC:tune]{randomForestSRC::tune.rfsrc}} if \code{tune=TRUE}; ignored otherwise. \code{doBest} should not be specified.}

\item{...}{ignored}
}
\value{
a \code{\link{pred_surv}} class containing fitted survival curves for individuals in \code{data}
}
\description{
Wrapper of \code{randomForestSRC::rfsrc}
}
