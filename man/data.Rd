\name{mu.qp}
\alias{mu.qp}
%\alias{mu.qp.jjas}
%\alias{mu.qp.rcm}
%\alias{mu.qp.escn}
%\alias{mu.qp.escn.jjas}
%\alias{mu.qp.res}
%\alias{mu.qp.lowq}
\alias{mu.qp.world}
\alias{mu.qp.world.jjas}
%\alias{mu.qp.met.no}
%\alias{mu.qp.met.no.jjas}
\alias{gdcn.inv}
\alias{f.wet.1990.1999}
\alias{mu.1990.1999}
\alias{f.wet.2080.2099}
\alias{mu.2080.2099}
\title{Pre-calculated results and GDCN meta-data}
\description{
\code{mu.qp} - contains the main results from \code{qplotGDCN}, and is a
list object containing matrices x1 and x2. The dimensions of the \code{x1} and
\code{x2} are \code{[N.p,N.l]}, where \code{N.p} is the number of
quantiles for each record (default: 19) and \code{N.l} the number of
locations (default: 13549). \code{mu.qp.jjas} is the same, but for the
"summer" months: June-September.

\code{mu.qp.escn} - the same as for \code{mu.qp}, but holding the
results for the ECA\&D data. Also \code{mu.qp.escn.jjas}.

\code{mu.qp.rcm}  - the same as for \code{mu.qp}, but holding the
results for the NARCCAP and ENSEMBLES RCM results. The different
projects are tagged through \code{as.character(attr(mu.qp.rcm$x1,"experiment"))}.

\code{mu.qp.res} - similar to \code{mu.qp.rcm}, but for the ENSEMBLES
HadRM3 for high and low spatial resolution.

\code{mu.qp.lowq} - the same as for \code{mu.qp}, but for
\code{prs=seq(0.05,0.95,by=0.05)}

\code{mu.qp.world} holds the extracted percentiles and mean from the
combined GDCN disks 1 \& 2, covering large portion of the world.

\code{f.wet.1990.1999} and \code{mu.1990.1999} hold the wet-day
frequency and the wet-day mean from the NCAR N_CCSM3.0 (CMIP3
20c3m.run1) for the period 1990--1999. \code{f.wet.2080.2099} and
\code{mu.2080.2099} hold corresponding values, but for the SRES A1b,
run2, and the period 2080--2099. These GCM results are used here merely
as a simple and crude demonstration for how \code{qqplotter} can be used
in downscaling (\code{qqscenario} in the 'qqplotter' package), but proper downscaling ought
to use results from  empirical-statistical downscaling or adjusted
results from regional climate models. Different runs were deliberately
chosen for the two time slices to avoid that these results are taken as
projections for the future. The intention here is only to provide these
data as a demonstration of the method. 
}
\usage{
data(mu.qp)
%data(mu.qp.jjas)
%data(mu.qp.lowq)
data(mu.qp.world)
data(mu.qp.world.jjas)
%data(mu.qp.met.no)
%data(mu.qp.met.no.jjas)
%data(mu.qp.rcm)
%data(mu.qp.escn)
%data(mu.qp.escn.jjas)
data(gdcn.inv)
data(f.wet.1990.1999)
data(mu.1990.1999)
data(f.wet.2080.2099)
data(mu.2080.2099)
}


\value{Matrices with  dimensions  \code{[N.p,N.l]}, where \code{N.p} is
  the number of quantiles for each record and \code{N.l} the number of
  locations. The matrices have a number of attributes with meta-data
}
\author{R.E. Benestad}
\examples{
data(mu.qp)
names(attributes(mu.qp))
}
\keyword{datasets}
