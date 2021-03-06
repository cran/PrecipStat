\name{readGDCN}
\alias{readGDCN}
\alias{readENSEMBLES}
\title{Read data}
\description{
\code{readGDCN} - reads daily data from the GDCN data set
\url{http://www.ncdc.noaa.gov/oa/climate/research/gdcn/gdcn.html}. The
script assumes called filenames with full path, and the function reads
the data in their original format.

\code{readENSEMBLES} - reads daily accumulated precipitation from the
ENSEMBLES RCMs in the netCDF format. 
}
\usage{
readGDCN(name,type="PRCP")
readENSEMBLES(name)
}

\arguments{
\item{name}{Name of the data files containing the GDCN data. Assumes
  full path.}
\item{type}{Type of variable. "PRCP"=precipitation.}
}

\value{\code{readGDCN} returns a vector and \code{readENSEMBLES} a 3D
  matrix, both with attributes holding meta-data. 
}
\author{R.E. Benestad}
\keyword{IO}
\keyword{file}
