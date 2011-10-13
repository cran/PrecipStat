# R.E. Benestad (rasmus.benestad@physics.org)
# R-script for reading CDCN daily precipitation data
# NCAR, Mesa Lab, Boulder March 31, 2011

# Look for similarity between exponential distribution and 24-hr precip world wide
# Look for systematic relationship between mu, meanP, and meanT, hence a dependence of shape of p.d.f. on
# local mean T(2m) and precip.




readGDCN <- function(name,type="PRCP") {
  # Reads the ASCII files with daily GDCN data:
  cnames <- c(paste("day.",1:31,sep=""),paste("flags.day.",1:31,sep=""))
  reshufle <- rep(0,62); ii <- 0
  for (i in seq(1,62,by=2)) {
    ii <- ii + 1
    reshufle[i:(i+1)] <- seq(ii,62,by=31)
  }
  cnames <- cnames[reshufle]
  x <- read.fwf(name,widths=c(3,8,4,2,4,rep(c(5,2),31)),
                col.names=c("country","stnr","year","month","type",cnames))
  ipick <- is.element(x$type,type)
  if (sum(ipick)>0) {
    dat <- as.matrix(x[ipick,seq(6,67,by=2)])*0.1
    dat[dat < -99] <- NA
    attr(dat,"Data source") <- "GDCN" 
    attr(dat,"year") <- x$year[ipick]
    attr(dat,"month") <- x$month[ipick]
    attr(dat,"URL") <- "http://www.ncdc.noaa.gov/oa/climate/research/gdcn/gdcn.html"
    attr(dat,"Station_number") <- x$stnr[1]
    attr(dat,"Country_code") <- x$country[1]
    attr(dat,"original data name") <- name
    attr(dat,"flags") <- x[ipick,seq(6,66,by=2)]
    attr(dat,"history") <- "read with readGDCN - R-script."
    attr(dat,"Observations") <- rownames(table(x$type))
  } else dat <- NULL
  invisible(dat)
}





readENSEMBLES <- function(name) {
  # Reads the netCDF files with daily RCM data from the ENSEMBLES project:
    print("readENSEMBLES")
    ncid <- open.ncdf(name)
    # Many of the ENSEMBLES files are organised differently with
    # different names for lon/lat dimensions etc.
    nv <- ncid$nvars;  cdfvars <- rep("-",nv)
    nd <- ncid$ndims;  cdfdims <- rep("-",nd)
    for (i in 1:nv) cdfvars[i] <- ncid$var[[i]]$name
    for (i in 1:nd) cdfdims[i] <- ncid$dim[[i]]$name
    print(lower.case(cdfvars))
    print(lower.case(cdfdims))
    ilon <- grep("lon",lower.case(cdfvars))
    if (length(ilon)>0) lons <- get.var.ncdf(ncid,cdfvars[ilon]) else {
      ilon <- grep("lon",lower.case(cdfdims))
      lons <- get.var.ncdf(ncid,cdfdims[ilon])
    }
    ilat <- grep("lat",lower.case(cdfvars))
    if (length(ilat)>0) lats <- get.var.ncdf(ncid,cdfvars[ilat]) else {
      ilat <- grep("lat",lower.case(cdfdims))
      lats <- get.var.ncdf(ncid,cdfdims[ilat])
    }
    time <- get.var.ncdf(ncid,"time"); nt <- length(time)
    if (length(dim(lons))==2) nx <- dim(lons)[1] else nx <- length(lons)
    if (length(dim(lats))==2) ny <- dim(lats)[2] else ny <- length(lats)
    print(paste("Dimensions=",nx,"x",ny,"x",nt))
    
    #The total interval is 1950:2100
    rcm <- get.var.ncdf(ncid,"pr",
                        start=c(round(nx/4),round(ny/4),1),
                        count=c(round(3*nx/4),round(3*ny/4),nt))
    close.ncdf(ncid)

    if (length(dim(lons))==2) {
      lons <- lons[round(nx/4):round(3*nx/4),round(ny/4):round(3*ny/4)]
      lats <- lats[round(nx/4):round(3*nx/4),round(ny/4):round(3*ny/4)]
    } else {
      lons <- lons[round(nx/4):round(3*nx/4)]
      lats <- lats[round(ny/4):round(3*ny/4)]
    }
    date <- caldat(julday(01,01,1950) + time/24)
    attr(rcm,"date") <- date
    attr(rcm,"longitude") <- lons
    attr(rcm,"latitude") <- lats
    invisible(rcm)
}

