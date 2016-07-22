# computation of multitemporal statistics
# John Truckenbrodt 2016
# this script takes a 3D raster stack and allows for parallel computing of user-defined multitemporal transformations (e.g. smoothing)
# along the  third dimension of the stack

# load required packages
library(abind)
library(raster)
library(parallel)
library(foreach)
library(doParallel)
################################################################

site="Jordan_Azraq"

# define raster stack
raster.name=sprintf("/geonfs01_vol1/ve39vem/swos_process/%s/%s_ASAR_ERS_VV_dB",site,site)

nodata=-99
# nodata=raster.ras[[1]]@file@nodatavalue

# define number of parallel processes
cores=10

# define output file parameters
out.name=sprintf("/geonfs01_vol1/ve39vem/swos_process/%s/%s_ASAR_ERS_VV_dB_filtered",site,site)
out.dtype="FLT4S"

# define multitemporal computation
worker = function(x) {
  require(timeSeries)
  if (length(x[is.finite(x)]) > 3) {
    smooth = timeSeries::smoothSpline(timeSeries::as.timeSeries(x[!is.na(x)]))
    x[!is.na(x)] = smooth$spline
    return(x)
  } else{
    return(x)
  }
}


# no user input beyond this point
##########################

# load data into raster object and 3D array
raster.ras=raster::stack(raster.name)
input=raster::as.array(raster.ras)
input[input==nodata]=NA

# register parallel computing backend
cl=parallel::makeCluster(cores)
doParallel::registerDoParallel(cl,cores)

# define indices for splitting the 3D array into equal parts for each of the parallel processes
strat=split(seq(nrow(input)),cut(seq(nrow(input)),breaks=cores,include.lowest=T))

# function for computing the defined measures
run=function(data,worker){
  return(aperm(apply(data,c(1,2),worker),c(2,3,1)))
}

# function for collecting and combining the subarray computations from the single processes into one final result
collector=function(...){abind::abind(...,along=1)}

# perform the actual computations
out=foreach(i=seq(strat),.combine=collector)%dopar%{
  sub=input[strat[[i]],,]
  run(sub,worker)
}

# unregister the parallel computing backend
unregister(cl)


# transform the computed array into a raster object with geo information
out=if(class(out)=="array") raster::brick(out) else raster::raster(out)
raster::extent(out)=raster::extent(raster.ras)
raster::projection(out)=raster::projection(raster.ras)

# define the NA flag of the output file
naflag=nodata

# write the raster to file
raster::rasterOptions(overwrite=T,datatype=out.dtype,setfileext=F)
raster::writeRaster(out,filename=out.name,format="ENVI",bandorder="BSQ",NAflag=naflag)

# edit the band names of the resulting ENVI file
hdrbands(paste(out.name,".hdr",sep=""),names(raster.ras))
