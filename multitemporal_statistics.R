# computation of multitemporal statistics
# John Truckenbrodt 2016
# this script takes a 3D raster stack and allows for parallel computing of user-defined multitemporal statistics 
# (i.e. along the  third dimension of the stack)

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

# define number of parallel processes
cores=5

# define the NA flag of the input and output files
naflag=-99

# define output file parameters

out.name=sprintf("/geonfs01_vol1/ve39vem/swos_process/%s/%s_VV_dB_stats",site,site)
out.dtype="FLT4S"

# define multitemporal statistics
workers=list(
  # minimum=function(x){return(min(x,na.rm=T))},
  # maximum=function(x){return(max(x,na.rm=T))},
  # av=function(x){return(mean(x,na.rm=T))},
  # sdev=function(x){return(sd(x,na.rm=T))},
  p05=function(x){return(quantile(x,probs=.05,na.rm=T,names=F))},
  med=function(x){return(median(x,na.rm=T))},
  p95=function(x){return(quantile(x,probs=.95,na.rm=T,names=F))}
  count=function(x){return(length(x[!is.na(x)]))}
)

# no user input beyond this point
##########################

# load data into raster object and 3D array
raster.ras=raster::stack(raster.name)

for(i in 1:dim(raster.ras)[3]){
  raster.ras[[i]]@file@nodatavalue=naflag
}

input=raster::as.array(raster.ras)

# register parallel computing backend
cl=parallel::makeCluster(cores)
doParallel::registerDoParallel(cl,cores)

# define indices for splitting the 3D array into equal parts for each of the parallel processes
strat=split(seq(nrow(input)),cut(seq(nrow(input)),breaks=cores,include.lowest=T))

# function for computing the defined measures
run=function(data,workers){
  if(length(workers)==1){
    return(apply(data,c(1,2),workers[[1]]))
  }else{
    result=lapply(workers,function(x)apply(data,c(1,2),x))
    return(array(unlist(result),dim=c(nrow(result[[1]]),ncol(result[[1]]),length(result))))
  }
}

# function for collecting and combining the subarray computations from the single processes into one final result
collector=function(...){abind::abind(...,along=1)}

# perform the actual computations
out=foreach(i=seq(strat),.combine=collector)%dopar%{
  sub=input[strat[[i]],,]
  run(sub,workers)
}

# unregister the parallel computing backend
unregister(cl)

# transform the computed array into a raster object with geo information
out=if(class(out)=="array") raster::brick(out) else raster::raster(out)
raster::extent(out)=raster::extent(raster.ras)
raster::projection(out)=raster::projection(raster.ras)

# write the raster to file
raster::rasterOptions(overwrite=T,datatype=out.dtype,setfileext=F)
raster::writeRaster(out,filename=out.name,format="ENVI",bandorder="BSQ",NAflag=naflag)

# edit the band names of the resulting ENVI file to carry information of the computed measures (e.g. minimum, maximum. p05, etc.)
hdrbands(paste(out.name,".hdr",sep=""),names(workers))
