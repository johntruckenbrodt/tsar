# Time Series computation Automation for Raster images
# John Truckenbrodt 2016-2017
# the function tsar of this script takes a 3D raster stack and 
# allows for parallelized computation of user-defined multitemporal statistics 
################################################################
# function to unregister parallel computing backend
unregister=function(cluster){
  require(parallel)
  require(foreach)
  parallel::stopCluster(cluster)
  env=foreach:::.foreachGlobals
  rm(list=ls(name=env),pos=env)
}
################################################################
# edit the band names of an ENVI hdr file
hdrbands=function(hdrfile, names){
  hdr=readLines(hdrfile)
  index=grep("band names",hdr)
  while(!grepl("}",hdr[index]))index=index+1
  bandnames=paste("band names = {",paste(names,collapse=",\n"),"}",sep="")
  hdr=c(hdr[-(grep("band names",hdr):index)],bandnames)
  writeLines(hdr,hdrfile)
}
################################################################
hms_span=function(start,end){
  #https://stackoverflow.com/questions/32100133/print-the-time-a-script-has-been-running-in-r
  dsec=as.numeric(difftime(end,start,unit="secs"))
  hours=floor(dsec/3600)
  minutes=floor((dsec-3600*hours)/60)
  seconds=dsec-3600*hours-60*minutes
  paste0(sapply(c(hours,minutes,seconds),function(x){formatC(x,width=2,format="d",flag="0")}),collapse=":")
}
################################################################
# main function for computations
#todo change na.in=NA to something else
#todo check Windows compatibility of environment passing (snow::clusterExport)
#todo input variable checks!
#todo compute .maxcombine from memory used
#todo option for defining mask file
#todo make every node combine by itself
#todo enable multiple na.out values
#todo memory monitoring
#todo test function raster::calc (in combination with directly writing with option bylayer=TRUE)

#' scalable time-series computations on 3D raster stacks
#' @param raster.name a 3D raster object with dimensions in order lines-samples-time
#' @param workers a named list containing the functions for computation
#' @param cores the number of parallel processes per node
#' @param out.name the name of the output. Either a single file or a folder of separate files (determined by parameter \code{separate})
#' @param out.bandnames (optional) the names of the output bands; names are determined from the function names in \code{workers} if left empty
#' @param out.dtype the datatype of the written files. This can either be a single value or a vector of values of same length as the files to be written.
#' See \code{\link[raster]{dataType}} for possible values.
#' @param separate should the resulting band be written to individual files? Otherwise a single ENVI block is written.
#' @param na.in the pixel value for NA in \code{raster.name}
#' @param na.out the pixel value for NA in the output files
#' @param overwrite should the output files be overwritten if they already exist? If \code{separate} all output files are checked
#' @param verbose write detailed information on the progress of function execution?
#' @param nodelist the names of additional server computing nodes accessible via SSH without password
#' @return None
#' @export
#' @seealso \code{\link[raster]{stack}}, \code{\link[foreach]{foreach}}, \code{\link[snow]{makeCluster}}

tsar=function(raster.name, workers, cores, out.name, out.bandnames=NULL, out.dtype="FLT4S", 
              separate=T, na.in=NA, na.out=-99, overwrite=T, verbose=T, nodelist=NULL){
  require(abind)
  require(raster)
  require(foreach)
  #require(parallel)
  #require(doParallel)
  require(doSNOW)
  
  # abort if files are to be written in a single file but multiple values for out.dtype are defined
  if(!separate&&length(out.dtype)>1){
    stop("multiple values for out.dtype, but only one file to be written")
  }
  
  if(!separate&&file.exists(out.name)){
    stop("target file already exists")
  }
  
  start.time=Sys.time()
  
  ###################################################
  # load data into raster object
  
  if(class(raster.name)=="character"){
    ras.in=raster::stack(raster.name)
  }else if(class(raster.name)=="RasterStack"){
    ras.in=raster.name
  }else{
    stop("input must be a filename, a vector of filenames or a raster stack")
  }
  # set the nodata flag
  if(!is.na(na.in))raster::NAvalue(ras.in)=na.in
  
  # define the input dimensions
  rows=dim(ras.in)[1]
  cols=dim(ras.in)[2]
  bands=dim(ras.in)[3]
  ###################################################
  # perform a test computation on a single pixel, assess the output length and parse the bandnames accordingly
  
  #assess the length of the worker output
  sample=ras.in[rows%/%2,cols%/%2,][1,]
  sample.out=lapply(workers,function(x)x(sample))
  out.nbands=length(unlist(sample.out))
  
  # define the bandnames based on the user input (if existent) or the worker names
  if(is.null(out.bandnames)){
    bandnames=c()
    for(name in names(workers)){
      len=length(sample.out[[name]])
      names=if(len==1) name else sapply(seq(len),function(x)paste(name,x,sep="_"))
      bandnames=c(bandnames,names)
    }
  }else{
    if(length(out.bandnames)==out.nbands){
      bandnames=out.bandnames
    }else{
      stop(sprintf("length mismatch of defined band names (%i) and sample output (%i)",length(out.bandnames),out.nbands))
    }
  }
  ###################################################
  # create names of the files to be produced and check whether any of them already exist
  
  if(separate){
    # create output names
    outdir=tools::file_path_sans_ext(out.name)
    dir.create(outdir,showWarnings=F)
    
    outnames=unname(sapply(bandnames,function(x)paste0(outdir,"/",x,".tif")))
    
    if(!overwrite){
      # group output names by responsible workers
      indices=unlist(lapply(names(workers),function(x)rep(x,length(sample.out[[x]]))))
      namegroups=split(outnames,indices)
      
      # check existence of output files and evaluate whether a worker has to be executed
      jobs=sapply(names(workers),function(x)!all(sapply(unlist(namegroups[[x]]),function(y)file.exists(y))))
      
      # reduce the list of workers to those that need to be executed
      workers=workers[jobs]
      
      # abort if no worker is left to be executed
      if(length(workers)==0){
        message("no file to be written")
        return()
      }
      
      # reduce the list of outnames and data types to those bands produced by the updated list of workers
      
      indices=which(indices %in% names(workers))
      
      outnames=outnames[indices]
      
      if(length(out.dtype)==out.nbands){
        out.dtype=out.dtype[indices]
      }
    }
    
    out.nfiles=length(outnames)
    
    # evaluate the data type(s) defined by out.dtype
    if(length(out.dtype)==1&&out.nfiles>1){
      out.dtype=rep(out.dtype,out.nfiles)
    }
    if(length(out.dtype)!=out.nfiles){
      stop(sprintf("length mismatch of defined data types in out.dtype (%i) and files to be written (%i)",length(out.dtype),out.nbands))
    }
  }
  ###################################################
  # define the functions for executing the workers
  
  run=function(x){
    return(unlist(lapply(workers,function(fun)fun(x))))
  }
  ###################################################
  # set up the environment to be passed to the parallel workers
  
  # this is likely not to work under Windows, see this link:
  # http://stackoverflow.com/questions/17345271/r-how-does-a-foreach-loop-find-a-function-that-should-be-invoked
  e=globalenv()
  functions=ls(e)[sapply(ls(e),function(x)class(e[[x]]))=="function"]
  new.env=new.env()
  for(fun in functions)assign(fun,value=e[[fun]],envir=new.env)
  
  # list all currently loaded packages (to be passed to the parallel workers)
  packages=gsub("package:","",grep("package",search(),value=T))

  ###################################################
  # register parallel computing backend and pass the newly created environment
  
  if(!is.null(nodelist)){
    # http://stackoverflow.com/questions/25370603/turn-on-all-cpus-for-all-nodes-on-a-cluster-snow-snowfall-package
    
    hosts=unlist(lapply(nodelist,function(x)rep(x,cores)))
    if(length(hosts)>128){
      message=sprintf("A maximum of 128 connections may be opened at once. %i nodes x %i cores = %i connections.",
                      length(nodelist),cores,length(nodelist)*cores)
      stop(message)
    }
    hostList=lapply(hosts,function(x)list(host=x))
    cl=snow::makeCluster(hostList,type="SOCK")
  }else{
    cl=snow::makeCluster(cores,type="SOCK")
  }
  doSNOW::registerDoSNOW(cl)
  snow::clusterExport(cl,list=ls(new.env),envir=new.env)
  ###################################################
  # execute the computations
  
  ras.out=raster::clusterR(ras.in, raster::calc, args=list(fun=run),cl=cl)
  ###################################################
  # unregister parallel computing backend
  snow::stopCluster(cl)
  ###################################################
  # write the results to disk
  
  if(verbose)cat("writing results to disk\n")
  
  # case I: multiple bands into a single ENVI block
  if(!separate&&out.nbands>1){
    
    raster::rasterOptions(overwrite=T,datatype=out.dtype,setfileext=F)
    raster::writeRaster(ras.out,filename=out.name,format="ENVI",bandorder="BSQ",NAflag=na.out)

    # edit the band names of the resulting ENVI file to carry information of the computed measures
    # (i.e. the names of the workers, e.g. minimum, maximum, p05, etc.)
    hdrbands(paste(out.name,".hdr",sep=""),bandnames)
    
  }else if(out.nbands==1){
    # case II: a single band written to a single GeoTiff
    
    raster::rasterOptions(overwrite=T,datatype=out.dtype,setfileext=T)
    raster::writeRaster(ras.out,filename=paste0(out.name,".tif"),format="GTiff",bandorder="BSQ",NAflag=na.out,options=c("COMPRESS=NONE"))
    
  }else{
    # case III: multiple bands each written to a single-band GeoTiff
    
    raster::rasterOptions(overwrite=T,setfileext=T)
    for(i in seq(out.nbands)){
      if(verbose)cat(sprintf("..%s\n",outnames[i]))
      raster::writeRaster(ras.out[[i]],filename=outnames[i],format="GTiff",bandorder="BSQ",
                          NAflag=na.out,options=c("COMPRESS=NONE"),datatype=out.dtype[i])
    }
  }
  ###################################################
  # clean up and finish
  
  rm(ras.out)
  gc(verbose=F)
  if(verbose)cat(sprintf("elapsed time: %s\n",hms_span(start.time,Sys.time())))
}

