# Time Series computation Automation for Raster images
# John Truckenbrodt 2016-2018
# the function tsar of this script takes a 3D raster stack and 
# allows for parallelized computation of user-defined multitemporal statistics 
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
#todo consider a check whether all files can be written
#todo investigate best block size configuration
#todo inform raster package developers of warning given by clusterR in case more than one filename is passed
#todo error handling

#' scalable time-series computations on 3D raster stacks
#' @param raster.name a 3D raster object with dimensions in order lines-samples-time
#' @param workers a named list containing the functions for computation
#' @param cores the number of parallel processes per node
#' @param out.name the name of the output. Either a single file or a folder of separate files 
#' (determined by parameter \code{separate})
#' @param out.bandnames (optional) the names of the output bands; names are determined from the function names 
#' in \code{workers} if left empty
#' @param out.dtype the datatype of the written files.
#' See \code{\link[raster]{dataType}} for options.
#' @param separate should the resulting bands be written to individual files? Otherwise a single ENVI block is written.
#' @param na.in the pixel value for NA in \code{raster.name}
#' @param na.out the pixel value for NA in the output files
#' @param overwrite should the output files be overwritten if they already exist? If \code{separate} all output files are checked
#' @param verbose write detailed information on the progress of function execution?
#' @param nodelist the names of additional server computing nodes accessible via SSH without password
#' @param bandorder the output file pixel arrangement,one of 'BIL', 'BIP' or 'BSQ'
#' @param maxmemory the maximum memory in Mb used per node
#' @param compress_tif should the written GeoTiff files be compressed?
#' @param mask an additional file or raster layer; computations on raster.name are only performed where mask is 1, 
#' otherwise NA is returned for all resulting layers
#' @return None
#' @export
#' @seealso \code{\link[raster]{stack}}, \code{\link[raster]{calc}}, \code{\link[raster]{clusterR}}, \code{\link[snow]{makeCluster}}

tsar=function(raster.name, workers, cores, out.name, out.bandnames=NULL, out.dtype="FLT4S", 
              separate=T, na.in=NA, na.out=-99, overwrite=F, verbose=T, nodelist=NULL, 
              bandorder="BSQ", maxmemory=100, compress_tif=F, mask=NULL){
  require(raster)
  require(snow)
  require(doSNOW)
  
  # abort if multiple values for out.dtype are defined
  if(length(out.dtype)>1){
    stop("multiple values for out.dtype")
  }
  
  if(!separate&&file.exists(out.name)){
    stop("target file/directory already exists")
  }
  
  #the following would be better than the above check since file.exists is also TRUE if its input is in fact a directory.
  #thus, it would not be possible to create a file with a name of an already existing directory.
  #problem is, that the raster package tries to delete any item of that name if overwrite=TRUE.
  #thus, it would try to delete a directory although a file is to be written.
  # if(!separate&&file.exists(out.name)&&!file.info(out.name)$isdir){
  #   stop("target file already exists")
  # }
  
  dir.create(dirname(out.name),showWarnings=F,recursive=T)
  
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
      stop(sprintf("length mismatch of defined band names (%i) and sample output (%i)"
                   ,length(out.bandnames),out.nbands))
    }
  }
  ###################################################
  #if a mask is defined, append it to the stack
  
  if(!is.null(mask)){
    
    if(class(mask)=="character"){
      mask.ras=raster::raster(mask)
    }else if(class(mask)=="RasterLayer"){
      mask.ras=mask
    }else{
      stop("mask input must be a filename or a raster stack")
    }
    ras.in=raster::addLayer(ras.in,mask.ras)
    apply_mask=T
  }else{
    apply_mask=F
  }
  ###################################################
  # create names of the files to be produced and check whether any of them already exist
  
  if(separate){
    # create output names
    outdir=tools::file_path_sans_ext(out.name)
    dir.create(outdir,showWarnings=F,recursive=T)
    
    out.name=unname(sapply(bandnames,function(x)paste0(outdir,"/",x,".tif")))
    
    if(!overwrite){
      # group output names by responsible workers
      indices=unlist(lapply(names(workers),function(x)rep(x,length(sample.out[[x]]))))
      namegroups=split(out.name,indices)
      
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
      
      out.name=out.name[indices]
    }
    
    out.nfiles=length(out.name)
  }
  ###################################################
  # define the function for executing the workers
  
  #further reading/improvement:
  #https://stackoverflow.com/questions/44593123/raster-calculation-on-rasterstack-only-if-not-na-in-other-rasterlayer
  
  run=function(x){
    if(apply_mask){
      if(!is.na(x[length(x)])&&x[length(x)]==1){
        result=unlist(lapply(workers,function(fun)fun(x[-length(x)])))
      }else{
        result=rep(NA,out.nbands)
      }
    }else{
      result=unlist(lapply(workers,function(fun)fun(x)))
    }
    return(result)
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
  # setup options for file writing
  
  if(!separate&&out.nbands>1){
    # case I: multiple bands into a single ENVI block
    
    raster::rasterOptions(setfileext=F,overwrite=T)
    format="ENVI"
    options=""
  }else{
    # case II: a single band written to a single GeoTiff
    # case III: multiple bands each written to a single-band GeoTiff
    
    raster::rasterOptions(setfileext=T,overwrite=T)
    format="GTiff"
    options=if (compress_tif) c("COMPRESS=DEFLATE", "PREDICTOR=2") else c("COMPRESS=NONE")
  }
  ###################################################
  # setup memory and execute the computations
  
  # compute the maximum number of cells which can be held in memory and pass it as raster package option
  cells=maxmemory/8*1024*1024
  raster::rasterOptions(maxmemory=cells, chunksize=cells/100)
  
  #add a progressbar if verbose=TRUE and prevent printing execution time in any case 
  #as this is done by custom function hms_span at the very end
  raster::rasterOptions(progress=if(verbose) "text" else "", timer=F)

  # run the processing
  ras.out=raster::clusterR(ras.in, raster::calc, args=list(fun=run), cl=cl, bylayer=separate,
                           filename=out.name, bandorder=bandorder, NAflag=na.out, 
                           format=format, datatype=out.dtype, options=options)
  
  # edit the band names of the resulting ENVI file to carry information of the computed measures
  if(format=="ENVI")hdrbands(paste0(out.name,".hdr"),bandnames)
  ###################################################
  # unregister parallel computing backend
  snow::stopCluster(cl)
  ###################################################
  # clean up and finish
  
  rm(ras.out)
  gc(verbose=F)
  if(verbose)cat(sprintf("elapsed time: %s\n",hms_span(start.time,Sys.time())))
}

