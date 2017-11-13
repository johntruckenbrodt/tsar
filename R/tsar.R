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
tsar=function(raster.name, workers, cores, out.name, out.bandnames=NULL, out.dtype="FLT4S", 
              separate=T, na.in=NA, na.out=-99, overwrite=T, verbose=T, nodelist=NULL){
  require(abind)
  require(raster)
  require(foreach)
  #require(parallel)
  #require(doParallel)
  require(doSNOW)
  
  start.time=Sys.time()
  
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
  
  #assess the length of the worker output
  sample=ras.in[rows%/%2,cols%/%2,][1,]
  sample.out=lapply(workers,function(x)x(sample))
  out.bands=length(unlist(sample.out))
  
  # define the bandnames based on the user input (if existent) or the worker names
  if(is.null(out.bandnames)){
    bandnames=c()
    for(name in names(workers)){
      len=length(sample.out[[name]])
      names=if(len==1) name else sapply(seq(len),function(x)paste(name,x,sep="_"))
      bandnames=c(bandnames,names)
    }
  }else{
    if(length(out.bandnames)==out.bands){
      bandnames=out.bandnames
    }else{
      stop(sprintf("length mismatch of defined band names (%i) and sample output (%i)",length(out.bandnames),out.bands))
    }
  }
  
  # check whether any of the files which will be produced already exist
  if(overwrite==F){
    if(separate){
      #create output names and group them by responsible workers
      outnames=unname(sapply(bandnames,function(x)paste0(tools::file_path_sans_ext(out.name),"_",x,".tif")))
      indices=lapply(names(workers),function(x)rep(x,length(sample.out[[x]])))
      namegroups=split(outnames,unlist(indices))
      
      #check existence of output files and evaluate whether a worker has to be started
      jobs=sapply(names(workers),function(x)!all(sapply(unlist(namegroups[[x]]),function(y)file.exists(y))))
      
      #reduce the list of workers to those that need to be executed
      workers=workers[jobs]
    }
  }
  
  # abort if no worker is left to executed or if results are to be written to one file which already exists
  if(length(workers)==0|(!separate&file.exists(out.name))){
    message("no file to be written")
    return()
  }
  
  # set up the environment to be passed to the parallel workers
  # this is likely not to work under Windows, see this link:
  # http://stackoverflow.com/questions/17345271/r-how-does-a-foreach-loop-find-a-function-that-should-be-invoked
  e=globalenv()
  functions=ls(e)[sapply(ls(e),function(x)class(e[[x]]))=="function"]
  new.env=new.env()
  for(fun in functions)assign(fun,value=e[[fun]],envir=new.env)
  
  # list all currently loaded packages (to be passed to the parallel workers)
  packages=gsub("package:","",grep("package",search(),value=T))
  
  # function for computing the defined measures
  run=function(data,workers){
    result=apply(data,c(1,2),function(x){return(unlist(lapply(workers,function(fun)fun(x))))})
    return(if(out.bands>1) aperm(result,c(2,3,1)) else result)
  }
  
  # function for collecting and combining the subarray computations from the single parallel processes into one final result
  # collector=function(...){abind::abind(...,along=1)}
  collector=list
  
  # register parallel computing backend and pass the newly created function environment
  #cl=parallel::makeCluster(cores)
  #doParallel::registerDoParallel(cl,cores)
  
  if(!is.null(nodelist)){
    # http://stackoverflow.com/questions/25370603/turn-on-all-cpus-for-all-nodes-on-a-cluster-snow-snowfall-package
    
    hosts=unlist(lapply(nodelist,function(x)rep(x,cores)))
    if(length(hosts)>128){
      message=sprintf("A maximum of 128 connections may be opened at once. %i nodes x %i cores = %i connections.",
                      length(nodelist),cores,length(nodelist)*cores)
      stop(message)
    }
    hostList=lapply(hosts,function(x)list(host=x))
    cl=makeCluster(hostList,type="SOCK")
    processes=length(hosts)
  }else{
    cl=snow::makeCluster(cores)
    processes=cores
  }
  
  doSNOW::registerDoSNOW(cl)
  snow::clusterExport(cl,list=ls(new.env),envir=new.env)
  
  # define indices for splitting the 3D array along the y-axis into equal parts for each of the parallel processes
  strat=split(seq(nrow(ras.in)),cut(seq(nrow(ras.in)),breaks=processes,include.lowest=T))
  
  # perform the actual computations
  # extra option for exporting the whole global enivronment: .export=ls(envir=globalenv())
  out.arr=foreach(i=seq(strat),.combine=collector,.verbose=verbose,.packages=packages,
                  .multicombine=T,.maxcombine=length(strat))%dopar%{
    indices=strat[[i]]
    rows.sub=length(indices)
    
    arr=raster::getValuesBlock(ras.in,row=min(indices),nrows=rows.sub)
    dim(arr)=c(cols,rows.sub,bands)
    arr=aperm(arr,c(2,1,3))
    
    result=run(arr,workers)
    return(result)
  }
  
  # unregister parallel computing backend
  #unregister(cl)
  snow::stopCluster(cl)
  
  if(verbose)cat("writing results to disk\n")
  
  if(!separate&&out.bands>1){
    
    raster::rasterOptions(overwrite=T,datatype=out.dtype,setfileext=F)
    out.ras=raster::brick(raster::extent(ras.in),nrows=rows,ncols=cols,crs=raster::projection(ras.in),nl=out.bands)

    out.ras=raster::writeStart(out.ras,filename=out.name,format="ENVI",bandorder="BSQ",NAflag=na.out)
    i=1
    for(sub in out.arr){
      rows.sub=dim(sub)[1]
      sub=aperm(sub,c(2,1,3))
      dim(sub)=c(rows.sub*cols,out.bands)
      out.ras=raster::writeValues(out.ras,sub,i)
      i=i+rows.sub
    }
    out.ras=raster::writeStop(out.ras)

    # edit the band names of the resulting ENVI file to carry information of the computed measures
    # (i.e. the names of the workers, e.g. minimum, maximum, p05, etc.)
    hdrbands(paste(out.name,".hdr",sep=""),bandnames)
  }else if(out.bands==1){
    raster::rasterOptions(overwrite=T,datatype=out.dtype,setfileext=T)
    out.arr=abind::abind(out.arr,along=1)
    out.ras=raster::raster(out.arr,template=ras.in)
    raster::writeRaster(out.ras,filename=outnames[1],format="GTiff",bandorder="BSQ",NAflag=na.out,options=c("COMPRESS=NONE"))
  }else{
    raster::rasterOptions(overwrite=T,datatype=out.dtype,setfileext=T)
    for(i in seq(out.bands)){
      out.arr.sub=abind::abind(lapply(out.arr,function(x)x[,,i]),along=1)
      out.ras=raster::raster(out.arr.sub,template=ras.in)
      raster::writeRaster(out.ras,filename=outnames[i],format="GTiff",bandorder="BSQ",NAflag=na.out,options=c("COMPRESS=NONE"))
    }
    rm(out.arr.sub,out.ras)
  }
  # # transform the computed array into a raster object with geo information
  # out.ras=if(class(out.arr)=="array") raster::brick(out.arr) else raster::raster(out.arr)
  # raster::extent(out.ras)=raster::extent(ras.in)
  # raster::projection(out.ras)=raster::projection(ras.in)
  # 
  # # write the raster to file
  # 
  # raster::rasterOptions(overwrite=overwrite,datatype=out.dtype,setfileext=F)
  # if(!separate&&out.bands>1){
  #   raster::writeRaster(out.ras,filename=out.name,format="ENVI",bandorder="BSQ",NAflag=na.out)
  #   
  #   # edit the band names of the resulting ENVI file to carry information of the computed measures 
  #   # (i.e. the names of the workers, e.g. minimum, maximum, p05, etc.)
  #   hdrbands(paste(out.name,".hdr",sep=""),bandnames)
  # }else{
  #   opt=c("COMPRESS=NONE")
  #   for(i in seq(out.bands)){
  #     out.name.band=paste(tools::file_path_sans_ext(out.name),"_", bandnames[i], ".tif", sep="")
  #     raster::writeRaster(out.ras[[i]],filename=out.name.band,format="GTiff",bandorder="BSQ",NAflag=na.out,options=opt)
  #   }
  # }
  rm(out.arr)
  gc(verbose=F)
  if(verbose)cat(sprintf("elapsed time: %s\n",hms_span(start.time,Sys.time())))
}

