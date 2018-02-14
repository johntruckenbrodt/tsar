# Author: Robert J. Hijmans
# Date :  November 2011
# Version 1.0
# Licence GPL v3

# modified by John Truckenbrodt 2018

clusterR <- function(x, fun, args=NULL, export=NULL, filename='', cl=NULL, m=2, ...) {
  require(raster)
  require(parallel)
  if (is.null(cl)) {
    cl <- raster::getCluster()
    on.exit( raster::returnCluster() )
  }
  if (!is.null(export)) {
    parallel::clusterExport(cl, export)	
  }
  
  .sendCall <- eval(parse(text="parallel:::sendCall"))
  .recvData <- eval(parse(text="parallel:::recvData"))
  
  nodes <- length(cl)
  
  out <- raster::raster(x)
  
  m <- max(1, round(m))
  tr <- raster::blockSize(x, minblocks=nodes*m )
  if (tr$n < nodes) {
    nodes <- tr$n
  }
  
  tr$row2 <- tr$row + tr$nrows - 1
  pb <- raster::pbCreate(tr$n, label='clusterR', ...)
  
  if (!is.null(args)){
    stopifnot(is.list(args))
    clusfun <- function(fun, i) {
      r <- raster::crop(x, raster::extent(out, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=raster::ncol(out)))
      r <- do.call(fun, c(r, args))
      raster::getValues(r)
    }
  } else {
    clusfun <- function(fun, i) {
      r <- raster::crop(x, raster::extent(out, r1=tr$row[i], r2=tr$row2[i], c1=1, c2=raster::ncol(out)))
      r <- fun(r)
      raster::getValues(r)
    }
  }
  #function for cleaning memory on individual nodes
  cleanup=function(){
    rm(list=ls())
    rm(list=ls(name=.GlobalEnv))
    gc()
  }
  
  cpim=raster::canProcessInMemory(x)

  maxnode=0        #the maximum node id which a job has been assigned to
  timeout=1        #the time in seconds to wait for finished results from the nodes
  submitted=0      #the number of submitted jobs
  received=0       #the number of received jobs
  queue=seq(nodes) #a simple queue for job assignment
  while(received<tr$n){
    socklist <- lapply(cl, function(x) x$con)
    ready=socketSelect(socklist,F,timeout=timeout)
    if (maxnode==nodes & sum(ready) > 0){
      n <- which.max(ready)
      d <- parallel:::recvData(cl[[n]])
      # print(sprintf('job %i finished by node %i', d$tag, n))
      if (! d$success ){
        print(d$value)
        stop('cluster error')
      }
      
      received=received+1
      raster::pbStep(pb, received)
      
      if(cpim){
        #initialize a brick and matrix to write the results to
        if (received==1) {
          nl <- NCOL(d$value)
          if (nl > 1) {
            out <- raster::brick(out, nl=nl)
          }
          res <- matrix(NA, nrow=raster::ncell(out), ncol=nl)
        }
        j <- d$tag
        #insert the values of the current process into the matrix
        i1=raster::cellFromRowCol(out, tr$row[j], 1)
        i2=raster::cellFromRowCol(out, tr$row2[j], raster::ncol(out))
        res[i1:i2, ] <- d$value
        
      }else{
        if(received==1){
          #nl is the number of layers here since each layer is written in one column
          nl <- NCOL(d$value) 
          outlist <- lapply(filename,function(x)raster::writeStart(out, filename=x, ...))
        }
        for(i in seq(nl)){
          # print(sprintf('writing result %i.%i to file %s', d$tag, i, filename[i]))
          outlist[[i]] <- raster::writeValues(outlist[[i]], d$value[,i], tr$row[d$tag])
        }
      }
      #perform a memory cleanup on the node
      parallel:::sendCall(cl[[n]],cleanup,list())
      d <- parallel:::recvData(cl[[n]])
      gc()
      
      #move the node id to the start of the processing queue and reset the timeout to 1 sec
      queue=c(n,queue[-n])
      timeout=1
    }else{
      i=submitted+1
      if(i<=tr$n){
        #pick the first node in the processing queue and assign a job to it
        node=queue[1]
        # print(sprintf('assigning job %i of %i to node %i', i, tr$n, node))
        parallel:::sendCall(cl[[node]], clusfun, list(fun, i), tag=i)
        
        #move the node to the last position in the queue and possibly increase the maximum node id used
        queue=c(queue[-1],node)
        maxnode=max(maxnode,node)
        #if all nodes have been asigned to a job increase the timeout to infinity
        #i.e. don't assign another job until one node has finished its job
        if(maxnode==nodes){
          timeout=NULL
        }
        submitted=i
      }
    }
  }
  if(cpim){
    #set the values of the matrix into the raster brick and write it to disk
    out <- raster::setValues(out, res)
    if(filename != ''){
      out <- raster::writeRaster(out, filename, ...)
    }
  }else{
    #stop writing the individual files and delete their handlers
    while(length(outlist)>0){
      outlist[[1]] <- raster::writeStop(outlist[[1]])
      outlist[[1]]=NULL
      gc()
    }
  }
  raster::pbClose(pb)
  return(NULL)
  # return(if(class(out)=="list") raster::brick(out) else out)
}
