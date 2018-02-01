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
  
  cpim=raster::canProcessInMemory(x)
  
  node=1
  maxnode=0
  timeout=1
  submitted=0
  received=0
  queue=seq(nodes)
  while(received<tr$n){
    socklist <- lapply(cl, function(x) x$con)
    ready=socketSelect(socklist,F,timeout=timeout)
    if (maxnode==nodes & sum(ready) > 0){
      n <- which.max(ready)
      d <- .recvData(cl[[n]])
      # print(sprintf('job %i finished by node %i', d$tag, n))
      if (! d$success ){ stop('cluster error') }
      
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
        res[i1:i2, ] <- d$value#d$value$value
        
      }else{
        if (received==1) {
          #nl is the number of layers here since each layer is written in one column
          nl <- NCOL(d$value) 
          if (nl > 1) {
            out <- raster::brick(out, nl=nl)
          }
          out <- raster::writeStart(out, filename=filename, ...)
        } 
        out <- raster::writeValues(out, d$value, tr$row[d$tag])
      }
      queue=c(queue[-n],n)
    }else{
      i=submitted+1
      if(i<=tr$n){
        node=queue[1]
        # print(sprintf('assigning job %i to node %i', i, node))
        .sendCall(cl[[node]], clusfun, list(fun, i), tag=i)
        queue=c(queue[-1],node)
        maxnode=max(maxnode,node)
        submitted=i
      }
    }
  }
  if(cpim){
    #set the values of the matrix into the raster brick and write it to disk
    out <- raster::setValues(out, res)
    if (filename != '') {
      out <- raster::writeRaster(out, filename, ...)
    }
  }else{
    out <- raster::writeStop(out)
  }
  raster::pbClose(pb)
  return(out)
}




clusterR2 <- function(x, fun, args=NULL, export=NULL, filename='', cl=NULL, m=2, ...) {
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
  
  
  if (!is.null(args)) {
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
  
  for (i in 1:nodes) {
    .sendCall(cl[[i]], clusfun, list(fun, i), tag=i)
  }
  
  #nothing is processed on memory here. This just decides whether the results are first gathered and then written or if the results are written to disk one by one
  if (raster::canProcessInMemory(x)) {
    
    index=1
    for (i in 1:tr$n) {
      raster::pbStep(pb, i)
      
      if(index>nodes)index=1
      
      #receive the result of one process; this blocks until the result is actually received
      #see here: https://stackoverflow.com/questions/39616127/checking-status-of-individual-node-of-cluster-created-with-parallelmakecluster
      d <- .recvData(cl[[index]])
      #d <- .recvOneData(cl)
      
      if (! d$success ) { 
        #if (! d$value$success ) { 
        print(d$value)
        #print(d$value$value)
        stop('cluster error') 
      }
      
      #initialize a brick and matrix to write the results to
      if (i ==1) {
        nl <- NCOL(d$value) 
        #nl <- NCOL(d$value$value)
        if (nl > 1) {
          out <- raster::brick(out, nl=nl)
        }
        res <- matrix(NA, nrow=raster::ncell(out), ncol=nl)
      } 
      
      j <- d$tag
      #j <- d$value$tag
      #insert the values of the current process into the matrix
      res[raster::cellFromRowCol(out, tr$row[j], 1):raster::cellFromRowCol(out, tr$row2[j], raster::ncol(out)), ] <- d$value#d$value$value
      ni <- nodes + i
      #print(sprintf('job %i finished by node %i', d$tag, i))
      #print(sprintf('job %i finished by node %i', d$value$tag, d$node))
      #trigger another process on the node which has just returned its result
      if (ni <= tr$n) {
        .sendCall(cl[[index]], clusfun, list(fun, ni), tag=ni)
      }
      index=index+1
    }
    #set the values of the matrix into the raster brick and write it to disk
    out <- raster::setValues(out, res)
    if (filename != '') {
      out <- raster::writeRaster(out, filename, ...)
    }
    raster::pbClose(pb)
    return(out)
    
  } else {
    
    index=1
    for (i in 1:tr$n) {
      raster::pbStep(pb, i)
      
      if(index>nodes)index=1
      
      d <- .recvData(cl[[index]])
      #d <- .recvOneData(cl)
      if (! d$success ) { stop('cluster error') }
      
      if (i ==1) {
        #nl is the number of layers here since each layer is written in one column
        nl <- NCOL(d$value) 
        if (nl > 1) {
          out <- raster::brick(out, nl=nl)
        }
        out <- raster::writeStart(out, filename=filename, ...)
      } 
      
      out <- raster::writeValues(out, d$value, tr$row[d$tag])
      ni <- nodes + i
      if (ni <= tr$n) {
        .sendCall(cl[[index]], clusfun, list(fun, ni), tag=ni)
      }
      index=index+1
    }
    out <- raster::writeStop(out)
    raster::pbClose(pb)
    return(out)
  }
}