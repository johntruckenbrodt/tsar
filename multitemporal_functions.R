################################################################
# function to unregister parallel computing backend
unregister=function(cluster){
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
  bandnames=paste("band names = {",paste(names,collapse=", "),"}",sep="")
  hdr=c(hdr[-(grep("band names",hdr):index)],bandnames)
  writeLines(hdr,hdrfile)
}
################################################################