## Class ##
setClass('SCEVAL', 
         slots = c(title = 'character', 
                   counts = 'list',
                   norm = 'list',
                   annot = 'data.frame',
                   redDim = 'list',
                   traject = 'list'))

## Methods ##
setGeneric("counts", function(x, ...){standardGeneric("counts")})
setGeneric("counts<-", function(x, ...){standardGeneric("counts<-")})
setGeneric("redDim", function(x, ...){standardGeneric("redDim")})
setGeneric("redDim<-", function(x, ...){standardGeneric("redDim<-")})
setGeneric("norm", function(x, ...){standardGeneric("norm")})
setGeneric("norm<-", function(x, ...){standardGeneric("norm<-")})
setGeneric("annot", function(x){standardGeneric("annot")})
setGeneric("traject", function(x, ...){standardGeneric("traject")})
setGeneric("traject<-", function(x, ...){standardGeneric("traject<-")})

setMethod("counts", "SCEVAL", function(x, tag=NULL){
  if(is.null(tag)){
    x@counts
  }else{
    x@counts[[tag]]  
  }
})
setMethod("counts<-", "SCEVAL", function(x, value=NULL){
  x@counts <- c(counts(x), value)
  return(x)
})
setMethod("redDim", "SCEVAL", function(x, tag=NULL){
  if(is.null(tag)){
    x@redDim
  }else{
    x@redDim[[tag]]  
  }
})
setMethod("redDim<-", "SCEVAL", function(x, value=NULL){
  x@redDim <- c(redDim(x), value)
  return(x)
})
setMethod("norm", "SCEVAL", function(x, tag=NULL){
  if(is.null(tag)){
    x@norm
  }else{
    x@norm[[tag]]  
  }
})
setMethod("norm<-", "SCEVAL", function(x, value=NULL){
  x@norm <- c(norm(x), value)
  return(x)
})
setMethod("annot", "SCEVAL", function(x){
  x@annot
})
setMethod("traject", "SCEVAL", function(x, tag=NULL){
  if(is.null(tag)){
    x@traject
  }else{
    x@traject[[tag]]  
  }
})
setMethod("traject<-", "SCEVAL", function(x, value=NULL){
  x@traject <- c(traject(x), value)
  return(x)
})

