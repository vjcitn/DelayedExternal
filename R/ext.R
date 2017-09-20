#> methods(class="DelayedArray")
# [1] -             !             [             [[            +            
# [6] acbind        anyNA         apply         arbind        as.array     
#[11] as.character  as.complex    as.data.frame as.integer    as.logical   
#[16] as.matrix     as.numeric    as.raw        as.vector     c            
#[21] cbind         coerce        DelayedArray  dim           dim<-        
#[26] dimnames      dimnames<-    dlogis        drop          is.finite    
#[31] is.infinite   is.na         is.nan        isEmpty       length       
#[36] Math          matrixClass   mean          names         names<-      
#[41] nchar         Ops           plogis        pmax2         pmin2        
#[46] qlogis        rbind         round         seed          show         
#[51] signif        split         splitAsList   Summary       t            
#[56] tolower       toupper       type          which        
#see '?methods' for accessing help and source code

# learning how to extend DelayedArray for external non-HDF5 store
# would likely move this code to restfulSE if successful

#' define container for H5S_dataset to use Delayed* interfaces
#' @import rhdf5client
#' @import DelayedArray
#' @slot seed instance of H5S_dataset
#' @slot dimnames list
#' @examples
#' library(rhdf5client)
#' if (!exists("bigec2")) example(H5S_source) # minimize hits to server
#' ds = bigec2[["assays"]] # methylation data
#' dem = DelayedHDF5SMatrix( ds )
#' dem
#' subset_seed_as_array(dem@seed, list(1:3, 1:4))
#' library(restfulSE)
#' data(banoSEMeta)
#' library(SummarizedExperiment)
#' assays(banoSEMeta)
#' assays(banoSEMeta)[[1]] = dem
#' assay(banoSEMeta[1:3,1:3])
#' assay(banoSEMeta[1:3,c("NA18498", "NA18501")])
#' @export
setClass("DelayedHDF5SMatrix", 
  representation(seed="H5S_dataset",
         dimnames="list"),
  contains="DelayedArray")
#' define a dim getter for H5S_dataset
#' @rdname DelayedHDF5SMatrix-class
#' @export
setMethod("dim", "H5S_dataset", function(x)
   rhdf5client:::internalDim(x))

#' define a dim getter for DelayedHDF5SMatrix
#' @rdname DelayedHDF5SMatrix-class
#' @export
setMethod("dim", "DelayedHDF5SMatrix", function(x)
    rev(rhdf5client:::internalDim(x@seed))) # to agree w/ R

#' define dimnames operations for DelayedHDF5SMatrix
#' @rdname DelayedHDF5SMatrix-class
#' @export
setMethod("dimnames", "DelayedHDF5SMatrix", function(x)
    x@dimnames)
#' define dimnames setter for DelayedHDF5SMatrix
#' @rdname DelayedHDF5SMatrix-class
#' @param value a list of dimnames
#' @export
setMethod("dimnames<-", "DelayedHDF5SMatrix", function(x,value)
    x@dimnames<-value)

#' define subset_seed_as_array, key infrastructure for Delayed* seed
#' @rdname DelayedHDF5SMatrix-class
#' @param seed instance of H5S_dataset
#' @param index list
#' @export
setMethod("subset_seed_as_array", "H5S_dataset",
   function(seed, index) {  # 'index' is expected to be an unnamed list of subscripts as positive integer
                            ### vectors, one vector per seed dimension. *Missing* list elements are allowed
                            ### and represented by NULLs.
   dims = rev(rhdf5client:::internalDim(seed))
   if (is.null(index[[1]])) index[[1]] = 1:dims[1]
   if (is.null(index[[2]])) index[[2]] = 1:dims[2]
   i = sproc(isplit(index[[1]])) # turn the request in R-idiom into a list of block-like queries
   j = sproc(isplit(index[[2]]))
   do.call(rbind, lapply(i, function(curi)  # could be improved -- collect column block for each row block
      do.call(cbind, lapply(j, function(curj) t(seed[curj, curi])))))
})
   
#' define bracket
#' @rdname DelayedHDF5SMatrix-class
#' @param i selection vector or missing
#' @param j selection vector or missing
#' @param drop logical ignored for now
#' @param \dots not used
#' @export
setMethod("[", "DelayedHDF5SMatrix", function(x, i, j, ..., drop=TRUE) {
#
# need to convert standard R idioms for subscripting to
# HDF5 slice specs -- can we reuse approach in rhdf5 -- seems we would need to use the C interface
#
    if (missing(i) && missing(j)) 
        return(x)
    if (!missing(i)) {
        if (is.character(i)) {
            if (is.null(dimnames(x)[[1]])) stop("need dimnames(x)[[1]] to resolve character row indexing")
            i = sproc(isplit(match(i, dimnames(x)[[1]])))
            }
        else if (is.numeric(i)) i = sproc(isplit(i))
        else stop("i must be character or numeric")
        }
    if (!missing(j)) {
        if (is.character(j)) {
            if (is.null(dimnames(x)[[2]])) stop("need dimnames(x)[[2]] to resolve character column indexing")
            j = sproc(isplit(match(j, dimnames(x)[[2]])))
            }
        else if (is.numeric(j)) j = sproc(isplit(j))
        else stop("j must be character or numeric")
        }
    do.call(rbind, lapply(i, function(curi)
      do.call(cbind, lapply(j, function(curj) t(x@seed[curj, curi])))))
})

#' constructor
#' @rdname DelayedHDF5SMatrix-class
#' @importFrom methods new
#' @param x instance of \code{\link[rhdf5client]{H5S_dataset-class}}
#' @param dimnames list or missing
#' @export
DelayedHDF5SMatrix = function(x, dimnames=list(NULL,NULL)) {
    new("DelayedHDF5SMatrix", seed=x, index=vector("list", length(dim(x))),
    metaindex=1:2,
    dimnames=dimnames)}

#' printer
#' @rdname DelayedHDF5SMatrix-class
#' @param object instance of \code{\link{DelayedHDF5SMatrix-class}}
#' @export
setMethod("show", "DelayedHDF5SMatrix", function(object) {
  cat("DelayedHDF5SMatrix of dim ", dim(object), "\n")
})


