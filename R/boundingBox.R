.boundingBox <- function(x, f) 
# As the name suggests, it computes a bounding box for groups of interactions in 'x'.
# Groups are specified by a factor 'f'.
#
# written by Aaron Lun
# created 2 December 2015
{
    if (missing(f)) { 
        f <- rep(1L, length(x))
    }
    if (length(f)!=length(x)) { 
        stop("length of 'f' must be equal to number of interactions")
    } 

    a1 <- anchors(x, type="first", id=TRUE)
    a2 <- anchors(x, type="second", id=TRUE)

    # Get the run lengths and values in a slightly convoluted way to handle factors (as rle() does not).
    by.f <- split(seq_along(f), f)
    o <- unlist(by.f)
    f.runs <- lengths(by.f)
    f.values <- names(by.f)

    # Extract the chromosomal information.
    x.reg <- regions(x)
    starts <- start(x.reg)
    ends <- end(x.reg)
    
    chrs <- seqnames(x.reg)
    ref.chr <- levels(chrs)
    chrs <- as.integer(chrs)

    # Using indices for easy comparison inside C++.
    bound1 <- .Call(cxx_get_box_bounds, f.runs, f.values, a1[o]-1L, chrs, starts, ends) 
    bound2 <- .Call(cxx_get_box_bounds, f.runs, f.values, a2[o]-1L, chrs, starts, ends) 

    gr1 <- GRanges(ref.chr[bound1[[1]]], IRanges(bound1[[2]], bound1[[3]]), seqinfo=seqinfo(x)) 
    gr2 <- GRanges(ref.chr[bound2[[1]]], IRanges(bound2[[2]], bound2[[3]]), seqinfo=seqinfo(x))
    out <- GInteractions(gr1, gr2)
    names(out) <- f.values
    return(out)
}

setMethod("boundingBox", "GInteractions", .boundingBox)
setMethod("boundingBox", "InteractionSet", .boundingBox)

