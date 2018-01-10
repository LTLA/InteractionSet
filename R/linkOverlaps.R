# This defines the linkOverlaps method. The idea is to use interactions
# to link any two sets of regions together; to figure out which ones are
# linked, and to identify the interactions linking them.

.linkOverlap <- function(query, subject1, subject2, ..., use.region="both") {
    a1 <- anchors(query, id=TRUE, type="first")
    a2 <- anchors(query, id=TRUE, type="second")
    nregs <- length(regions(query))

    olap1 <- .fast_overlap(query, subject1, ..., gi.is.query=TRUE)
    bounds1 <- .get_olap_bounds(olap1$gi.dex, nregs)
    nregs1 <- length(subject1)
    if (missing(subject2)) { 
        olap2 <- olap1
        bounds2 <- bounds1
        is.same <- TRUE
        nregs2 <- nregs1
    } else {
        olap2 <- .fast_overlap(query, subject2, ..., gi.is.query=TRUE)
        bounds2 <- .get_olap_bounds(olap2$gi.dex, nregs)
        is.same <- FALSE
        nregs2 <- length(subject2)
    }

    out <- .Call(cxx_expand_pair_links, a1 - 1L, a2 - 1L, 
                 bounds1$first - 1L, bounds1$last, olap1$ranges.dex - 1L, nregs1,
                 bounds2$first - 1L, bounds2$last, olap2$ranges.dex - 1L, nregs2,
                 is.same, .decode_region_mode(use.region, c("both", "same", "reverse")))
    return(DataFrame(query=out[[1]]+1L, subject1=out[[2]]+1L, subject2=out[[3]]+1L))    
}

setMethod("linkOverlaps", c("GInteractions", "Vector", "Vector"), .linkOverlap)
setMethod("linkOverlaps", c("GInteractions", "Vector", "missing"), .linkOverlap)
setMethod("linkOverlaps", c("InteractionSet", "Vector", "Vector"), .linkOverlap)
setMethod("linkOverlaps", c("InteractionSet", "Vector", "missing"), .linkOverlap)

# Defining another function for Hits inputs. This is more efficient in cases
# where overlaps between 'regions(query)' and the subjects have already been 
# computed, and you want to identify linkages based on pairing of regions,
# e.g., if the intervals being overlapped are constant across calls.

.linkOverlap_Hits <- function(query, subject1, subject2, use.region="both") {
    a1 <- anchors(query, id=TRUE, type="first")
    a2 <- anchors(query, id=TRUE, type="second")
    nregs <- length(regions(query))
    if (nLnode(subject1)!=nregs) {
        stop("'nLnode' of 'subject1' should be equal to 'length(region(query))'")
    }
 
    bounds1 <- .get_olap_bounds(queryHits(subject1), nregs)
    nregs1 <- nRnode(subject1)
    if (missing(subject2)) { 
        subject2 <- subject1
        bounds2 <- bounds1
        is.same <- TRUE
        nregs2 <- nregs1
    } else {
        if (nLnode(subject2)!=nregs) {
            stop("'nLnode' of 'subject2' should be equal to 'length(region(query))'")
        }
        bounds2 <- .get_olap_bounds(queryHits(subject2), nregs)
        is.same <- FALSE
        nregs2 <- nRnode(subject2)
    }

    out <- .Call(cxx_expand_pair_links, a1 - 1L, a2 - 1L, 
                 bounds1$first - 1L, bounds1$last, subjectHits(subject1) - 1L, nregs1,
                 bounds2$first - 1L, bounds2$last, subjectHits(subject2) - 1L, nregs2,
                 is.same, .decode_region_mode(use.region, c("both", "same", "reverse")))
    return(DataFrame(query=out[[1]]+1L, subject1=out[[2]]+1L, subject2=out[[3]]+1L))    
}

setMethod("linkOverlaps", c("GInteractions", "Hits", "Hits"), .linkOverlap_Hits)
setMethod("linkOverlaps", c("GInteractions", "Hits", "missing"), .linkOverlap_Hits)
setMethod("linkOverlaps", c("InteractionSet", "Hits", "Hits"), .linkOverlap_Hits)
setMethod("linkOverlaps", c("InteractionSet", "Hits", "missing"), .linkOverlap_Hits)

