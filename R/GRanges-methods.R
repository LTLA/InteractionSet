#######################################
# Internal functions for use in GRanges methods for paired regions.

.generate_regions <- function(x, FUN, args, other.args) {
    a1 <- anchor1(x)
    a2 <- anchor2(x)
    N1 <- length(a1)
    N2 <- length(a2)

    # Expanding all arguments to the requested length.
    arg1 <- arg2 <- args
    for (arg in names(args)) { 
        current <- args[[arg]]

        if (is.list(current)) { 
            c1 <- current[[1]]
            c2 <- current[[2]]            
        } else {
            c1 <- c2 <- current
        }

        c1 <- .expand_to_length(c1, N1)
        c2 <- .expand_to_length(c2, N2)
        arg1[[arg]] <- c1
        arg2[[arg]] <- c2
    }

    # Minimizing expansion by identifying unique modifications for each anchor region.
    reg <- regions(x)
    reg1 <- .apply_unique_mods(reg, a1, FUN, arg1, other.args)
    reg2 <- .apply_unique_mods(reg, a2, FUN, arg2, other.args)
    mod.reg <- .collate_GRanges(reg1, reg2)

    # Fiddling with the input/output ranges.
    unchecked_regions(x) <- mod.reg$ranges
    unchecked_anchor1(x) <- mod.reg$indices[[1]][a1]
    unchecked_anchor2(x) <- mod.reg$indices[[2]][a2]
    return(x)
}

.expand_to_length <- function(v, N) {
    tmp <- rep(v, length.out=N)
    tmp[] <- v # to trigger warning upon recycling a non-multiple vector.
    return(tmp)
}

.apply_unique_mods <- function(regions, index, FUN, args, other.args) {
    to.order <- c(list(index), args)
    o <- do.call(order, to.order)

    # Identifying unique modifications.
    mod <- integer(0) 
    if (length(o)) { 
        is.uniq <- FALSE
        for (element in to.order) {
            element <- element[o]
            is.uniq <- is.uniq | (element[-1]!=element[-length(o)])
        }
        mod <- o[c(TRUE, is.uniq)]
    }

    # Applying the modifications.
    mod.args <- lapply(args, "[", mod)
    to.mod <- index[mod]
    regions[to.mod] <- do.call(FUN, c(list(x=regions[to.mod]), mod.args, other.args))
    return(regions)
}

#######################################

for (siglist in c("GInteractions", "ContactMatrix")) {
    setMethod("trim", siglist, function(x, use.names=TRUE) {
        regions(x) <- trim(regions(x), use.names=use.names)
        return(x)
    })

    setMethod("resize", siglist, function(x, width, fix="start", use.names=TRUE, ...) {
        .generate_regions(x, FUN=resize, args=list(width=width, fix=fix), other.args=list(use.names=use.names, ...))
    })

    setMethod("narrow", siglist, function(x, start=NA, end=NA, width=NA, use.names=TRUE) {
        .generate_regions(x, FUN=narrow, args=list(start=start, end=end, width=width), other.args=list(use.names=use.names))
    })  

    setMethod("shift", siglist, function(x, shift=0L, use.names=TRUE) {
        .generate_regions(x, FUN=IRanges::shift, args=list(shift=shift), other.args=list(use.names=use.names))
    })

    setMethod("flank", siglist, function(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) {
        .generate_regions(x, FUN=flank, args=list(width=width, start=start), 
                          other.args=list(use.names=use.names, ignore.strand=ignore.strand))
    })
}

setMethod("width", "GInteractions", function(x) {
    w <- width(regions(x))          
    DataFrame(anchor1=w[anchor1(x)], anchor2=w[anchor2(x)])
})

setMethod("width", "ContactMatrix", function(x) {
    w <- width(regions(x))          
    list(anchor1=w[anchor1(x)], anchor2=w[anchor2(x)])
})

# Same methods for the interaction set.

setMethod("trim", "InteractionSet", function(x, use.names=TRUE) { 
    interactions(x) <- trim(interactions(x), use.names=use.names) 
    return(x)
})

setMethod("resize", "InteractionSet", function(x, width, fix="start", use.names=TRUE, ...) {
    interactions(x) <- resize(interactions(x), width=width, fix=fix, use.names=use.names, ...)
    return(x)
})

setMethod("narrow", "InteractionSet", function(x, start=NA, end=NA, width=NA, use.names=TRUE) {
    interactions(x) <- narrow(interactions(x), start=start, end=end, width=width, use.names=use.names)
    return(x)        
})  

setMethod("shift", "InteractionSet", function(x, shift=0L, use.names=TRUE) {
    interactions(x) <- shift(interactions(x), shift=shift, use.names=use.names)
    return(x)
})

setMethod("flank", "InteractionSet", function(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE) {
    interactions(x) <- flank(interactions(x), width=width, start=start, both=both, 
                             use.names=use.names, ignore.strand=ignore.strand)
    return(x)
})

setMethod("width", "InteractionSet", function(x) { width(interactions(x)) })
