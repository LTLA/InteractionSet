# Defines the match function.

#' @importFrom BiocGenerics match
#' @importFrom S4Vectors Hits
#' @importMethodsFrom S4Vectors match
setMethod("match", c("GInteractions", "GInteractions"), 
    function(x, table, nomatch = NA_integer_, incomparables = NULL, ...) {
        .strict_check(x, table)

        rx <- regions(x)
        rtab <- regions(table)
        if (length(rx)!=length(rtab) || any(rx!=rtab)) { 
            # Diverting to findOverlaps if the regions are not equal.
            return(findOverlaps(x, table, type="equal", select="first", use.region="same"))
        }

        # Using the Hits method, for convenience.
        nR <- length(rx)
        a1 <- Hits(anchor1(x), anchor2(x), nR, nR, sort.by.query=FALSE)
        a2 <- Hits(anchor1(table), anchor2(table), nR, nR, sort.by.query=FALSE)
        match(a1, a2, nomatch=nomatch, incomparables=incomparables, ...)
    }
)

setMethod("match", c("GInteractions", "InteractionSet"), 
    function(x, table, nomatch = NA_integer_, incomparables = NULL, ...) {
        match(x, interactions(table), nomatch=nomatch, incomparables=incomparables, ...)
    }
)

setMethod("match", c("InteractionSet", "GInteractions"), 
    function(x, table, nomatch = NA_integer_, incomparables = NULL, ...) {
        match(interactions(x), table, nomatch=nomatch, incomparables=incomparables, ...)
    }
)

setMethod("match", c("InteractionSet", "InteractionSet"), 
    function(x, table, nomatch = NA_integer_, incomparables = NULL, ...) {
        match(interactions(x), table, nomatch=nomatch, incomparables=incomparables, ...)
    }
)

# Defining the compare function (only for two GInteractions, it wouldn't make sense otherwise).

setMethod("pcompare",  c("GInteractions", "GInteractions"), function(x, y) { 
    .strict_check(x, y)
    a1.x <- anchor1(x)
    a2.x <- anchor2(x)
    a1.y <- anchor1(y)
    a2.y <- anchor2(y)

    rx <- regions(x)
    ry <- regions(y)
    if (length(rx)!=length(ry) || any(rx!=ry)) { # Coercing them to the same system, if they're not equal.
        collated <- .collate_GRanges(rx, ry)
        a1.x <- collated$indices[[1]][a1.x]
        a2.x <- collated$indices[[1]][a2.x]
        a1.y <- collated$indices[[2]][a1.y]
        a2.y <- collated$indices[[2]][a2.y]        
    }

    output1 <- a1.x - a1.y
    output2 <- a2.x - a2.y # explicitly formed for proper recycling.
    tied1 <- output1==0L
    output1[tied1] <- output2[tied1]
    return(output1)
})

.strict_check <- function(x, y) {
    if (is(x, "StrictGInteractions")!=is(y, "StrictGInteractions") 
        || is(x, "ReverseStrictGInteractions")!=is(y, "ReverseStrictGInteractions")) { 
        warning("comparison between GInteractions objects of different strictness")
    }
}
