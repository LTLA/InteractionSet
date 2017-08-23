for (siglist in c("GInteractions", "ContactMatrix")) {
    # We use regions<- to force reordering if necessary.
    setMethod("trim", siglist, function(x, use.names=TRUE) {
        regions(x) <- trim(regions(x), use.names=use.names)
        return(x)
    })

    setMethod("resize", siglist, function(x, width, fix="start", use.names=TRUE, ...) {
        regions(x) <- resize(regions(x), width=width, fix=fix, use.names=use.names, ...)
        return(x)
    })

    setMethod("narrow", siglist, function(x, start=NA, end=NA, width=NA, use.names=TRUE) {
        regions(x) <- narrow(regions(x), start=start, end=end, width=width, use.names=use.names)
        return(x)        
    })  

    setMethod("shift", siglist, function(x, shift=0L, use.names=TRUE) {
        regions(x) <- shift(regions(x), shift=shift, use.names=use.names)
        return(x)
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

setMethod("width", "InteractionSet", function(x) { width(interactions(x)) })
