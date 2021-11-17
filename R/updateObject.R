
setMethod("updateObject", "GInteractions",
    function(object, ..., verbose=FALSE) {
        object@regions <- updateObject(object@regions, ..., verbose=verbose)
        callNextMethod()
    }
)

