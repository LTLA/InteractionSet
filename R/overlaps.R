###############################################################
# This defines the findOverlaps method for GInteractions objects.

.get_used <- function(gi)
# Gets the indices for all regions(gi) that are actually used as anchors.
{
    all1 <- anchors(gi, type="first", id=TRUE)
    all2 <- anchors(gi, type="second", id=TRUE)
    used <- logical(length(regions(gi)))
    used[all1] <- TRUE
    used[all2] <- TRUE
    which(used)
}

.fast_overlap <- function(gi, ranges, ..., gi.is.query=TRUE)
# Overlaps regions(gi) with ranges, but only for those regions used as anchors.
{
    regs <- regions(gi)
    subset <- .get_used(gi)
    if (length(subset)!=length(regs)) { regs <- regs[subset] }

    # Need this switch, as type="within" will vary depending on manner of query vs. subject.
    if (gi.is.query) {
        olap <- findOverlaps(regs, ranges, select="all", ...)
        gi.dex <- subset[queryHits(olap)]
        ranges.dex <- subjectHits(olap)
    } else {
        olap <- findOverlaps(ranges, regs, select="all", ...)
        gi.dex <- subset[subjectHits(olap)]
        ranges.dex <- queryHits(olap)
        o <- order(gi.dex, ranges.dex)
        gi.dex <- gi.dex[o]
        ranges.dex <- ranges.dex[o]
    }

    return(list(gi.dex=gi.dex, ranges.dex=ranges.dex))
}

.get_olap_bounds <- function(hits, N)
# Gets the first and last index of the overlap vector for each GI index.
{
    current.rle <- rle(hits)
    first.in.rle <- rep(1L, N)
    last.in.rle <- integer(N)
    cum.end <- cumsum(current.rle$lengths)
    first.in.rle[current.rle$values] <- cum.end - current.rle$lengths + 1L
    last.in.rle[current.rle$values] <- cum.end
    return(list(first=first.in.rle, last=last.in.rle))
}

.decode_region_mode <- function(use.region, possibilities=c("both", "first", "second"))
# Gets the index of possibilities, for easy entry into C++.
# 1 -> both, 2 -> first, 3 -> second.
# This needs to be sync'd with set_mode_values in overlaps.cpp.
{
    use.region <- match.arg(use.region, possibilities)
    match(use.region, possibilities)
}

.linear_olap_finder <- function(gi, ranges, ..., select, gi.is.query=TRUE, use.region="both")
# Identifies linear overlaps, with differing C++ function depending on whether
# all overlaps are desired, or counting should be performed, etc.
{
    olap <- .fast_overlap(gi, ranges, ..., gi.is.query=gi.is.query)
    a1 <- anchors(gi, type="first", id=TRUE)
    a2 <- anchors(gi, type="second", id=TRUE)

    # Getting all combinations of overlaps (zero-indexing for C code).
    bounds <- .get_olap_bounds(olap$gi.dex, length(regions(gi)))
    out <- .Call(cxx_linear_olaps, a1 - 1L, a2 - 1L, bounds$first - 1L, bounds$last,
                 olap$ranges.dex - 1L, length(ranges), .decode_region_mode(use.region),
                 select, gi.is.query)

    # Processing into a Hits object if required.
    final <- out
    if (select=="all") {
        if (!gi.is.query) {
            final <- Hits(out[[2]]+1L, out[[1]]+1L, length(ranges), length(gi), sort.by.query=TRUE)
        } else {
            final <- Hits(out[[1]]+1L, out[[2]]+1L, length(gi), length(ranges), sort.by.query=TRUE)
        }
    } else if (select=="count") {
        if (gi.is.query) {
            names(final) <- names(gi)
        } else {
            names(final) <- names(ranges)
        }
    } else {
        final <- final + 1L # get rid of zero-indexing for first/last/arbitrary.
    }
    return(final)
}

setMethod("findOverlaps", c(query="GInteractions", subject="Vector"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=TRUE, ..., use.region="both") {

        type <- match.arg(type)
        select <- match.arg(select)
        .linear_olap_finder(query, subject,
                    maxgap=maxgap, minoverlap=minoverlap, type=type,
                    select=select, ignore.strand=ignore.strand, ...,
                    gi.is.query=TRUE, use.region=use.region)
    }
)

setMethod("findOverlaps", c(query="Vector", subject="GInteractions"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=TRUE, ..., use.region="both") {

        type <- match.arg(type)
        select <- match.arg(select)
        .linear_olap_finder(subject, query,
                    maxgap=maxgap, minoverlap=minoverlap, type=type,
                    select=select, ignore.strand=ignore.strand, ...,
                    gi.is.query=FALSE, use.region=use.region)
    }
)

###############################################################

.paired_overlap_finder2 <- function(gi.left, gi.right, ..., select, use.region="both")
# Identifies overlaps between two GI objects (left and right),
# with different C++ functions for different behaviours as before.
{
    used2 <- .get_used(gi.right)
    olap <- .fast_overlap(gi.left, regions(gi.right)[used2], ..., gi.is.query=TRUE)
    olap$ranges.dex <- used2[olap$ranges.dex]

    left.a1 <- anchors(gi.left, type="first", id=TRUE)
    left.a2 <- anchors(gi.left, type="second", id=TRUE)
    left.bounds <- .get_olap_bounds(olap$gi.dex, length(regions(gi.left)))

    right.a1 <- anchors(gi.right, type="first", id=TRUE)
    o1 <- order(right.a1)
    right.bounds1 <- .get_olap_bounds(right.a1[o1], length(regions(gi.right)))

    right.a2 <- anchors(gi.right, type="second", id=TRUE)
    o2 <- order(right.a2)
    right.bounds2 <- .get_olap_bounds(right.a2[o2], length(regions(gi.right)))

    # Getting all 2D overlaps.
    npairs <- length(gi.right)
    out <- .Call(cxx_paired_olaps, left.a1 - 1L, left.a2 - 1L,
                 left.bounds$first - 1L, left.bounds$last, olap$ranges.dex - 1L,
                 right.bounds1$first - 1L, right.bounds1$last, o1 - 1L,
                 right.bounds2$first - 1L, right.bounds2$last, o2 - 1L,
                 .decode_region_mode(use.region, c("both", "same", "reverse")), select)

    # Deciding what output to return.
    final <- out
    if (select=="all") {
        final <- Hits(out[[1]]+1L, out[[2]]+1L, length(gi.left), length(gi.right), sort.by.query=TRUE)
    } else if (select!="count") {
        final <- final + 1L
        names(final) <- names(gi.left)
    }
    return(final)
}

setMethod("findOverlaps", c(query="GInteractions", subject="GInteractions"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=TRUE, ..., use.region="both") {

        type <- match.arg(type)
        select <- match.arg(select)
        .paired_overlap_finder2(query, subject, select=select,
                    maxgap=maxgap, minoverlap=minoverlap, type=type,
                    ignore.strand=ignore.strand, ..., use.region=use.region)
    }
)

setMethod("findOverlaps", c(query="GInteractions", subject="missing"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=TRUE, ...,
             drop.self=FALSE, drop.redundant=FALSE,
             use.region="both") {

        type <- match.arg(type)
        out <- findOverlaps(query, query, maxgap=maxgap, minoverlap=minoverlap,
                            type=type, select="all", ignore.strand=ignore.strand,
                            ..., use.region=use.region)
        select <- match.arg(select)
        IRanges:::process_self_hits(out, select, drop.self, drop.redundant)
    }
)

###############################################################
# This defines the countOverlaps method.

setMethod("countOverlaps", c(query="GInteractions", subject="Vector"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=TRUE, ..., use.region="both") {

        type <- match.arg(type)
        .linear_olap_finder(query, subject,
            maxgap=maxgap, minoverlap=minoverlap, type=type,
            ignore.strand=ignore.strand, ..., select="count",
            gi.is.query=TRUE, use.region=use.region)
    }
)

setMethod("countOverlaps", c(query="Vector", subject="GInteractions"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=TRUE, ..., use.region="both") {

        type <- match.arg(type)
        .linear_olap_finder(subject, query,
            maxgap=maxgap, minoverlap=minoverlap, type=type,
            ignore.strand=ignore.strand, ..., select="count",
            gi.is.query=FALSE, use.region=use.region)
    }
)

setMethod("countOverlaps", c(query="GInteractions", subject="GInteractions"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=TRUE, ..., use.region="both") {

        type <- match.arg(type)
        .paired_overlap_finder2(query, subject, maxgap=maxgap, minoverlap=minoverlap, type=type,
                    ignore.strand=ignore.strand, ..., select="count", use.region=use.region)
    }
)

setMethod("countOverlaps", c(query="GInteractions", subject="missing"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=TRUE, ..., use.region='both') {

        type <- match.arg(type)
        hits <- findOverlaps(query, maxgap=maxgap, minoverlap=minoverlap, type=type,
                             ignore.strand=ignore.strand, ..., use.region=use.region)
        ans <- countQueryHits(hits)
        names(ans) <- names(query)
        ans
    }
)

###############################################################
# Defining corresponding functions for InteractionSet objects.

for (sig in c("Vector", "GInteractions")) { # Need to specify GInteractions to avoid ambiguous redirects with Vector.

    # Query is an InteractionSet.
    setMethod("countOverlaps", c("InteractionSet", sig),
        function(query, subject, maxgap=-1L, minoverlap=0L,
                 type=c("any", "start", "end", "within", "equal"),
                 ignore.strand=TRUE, ..., use.region='both') {
            type <- match.arg(type)
            countOverlaps(interactions(query), subject, maxgap=maxgap, minoverlap=minoverlap,
                          type=type, ignore.strand=ignore.strand, ..., use.region=use.region)
        }
    )

    setMethod("findOverlaps", c("InteractionSet", sig),
        function(query, subject, maxgap=-1L, minoverlap=0L,
                 type=c("any", "start", "end", "within", "equal"),
                 select=c("all", "first", "last", "arbitrary"),
                 ignore.strand=TRUE, ..., use.region='both') {
            select <- match.arg(select)
            type <- match.arg(type)
            findOverlaps(interactions(query), subject, maxgap=maxgap,
                minoverlap=minoverlap, type=type, select=select,
                ignore.strand=ignore.strand, ..., use.region=use.region)
        }
    )

    # Subject is an InteractionSet.
    setMethod("countOverlaps", c(sig, "InteractionSet"),
        function(query, subject, maxgap=-1L, minoverlap=0L,
                 type=c("any", "start", "end", "within", "equal"),
                 ignore.strand=TRUE, ..., use.region='both') {
            type <- match.arg(type)
            countOverlaps(query, interactions(subject), maxgap=maxgap, minoverlap=minoverlap,
                          type=type, ignore.strand=ignore.strand, ..., use.region=use.region)
        }
    )

    setMethod("findOverlaps", c(sig, "InteractionSet"),
        function(query, subject, maxgap=-1L, minoverlap=0L,
                 type=c("any", "start", "end", "within", "equal"),
                 select=c("all", "first", "last", "arbitrary"),
                 ignore.strand=TRUE, ..., use.region='both') {
            select <- match.arg(select)
            type <- match.arg(type)
            findOverlaps(query, interactions(subject), maxgap=maxgap,
                         minoverlap=minoverlap, type=type, select=select,
                         ignore.strand=ignore.strand, ..., use.region=use.region)
        }
    )
}

# Both arguments are InteractionSet objects.
setMethod("countOverlaps", c("InteractionSet", "InteractionSet"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=TRUE, ..., use.region='both') {
        type <- match.arg(type)
        countOverlaps(interactions(query), interactions(subject), maxgap=maxgap, minoverlap=minoverlap,
                      type=type, ignore.strand=ignore.strand, ..., use.region=use.region)
    }
)

setMethod("findOverlaps", c("InteractionSet", "InteractionSet"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=TRUE, ..., use.region='both') {
        select <- match.arg(select)
        type <- match.arg(type)
        findOverlaps(interactions(query), interactions(subject), maxgap=maxgap, minoverlap=minoverlap,
                     type=type, select=select, ignore.strand=ignore.strand, ..., use.region=use.region)
    }
)

# Missing the subject.
setMethod("countOverlaps", c("InteractionSet", "missing"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=TRUE, ..., use.region='both') {
        type <- match.arg(type)
        countOverlaps(interactions(query),maxgap=maxgap, minoverlap=minoverlap, type=type,
                      ignore.strand=ignore.strand, ..., use.region=use.region)
    }
)

setMethod("findOverlaps", c("InteractionSet", "missing"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=TRUE, ..., use.region='both') {
        select <- match.arg(select)
        type <- match.arg(type)
        findOverlaps(interactions(query), maxgap=maxgap, minoverlap=minoverlap, type=type,
                     select=select, ignore.strand=ignore.strand, ..., use.region=use.region)
    }
)

###############################################################
# Defining overlapsAny for ContactMatrix objects.

setMethod("overlapsAny", c("ContactMatrix", "Vector"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
        type=c("any", "start", "end", "within", "equal"),
        ignore.strand=TRUE, ...) {
        a1 <- anchors(query, id=TRUE, type="row")
        a2 <- anchors(query, id=TRUE, type="column")

        is.used <- union(a1, a2)
        is.overlapped <- logical(length(regions(query)))
        is.overlapped[is.used] <- overlapsAny(regions(query)[is.used], subject, maxgap=maxgap, minoverlap=minoverlap,
                                              type=type, ignore.strand=ignore.strand, ...)
        return(list(row=is.overlapped[a1], column=is.overlapped[a2]))
})

# Use outer(output$row, output$column, "|" or "&") to get the logical area in the interaction space.
# Not sure it makes a great deal of sense to define 'findOverlaps' here.

for (sig in c("GInteractions", "InteractionSet")) {
    setMethod("overlapsAny", c("ContactMatrix", sig),
        function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=TRUE, ..., use.region='both') {

        # It's possible to do this more efficiently by avoiding instantiation of the full object.
        # But it would require a total re-implementation at the C++ level, which is a pain.
        row.a <- rep(anchors(query, type="row", id=TRUE), ncol(query))
        col.a <- rep(anchors(query, type="column", id=TRUE), each=nrow(query))
        new.query <- GInteractions(row.a, col.a, regions(query))
        out <- overlapsAny(new.query, subject, maxgap=maxgap, minoverlap=minoverlap, type=type,
                           ignore.strand=ignore.strand, ..., use.region=use.region)
        dim(out) <- dim(query)
        return(out)
    })
}

# Haven't defined the converse methods, as it's not clear whether you want to consider the entire
# interaction space in the ContactMatrix, or just the non-NA entries.

###############################################################
# End
