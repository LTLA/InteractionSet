# Testing the various overlap methods for InteractionSet objects.
# This implicitly tests the GInteractions methods, because one calls the other.
# library(InteractionSet); library(testthat); source("test-olap.R")

set.seed(300)
N <- 50
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
stranded <- sample(c("-", "+", "*"), N, replace=TRUE)
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-20, 20)), IRanges(all.starts, all.ends), strand=stranded)
all.regions <- unique(all.regions)
N <- length(all.regions)

Np <- 100
all.anchor1 <- sample(N, Np, replace=TRUE)
all.anchor2 <- sample(N, Np, replace=TRUE)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
colnames(counts) <- seq_len(Nlibs)
x <- InteractionSet(counts, GInteractions(all.anchor1, all.anchor2, all.regions))

#######################################################
#######################################################
# GInterations: 1D overlaps

### LOOP START
for (object in c("GRanges", "GRangesList")) {
###

if (object=="GRanges") {
    # GRanges.
    set.seed(301)
    Nq <- 10
    query.starts <- round(runif(Nq, 1, 100))
    query.ends <- query.starts + round(runif(Nq, 5, 20))
    query.strand <- sample(c("-", "+", "*"), Nq, replace=TRUE)
    query.regions <- GRanges(rep(c("chrA", "chrB"), Nq/2), IRanges(query.starts, query.ends), strand=query.strand)

} else if (object=="GRangesList") {
    # Paired regions in a GRangesList (don't be fooled by the GInteractions below!)
    set.seed(302)
    
    Nq2 <- 100
    query.starts <- round(runif(Nq2, 1, 100))
    query.ends <- query.starts + round(runif(Nq2, 5, 20))
    query.strand <- sample(c("-", "+", "*"), Nq2, replace=TRUE)
    query.regions1 <- GRanges(rep(c("chrA", "chrB"), Nq2/2), IRanges(query.starts, query.ends), strand=query.strand)

    query.starts <- round(runif(Nq2, 1, 100))
    query.ends <- query.starts + round(runif(Nq2, 5, 20))
    query.strand <- sample(c("-", "+", "*"), Nq2, replace=TRUE)
    query.regions2 <- GRanges(rep(c("chrA", "chrB"), Nq2/2), IRanges(query.starts, query.ends), strand=query.strand)

    query.regions <- pairs(GInteractions(query.regions1, query.regions2), as.grlist=TRUE)
}

for (param in seq_len(7)) {
    type <- "any"
    maxgap <- -1L
    minoverlap <- 0L
    use.region <- "both"
    ignore.strand <- TRUE

    if (param==2L) { 
        maxgap <- 10L
    } else if (param==3L) {
        minoverlap <- 10L
    } else if (param==4L) {
        type <- "within"
    } else if (param==5L) {
        use.region <- "first"
    } else if(param==6L) {
        use.region <- "second"
    } else if (param==7L) {
        ignore.strand <- FALSE
    }

    test_that(sprintf("findOverlaps works with 1D overlaps to %s (%i)", object, param), {
        # Manual calculation of overlaps from anchors().              
        expected1 <- findOverlaps(anchors(x, type="first"), query.regions, type=type, 
                                  maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
        expected2 <- findOverlaps(anchors(x, type="second"), query.regions, type=type, 
                                  maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)

        if (use.region=="both") {
            ref <- Hits(c(queryHits(expected1), queryHits(expected2)), c(subjectHits(expected1), subjectHits(expected2)),
                        nLnode=length(x), nRnode=length(query.regions), sort.by.query=TRUE)
        } else if (use.region=="first") {
            ref <- expected1
        } else if (use.region=="second") { 
            ref <- expected2
        }
       
        olap <- findOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                             use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(olap, unique(sort(ref)))

        # Checking how it behaves with different 'select' specifications. 
        for (select in c("first", "last", "arbitrary")) { 
            selected <- findOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                     select=select, use.region=use.region, ignore.strand=ignore.strand)
            if (select!="arbitrary") {
                expect_identical(selected, selectHits(olap, select))
            } else {
                is.okay <- !is.na(selected)
                expect_identical(is.okay, !is.na(selectHits(olap, select)))
                arbiter <- Hits(which(is.okay), selected[is.okay], length(x), length(query.regions), sort.by.query=TRUE)
                expect_true(all(!is.na(match(arbiter, olap))))
            }
        }

        # Checking that 'overlapsAny', 'countOverlaps' and 'subsetByOverlaps' works.
        count.lap <- countOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                   use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(count.lap, selectHits(olap, "count"))
        
        out <- overlapsAny(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                           use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(out, !is.na(selected))

        expect_equal(subsetByOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                      use.region=use.region, ignore.strand=ignore.strand), x[out,])

        # Checking the countOverlaps preserves names.
        x2 <- x
        rownames(x2) <- paste0("Y", seq_along(x2))
        count.lap2 <- countOverlaps(x2, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                    use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(count.lap2, setNames(count.lap, names(x2)))
    })

    test_that(sprintf("findOverlaps works with 1D overlaps to %s (reversed, %i)", object, param), {
        # Manual calculation of overlaps from anchors(), after swapping x and query.regions.
        rexpected1 <- findOverlaps(query.regions, anchors(x, type="first"), type=type, maxgap=maxgap, 
                                   minoverlap=minoverlap, ignore.strand=ignore.strand)
        rexpected2 <- findOverlaps(query.regions, anchors(x, type="second"), type=type, maxgap=maxgap, 
                                   minoverlap=minoverlap, ignore.strand=ignore.strand)

        if (use.region=="both") {
            rref <- Hits(c(queryHits(rexpected1), queryHits(rexpected2)), c(subjectHits(rexpected1), subjectHits(rexpected2)),
                         nLnode=length(query.regions), nRnode=length(x), sort.by.query=TRUE)
        } else if (use.region=="first") {
            rref <- rexpected1
        } else if (use.region=="second") {
            rref <- rexpected2
        }
        
        rolap <- findOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                              use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(rolap, unique(sort(rref)))
        if (type!="within") { 
            olap <- findOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                 use.region=use.region, ignore.strand=ignore.strand)
            expect_identical(rolap, t(olap))
        } 
   
        # Seeing how it behaves with 'select' 
        for (select in c("first", "last", "arbitrary")) { 
            selected <- findOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                     select=select, use.region=use.region, ignore.strand=ignore.strand)

            if (select!="arbitrary") {
                expect_identical(selected, selectHits(rolap, select))
            } else {
                is.okay <- !is.na(selected)
                expect_identical(is.okay, !is.na(selectHits(rolap, select)))
                arbiter <- Hits(which(is.okay), selected[is.okay], length(query.regions), length(x), sort.by.query=TRUE)
                expect_true(all(!is.na(match(arbiter, rolap))))
            }
        }
        
        # Seeing how the other *Overlaps functions behave. 
        count.lap <- countOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                   use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(count.lap, selectHits(rolap, "count"))

        out <- overlapsAny(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                           use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(out, !is.na(selected))
        expect_identical(subsetByOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                          use.region=use.region, ignore.strand=ignore.strand), query.regions[out])

        # Checking the countOverlaps preserves names.
        q2 <- query.regions
        names(q2) <- paste0("G", seq_along(q2))
        count.lap2 <- countOverlaps(q2, x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                    use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(count.lap2, setNames(count.lap, names(q2)))
    })

    if (!ignore.strand) {
        # Just checking that ignore.strand for findOverlaps,GenomicRanges does as expected.
        test_that("findOverlaps for 1D overlaps works as expected without strandedness", {
            ref.olap <- findOverlaps(x, query.regions, ignore.strand=TRUE)
            rev.olap <- findOverlaps(query.regions, x, ignore.strand=TRUE)
            strand(regions(x)) <- "*"
            strand(query.regions) <- "*"
            expect_identical(findOverlaps(x, query.regions, ignore.strand=FALSE), ref.olap)
            expect_identical(findOverlaps(query.regions, x, ignore.strand=FALSE), rev.olap)
        })
    }

    test_that(sprintf("findOverlaps behaves sensibly with silly %s", object), {
        expect_equal(findOverlaps(x[0,], query.regions), Hits(nLnode=0, nRnode=length(query.regions), sort.by.query=TRUE))
        expect_equal(findOverlaps(x, query.regions[0]), Hits(nLnode=length(x), nRnode=0L, sort.by.query=TRUE))
        expect_equal(findOverlaps(query.regions, x[0]), Hits(nRnode=0, nLnode=length(query.regions), sort.by.query=TRUE))
        expect_equal(findOverlaps(query.regions[0], x), Hits(nRnode=length(x), nLnode=0L, sort.by.query=TRUE))
        expect_equal(findOverlaps(x[0,], query.regions[0]), Hits(nLnode=0L, nRnode=0L, sort.by.query=TRUE))
        expect_equal(findOverlaps(query.regions[0,], x[0]), Hits(nLnode=0L, nRnode=0L, sort.by.query=TRUE))
        
        expect_equal(findOverlaps(x[0,], query.regions, select="first"), integer(0))
        expect_equal(findOverlaps(query.regions[0], x, select="first"), integer(0))
        expect_equal(findOverlaps(x, query.regions[0], select="first"), rep(as.integer(NA), nrow(x)))
        expect_equal(findOverlaps(query.regions, x[0], select="first"), rep(as.integer(NA), length(query.regions)))
        expect_equal(overlapsAny(x[0,], query.regions), logical(0))
        expect_equal(overlapsAny(query.regions[0]), logical(0))
        expect_equal(overlapsAny(x, query.regions[0]), logical(nrow(x)))
        expect_equal(overlapsAny(query.regions, x[0]), logical(length(query.regions)))
    })
}

### LOOP END
}
###

#######################################################
#######################################################
# GInterations/InteractionSet: 2D overlaps

### LOOP START
for (object in c("InteractionSet", "GInteractions")) {
### 

set.seed(303)
N2 <- 75
next.starts <- round(runif(N2, 1, 100))
next.ends <- next.starts + round(runif(N2, 5, 20))
next.stranded <- sample(c("-", "+", "*"), N2, replace=TRUE)
next.regions <- GRanges(rep(c("chrA", "chrB"), c(N2-20, 20)), IRanges(next.starts, next.ends), strand=next.stranded)
next.regions <- unique(next.regions)
N2 <- length(next.regions)

Np2 <- 120
next.anchor1 <- sample(N2, Np2, replace=TRUE)
next.anchor2 <- sample(N2, Np2, replace=TRUE)
counts <- matrix(rpois(Np2*Nlibs, lambda=10), ncol=Nlibs)
colnames(counts) <- seq_len(Nlibs)

if (object=="InteractionSet") {
    x2 <- InteractionSet(counts, GInteractions(next.anchor1, next.anchor2, next.regions))
} else {
    x2 <- GInteractions(next.anchor1, next.anchor2, next.regions)
}

for (param in seq_len(6)) {
    type <- "any"
    maxgap <- -1L
    minoverlap <- 0L
    use.region <- "both"
    ignore.strand <- TRUE

    if (param==2L) { 
        maxgap <- 10L
    } else if (param==3L) {
        minoverlap <- 10L
    } else if (param==4L) {
        type <- "within"
    } else if (param==5L) {
        use.region <- "same"
    } else if(param==6L) {
        use.region <- "reverse"
    } else if (param==7L) {
        ignore.strand <- FALSE
    }

    test_that(sprintf("findOverlaps works with 2D overlaps to %s (%i)", object, param), {
        # Manually calculating overlaps based on anchors().
        expected1.A <- findOverlaps(anchors(x, type="first"), anchors(x2, type="first"), type=type, 
                                    maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
        expected1.B <- findOverlaps(anchors(x, type="first"), anchors(x2, type="second"), type=type, 
                                    maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
        expected2.B <- findOverlaps(anchors(x, type="second"),anchors(x2, type="first"), type=type, # Yes, the A/B switch in order is deliberate.
                                    maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand) 
        expected2.A <- findOverlaps(anchors(x, type="second"), anchors(x2, type="second"), type=type, 
                                    maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
    
        expected1.A <- paste0(queryHits(expected1.A), ".", subjectHits(expected1.A), ".A")
        expected1.B <- paste0(queryHits(expected1.B), ".", subjectHits(expected1.B), ".B")
        expected2.A <- paste0(queryHits(expected2.A), ".", subjectHits(expected2.A), ".A")
        expected2.B <- paste0(queryHits(expected2.B), ".", subjectHits(expected2.B), ".B")
        if (use.region=="both") {
            expected1 <- c(expected1.A, expected1.B)
            expected2 <- c(expected2.A, expected2.B)
        } else if (use.region=="same") {
            expected1 <- expected1.A
            expected2 <- expected2.A
        } else {
            expected1 <- expected1.B
            expected2 <- expected2.B
        }
    
        expected <- intersect(expected1, expected2)
        harvest <- do.call(rbind, strsplit(expected, "\\."))
        ref <- Hits(as.integer(harvest[,1]), as.integer(harvest[,2]), nLnode=length(x), nRnode=length(x2), sort.by.query=TRUE)
        ref <- sort(unique(ref))
        
        olap <- findOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                             use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(olap, ref)
    
        # Checking 'select' arguments.
        for (select in c("first", "last", "arbitrary")) { 
            selected <- findOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                     use.region=use.region, select=select, ignore.strand=ignore.strand)

            if (select!="arbitrary") {
                expect_identical(selected, selectHits(olap, select))
            } else {
                is.okay <- !is.na(selected)
                expect_identical(is.okay, !is.na(selectHits(olap, select)))
                arbiter <- Hits(which(is.okay), selected[is.okay], length(x), length(x2), sort.by.query=TRUE)
                expect_true(all(!is.na(match(arbiter, olap))))
            }
        }
    
        # Checking that 'overlapsAny', 'countOverlaps' and 'subsetByOverlaps' works.
        count.lap <- countOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                   use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(count.lap, selectHits(olap, "count"))

        out <- overlapsAny(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                           use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(out, !is.na(selected))
        expect_equal(subsetByOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                      use.region=use.region, ignore.strand=ignore.strand), x[out,])
    })

    # Note: No need to flip, it's the same method.
    # However, we do test self-overlaps.

    test_that(sprintf("findOverlaps works with 2D self-overlaps (%i)", param), {
        # Manually calculating overlaps based on anchors().
        expected1.A <- findOverlaps(anchors(x, type="first"), anchors(x, type="first"), type=type, 
                                    maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
        expected1.B <- findOverlaps(anchors(x, type="first"), anchors(x, type="second"), type=type, 
                                    maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
        expected2.B <- findOverlaps(anchors(x, type="second"), anchors(x, type="first"), type=type, # Again; yes, the A/B switch is deliberate.
                                    maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand) 
        expected2.A <- findOverlaps(anchors(x, type="second"), anchors(x, type="second"), type=type, 
                                    maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
    
        expected1.A <- paste0(queryHits(expected1.A), ".", subjectHits(expected1.A), ".A")
        expected1.B <- paste0(queryHits(expected1.B), ".", subjectHits(expected1.B), ".B")
        expected2.A <- paste0(queryHits(expected2.A), ".", subjectHits(expected2.A), ".A")
        expected2.B <- paste0(queryHits(expected2.B), ".", subjectHits(expected2.B), ".B")
        if (use.region=="both") {
            expected1 <- c(expected1.A, expected1.B)
            expected2 <- c(expected2.A, expected2.B)
        } else if (use.region=="same") {
            expected1 <- expected1.A
            expected2 <- expected2.A
        } else {
            expected1 <- expected1.B
            expected2 <- expected2.B
        }
    
        expected <- intersect(expected1, expected2)
        harvest <- do.call(rbind, strsplit(expected, "\\."))
        ref <- SelfHits(as.integer(harvest[,1]), as.integer(harvest[,2]), nnode=length(x), sort.by.query=TRUE)
        ref <- sort(unique(ref))
        
        self.olap <- findOverlaps(x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                  use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(self.olap, ref)

        # Checking that the function responds to the various drop.* arguments 
        self.olap.noself <- findOverlaps(x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region,
                                         drop.self=TRUE, ignore.strand=ignore.strand)
        expect_identical(self.olap.noself, self.olap[!isSelfHit(self.olap)])

        self.olap.nored <- findOverlaps(x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region,
                                        drop.redundant=TRUE, ignore.strand=ignore.strand)
        expect_identical(self.olap.nored, self.olap[!isRedundantHit(self.olap)])
    
        # Checking the behaviour of 'select'.
        for (select in c("first", "last", "arbitrary")) { 
            selected <- findOverlaps(x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                     use.region=use.region, select=select, ignore.strand=ignore.strand)

            if (select!="arbitrary") {
                expect_identical(selected, selectHits(self.olap, select))
            } else {
                is.okay <- !is.na(selected)
                expect_identical(is.okay, !is.na(selectHits(self.olap, select)))
                arbiter <- Hits(which(is.okay), selected[is.okay], length(x), length(x), sort.by.query=TRUE)
                expect_true(all(!is.na(match(arbiter, self.olap))))
            }
        }
    
        # Checking that 'overlapsAny', 'countOverlaps' and 'subsetByOverlaps' works.
        count.lap <- countOverlaps(x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                   use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(count.lap, selectHits(self.olap, "count"))

        out <- overlapsAny(x, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                           use.region=use.region, ignore.strand=ignore.strand)
        expect_identical(out, !is.na(selected))
    })

    if (!ignore.strand) {
        # Just checking that ignore.strand for findOverlaps,GenomicRanges does as expected.
        test_that("findOverlaps for 2D overlaps works as expected without strandedness", {
            ref.olap <- findOverlaps(x, x2, ignore.strand=TRUE)
            self.olap <- findOverlaps(x, ignore.strand=TRUE)
            strand(regions(x)) <- "*"
            strand(regions(x2)) <- "*"
            expect_identical(findOverlaps(x, x2, ignore.strand=FALSE), ref.olap)
            expect_identical(findOverlaps(x, ignore.strand=FALSE), self.olap)
        })
    }
}

### LOOP END
}
###

test_that("findOverlaps in 2D-mode works with silly inputs", {
    expect_equal(findOverlaps(x[0], x2), Hits(nLnode=0, nRnode=Np2, sort.by.query=TRUE))
    expect_equal(findOverlaps(x, x2[0]), Hits(nRnode=0, nLnode=nrow(x), sort.by.query=TRUE))
    expect_equal(findOverlaps(x[0], x2[0]), Hits(nLnode=0L, nRnode=0L, sort.by.query=TRUE))
    
    expect_equal(findOverlaps(x[0], x2, select="first"), integer(0))
    expect_equal(findOverlaps(x, x2[0], select="first"), rep(as.integer(NA), nrow(x)))
    expect_equal(overlapsAny(x[0], x2), logical(0))
    expect_equal(overlapsAny(x, x2[0]), logical(nrow(x)))
})

#######################################################
#######################################################
# ContactMatrix overlaps

set.seed(304)
Nr <- 100
Nc <- 200
all.anchor1 <- sample(N, Nr, replace=TRUE)
all.anchor2 <- sample(N, Nc, replace=TRUE)
counts <- matrix(rpois(Nr*Nc, lambda=10), Nr, Nc)
x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)

Nq <- 6
query.starts <- round(runif(Nq, 1, 100))
query.ends <- query.starts + round(runif(Nq, 5, 20))
query.regions <- GRanges(rep(c("chrA", "chrB"), Nq/2), IRanges(query.starts, query.ends))

test_that("overlapsAny for CM objects works with 1D overlaps", {
    olap <- overlapsAny(x, query.regions)
    expect_identical(olap, list(row=overlapsAny(anchors(x, type="row"), query.regions),
                                column=overlapsAny(anchors(x, type="column"), query.regions)))
    
    olap <- overlapsAny(x, query.regions, type="within")
    expect_identical(olap, list(row=overlapsAny(anchors(x, type="row"), query.regions, type="within"),
                                column=overlapsAny(anchors(x, type="column"), query.regions, type="within")))

    olap <- overlapsAny(x, query.regions, ignore.strand=FALSE)
    expect_identical(olap, list(row=overlapsAny(anchors(x, type="row"), query.regions, ignore.strand=FALSE),
                                column=overlapsAny(anchors(x, type="column"), query.regions, ignore.strand=FALSE)))
    
    expect_equal(overlapsAny(x[0,], query.regions), list(row=logical(0), column=overlapsAny(anchors(x, type="column"), query.regions)))
    expect_equal(overlapsAny(x[,0], query.regions), list(row=overlapsAny(anchors(x, type="row"), query.regions), column=logical(0)))
    expect_equal(overlapsAny(x[0,0], query.regions), list(row=logical(0), column=logical(0)))
})

test_that("overlapsAny for CM objects works with 2D overlaps", {
    Nq3 <- 10
    query.starts <- round(runif(Nq3, 1, 100))
    query.ends <- query.starts + round(runif(Nq3, 5, 20))
    query.regions1 <- GRanges(rep(c("chrA", "chrB"), Nq3/2), IRanges(query.starts, query.ends))
    query.starts <- round(runif(Nq3, 1, 100))
    query.ends <- query.starts + round(runif(Nq3, 5, 20))
    query.regions2 <- GRanges(rep(c("chrA", "chrB"), Nq3/2), IRanges(query.starts, query.ends))
    pairing <- GInteractions(query.regions1, query.regions2)
    
    olap <- overlapsAny(x, pairing)
    temp.iset <- deflate(x, collapse=TRUE)
    ref <- overlapsAny(temp.iset, pairing) # no NAs, everyone's represented here.
    ref <- inflate(temp.iset, anchors(x, type="row", id=TRUE), anchors(x, type="column", id=TRUE), fill=ref)
    expect_identical(olap, unname(as.matrix(as.matrix(ref))))
    
    olap <- overlapsAny(x, pairing, type="within")
    ref <- overlapsAny(temp.iset, pairing, type="within")
    ref <- inflate(temp.iset, anchors(x, type="row", id=TRUE), anchors(x, type="column", id=TRUE), fill=ref)
    expect_identical(olap, unname(as.matrix(as.matrix(ref))))

    olap <- overlapsAny(x, pairing, ignore.strand=FALSE)
    ref <- overlapsAny(temp.iset, pairing, ignore.strand=FALSE)
    ref <- inflate(temp.iset, anchors(x, type="row", id=TRUE), anchors(x, type="column", id=TRUE), fill=ref)
    expect_identical(olap, unname(as.matrix(as.matrix(ref))))
    
    olap <- overlapsAny(x, pairing)
    new.gi <- GInteractions(query.regions1, query.regions2)
    expect_identical(olap, overlapsAny(x, new.gi))
    
    new.iset <- InteractionSet(matrix(runif(10), dimnames=list(NULL, 1)), new.gi)
    expect_identical(olap, overlapsAny(x, new.iset))
    
    # Checking with order enforcement.
    olap <- overlapsAny(x, pairing, use.region="same")
    expect_identical(olap, overlapsAny(x, new.gi, use.region="same"))
    
    olap <- overlapsAny(x, pairing, use.region="reverse")
    expect_identical(olap, overlapsAny(x, new.gi, use.region="reverse"))
})

