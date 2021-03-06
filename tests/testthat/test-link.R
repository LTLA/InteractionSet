# Tests the linkOverlaps method
# library(InteractionSet); library(testthat); source("test-link.R")

set.seed(9000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.strand <- sample(c("*", "+", "-"), N, replace=TRUE)
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends), strand=all.strand)

Np <- 100
all.anchor1 <- sample(N, Np, replace=TRUE)
all.anchor2 <- sample(N, Np, replace=TRUE)
x <- GInteractions(all.anchor1, all.anchor2, all.regions)

# Generating some random regions to test against.

Ngenes <- 10
gene.starts <- round(runif(Ngenes, 1, 100))
gene.ends <- gene.starts + round(runif(Ngenes, 5, 20))
gene.strand <- sample(c("*", "+", "-"), Ngenes, replace=TRUE)
gene.regions <- GRanges(rep(c("chrA", "chrB"), Ngenes/2), IRanges(gene.starts, gene.ends), strand=gene.strand)

Nenh <- 8
enh.starts <- round(runif(Nenh, 1, 100))
enh.ends <- enh.starts + round(runif(Nenh, 5, 20))
enh.strand <- sample(c("*", "+", "-"), Nenh, replace=TRUE)
enh.regions <- GRanges(rep(c("chrA", "chrB"), Nenh/2), IRanges(enh.starts, enh.ends), strand=enh.strand)

# Testing for ISets and GIs, in a range of scenarios.

for (cls in 1:2) { 
    if (cls==1L) {
        obj <- x
    } else {
        obj <- InteractionSet(matrix(0, Np, 4, dimnames=list(NULL, seq_len(4))), x)
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
            use.region <- "same"
        } else if (param==6L) {
            use.region <- "reverse"
        } else if (param==7L) {
            ignore.strand <- FALSE
        }

        test_that(sprintf("linking overlaps works between two regions (%i, %i)", cls, param), {
            # Getting the reference results by doing it manually via 'merge'.
            olap1.g <- findOverlaps(anchors(obj, type="first"), gene.regions, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, type=type)
            olap2.g <- findOverlaps(anchors(obj, type="second"), gene.regions, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, type=type)
            olap1.e <- findOverlaps(anchors(obj, type="first"), enh.regions, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, type=type)
            olap2.e <- findOverlaps(anchors(obj, type="second"), enh.regions, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, type=type)

            combo1 <- base::merge(olap1.g, olap2.e, by.x=1, by.y=1)
            combo2 <- base::merge(olap2.g, olap1.e, by.x=1, by.y=1)
            if (use.region=="both") { 
                combo <- rbind(combo1, combo2)
            } else if (use.region=="same") {
                combo <- combo1
            } else {
                combo <- combo2
            }

            combo <- DataFrame(combo)
            colnames(combo) <- c("query", "subject1", "subject2")
            rownames(combo) <- NULL

            is.dup <- duplicated(paste0(combo$query, ".", combo$subject1, ".", combo$subject2))
            combo <- combo[!is.dup,]
            o <- order(combo$query, combo$subject1, combo$subject2)
            combo <- combo[o,]
    
            expect_identical(combo, linkOverlaps(obj, gene.regions, enh.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                                 ignore.strand=ignore.strand, use.region=use.region))
        })

        test_that(sprintf("linking overlaps works between the same regions (%i, %i)", cls, param), {
            # More manually merging, but just with the gene regions.
            olap1.g <- findOverlaps(anchors(obj, type="first"), gene.regions, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, type=type)
            olap2.g <- findOverlaps(anchors(obj, type="second"), gene.regions, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, type=type)

            combo.S <- base::merge(olap1.g, olap2.g, by.x=1, by.y=1)
            combo.S <- DataFrame(combo.S)
            colnames(combo.S) <- c("query", "subject1", "subject2")
            rownames(combo.S) <- NULL

            new.s1 <- pmax(combo.S$subject1, combo.S$subject2)
            new.s2 <- pmin(combo.S$subject1, combo.S$subject2)
            combo.S$subject1 <- new.s1
            combo.S$subject2 <- new.s2
    
            is.dup <- duplicated(paste0(combo.S$query, ".", combo.S$subject1, ".", combo.S$subject2)) | is.na(combo.S$query)
            combo.S <- combo.S[!is.dup,]
            o <- order(combo.S$query, combo.S$subject1, combo.S$subject2)
            combo.S <- combo.S[o,]
    
            expect_identical(combo.S, linkOverlaps(obj, gene.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, 
                                                   ignore.strand=ignore.strand, use.region=use.region))
        })

        test_that(sprintf("linking overlaps works with Hits (%i, %i)", cls, param), {
            ref1 <- linkOverlaps(obj, gene.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region, ignore.strand=ignore.strand)
            olap1 <- findOverlaps(regions(obj), gene.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
            expect_identical(ref1, linkOverlaps(obj, olap1, use.region=use.region))
            
            ref2 <- linkOverlaps(obj, gene.regions, enh.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region, ignore.strand=ignore.strand)
            olap2 <- findOverlaps(regions(obj), enh.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand)
            expect_identical(ref2, linkOverlaps(obj, olap1, olap2, use.region=use.region))
        })

        if (!ignore.strand) {
            # Just checking that ignore.strand works as expected.
            test_that("linkOverlaps works as expected without strandedness", {
                olap1 <- linkOverlaps(obj, gene.regions, ignore.strand=TRUE)
                olap2 <- linkOverlaps(obj, gene.regions, enh.regions, ignore.strand=TRUE)

                strand(regions(obj)) <- "*"
                strand(gene.regions) <- "*"
                strand(enh.regions) <- "*"

                expect_identical(olap1, linkOverlaps(obj, gene.regions, ignore.strand=FALSE))
                expect_identical(olap2, linkOverlaps(obj, gene.regions, enh.regions, ignore.strand=FALSE))
           })
        }
    }

    # Testing against empty slots.
    test_that("linking overlaps behaves with empty inputs", {
        expect_identical(DataFrame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions[0], gene.regions))
        expect_identical(DataFrame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions, gene.regions[0]))
        expect_identical(DataFrame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions[0], gene.regions[0]))

        expect_error(linkOverlaps(obj, Hits()), "'nLnode' of 'subject1'")
        okay.hits <- Hits(nLnode=length(regions(obj)))
        expect_identical(DataFrame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, okay.hits))
        expect_error(linkOverlaps(obj, okay.hits, Hits()), "'nLnode' of 'subject2'")
        expect_identical(DataFrame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, okay.hits, okay.hits))
    })
}


