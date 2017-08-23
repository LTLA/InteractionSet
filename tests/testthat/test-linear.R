# Testing the linearize function for InteractionSet objects.
# library(InteractionSet); library(testthat); source("test-linear.R")

set.seed(200)
N <- 50
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))
all.regions <- unique(all.regions)
N <- length(all.regions)

Np <- 100
all.anchor1 <- sample(N, Np, replace=TRUE)
all.anchor2 <- sample(N, Np, replace=TRUE)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
colnames(counts) <- seq_len(Nlibs)
x <- InteractionSet(counts, GInteractions(all.anchor1, all.anchor2, all.regions))

####################################

test_that("conversion via linearization is working", {
    o <- order(all.regions)
    new.regions <- all.regions[o]
    new.pos <- integer(length(o))
    new.pos[o] <- seq_along(new.pos)
    new.anchor1 <- new.pos[all.anchor1]
    new.anchor2 <- new.pos[all.anchor2]

    for (interest in seq_len(N)) {
        cur.reg <- regions(x)[interest]
        out <- linearize(x, interest)
        chosen <- new.anchor1==interest | new.anchor2==interest
        expect_identical(assay(out), assay(x[chosen,]))
        expect_output(show(out), sprintf("class: RangedSummarizedExperiment 
dim: %i 4 
metadata(0):
assays(1): ''
rownames: NULL
rowData names(0):
colnames(4): 1 2 3 4
colData names(0):", sum(chosen)), fixed=TRUE)
    
        new.ranges <- anchors(x, type="first")
        not1 <- new.ranges!=cur.reg
        new.ranges[!not1] <- anchors(x, type="second")[!not1]
        new.ranges <- new.ranges[chosen]
        expect_identical(rowRanges(out), new.ranges)
    
        # Comparing with equivalent GRanges method.
        out2 <- linearize(x, cur.reg, type="equal")
        expect_equal(out, out2)
    
        # Also checking with GInteractions methods.
        out3 <- linearize(interactions(x), cur.reg, type="equal")
        expect_identical(out3, rowRanges(out2))
        out3 <- linearize(interactions(x), interest)
        expect_identical(out3, rowRanges(out2))
    }
})

test_that("linearization behaves correctly with silly inputs", {
    expect_identical(nrow(linearize(x, 0)), 0L)
    expect_identical(nrow(linearize(x[0,], 1)), 0L)
    expect_identical(nrow(suppressWarnings(linearize(x, GRanges("chrC", IRanges(1,1))))), 0L)
})

####################################

test_that("linearization works with multiple overlaps", {
    lenA <- max(end(regions(x)[seqnames(regions(x))=="chrA"]))
    x2 <- x
    mcols(x2)$Index <- seq_len(nrow(x2))
    out <- linearize(x2, GRanges("chrA", IRanges(1, lenA)))
    all.a <- anchors(x2[mcols(out)$Index])
    mcols(out)$Index <- NULL
   
    # Testing for proper merging of internal ranges. 
    is.same <- as.logical(seqnames(all.a$first)==seqnames(all.a$second))
    expect_identical(pmin(start(all.a$first), start(all.a$second))[is.same], start(rowRanges(out))[is.same])
    expect_identical(pmax(end(all.a$first), end(all.a$second))[is.same], end(rowRanges(out))[is.same])
    expect_identical(seqnames(all.a$first)[is.same], seqnames(rowRanges(out))[is.same])
    
    is.B1 <- as.logical(seqnames(all.a$first)=="chrB")
    expect_identical(all.a$first[is.B1], rowRanges(out)[is.B1])
    is.B2 <- as.logical(seqnames(all.a$second)=="chrB")
    expect_identical(all.a$second[is.B2], rowRanges(out)[is.B2])
    
    lenB <- max(end(regions(x)[seqnames(regions(x))=="chrB"]))
    expect_error(linearize(x2, GRanges(c("chrA", "chrB"), IRanges(1, c(lenA, lenB)))), "multi-chromosome sets of 'ref' are not supported")
})

test_that("linearization correctly ignores internal interactions", {
    lenA <- max(end(regions(x)[seqnames(regions(x))=="chrA"]))
    x2 <- x
    mcols(x2)$Index <- seq_len(nrow(x2))
    out <- linearize(x2, GRanges("chrA", IRanges(1, lenA)), internal=FALSE)
    all.a <- anchors(x2[mcols(out)$Index])
    mcols(out)$Index <- NULL
    
    is.diff <- as.logical(seqnames(all.a$first)!=seqnames(all.a$second))
    expect_true(all(is.diff))
    is.B1 <- as.logical(seqnames(all.a$first)=="chrB")
    expect_identical(all.a$first[is.B1], rowRanges(out)[is.B1])
    is.B2 <- as.logical(seqnames(all.a$second)=="chrB")
    expect_identical(all.a$second[is.B2], rowRanges(out)[is.B2])
})

# Testing what happens with region and interaction-specific metadata.

set.seed(100)
test_that("linearization preserves metadata", {
    out <- linearize(x, 1)
    mcols(x) <- DataFrame(Blah=runif(nrow(x)), Index=seq_len(nrow(x)))
    mcols(regions(x)) <- DataFrame(Foo=runif(length(regions(x))))

    out2 <- linearize(x, 1)
    expect_identical(colnames(mcols(out2)), c("Foo", "Blah", "Index"))
    expect_identical(mcols(out2)$Blah, mcols(x)$Blah[mcols(out2)$Index])

    a1 <- anchors(x, type="first", id=TRUE)
    a2 <- anchors(x, type="second", id=TRUE)
    index <- ifelse(a1!=1L, a1, a2)[mcols(out2)$Index]
    expect_identical(mcols(out2)$Foo, mcols(regions(x))$Foo[index])
})
