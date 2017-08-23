# Tests the trimming of GInteractions objects.
# library(InteractionSet); library(testthat); source("test-grmeth.R")

set.seed(10000)

N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 50))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
x <- GInteractions(all.anchor1, all.anchor2, all.regions)

test_that("trimming works for all objects", {
    suppressWarnings(seqlengths(x) <- c(chrA=50, chrB=100))
    x2 <- trim(x)
    ref <- trim(regions(x))

    expect_identical(anchors(x2, type="first"), ref[anchors(x, type="first", id=TRUE)])
    expect_identical(anchors(x2, type="second"), ref[anchors(x, type="second", id=TRUE)])
    expect_false(is.unsorted(regions(x2)))
    
    iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
    iset2 <- trim(iset)
    expect_identical(class(iset2), class(iset))
    expect_identical(interactions(iset2), x2)

    # Now for CM objects.
    y <- inflate(x, "chrA", "chrB")
    y2 <- trim(y)

    expect_identical(class(y), class(y2))
    expect_identical(anchors(y2, type="row"), ref[anchors(y, type="row", id=TRUE)])
    expect_identical(anchors(y2, type="column"), ref[anchors(y, type="column", id=TRUE)])
})

test_that("resizing works for all objects", {
    new.size <- round(runif(N, 10, 50))
    suppressWarnings(x2 <- resize(x, fix="center", width=new.size))
    ref <- resize(regions(x), fix="center", width=new.size)

    expect_identical(anchors(x2, type="first"), ref[anchors(x, type="first", id=TRUE)])
    expect_identical(anchors(x2, type="second"), ref[anchors(x, type="second", id=TRUE)])
    expect_false(is.unsorted(regions(x2)))
    
    iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
    iset2 <- resize(iset, fix="center", width=new.size)
    expect_identical(class(iset2), class(iset))
    expect_identical(interactions(iset2), x2)

    # Now for CM objects.
    y <- inflate(x, "chrA", "chrB")
    y2 <- resize(y, fix="center", width=new.size)

    expect_identical(class(y), class(y2))
    expect_identical(anchors(y2, type="row"), ref[anchors(y, type="row", id=TRUE)])
    expect_identical(anchors(y2, type="column"), ref[anchors(y, type="column", id=TRUE)])
})

test_that("narrowing works for all objects", {
    suppressWarnings(x2 <- narrow(x, start=3))
    ref <- narrow(regions(x), start=3)

    expect_identical(anchors(x2, type="first"), ref[anchors(x, type="first", id=TRUE)])
    expect_identical(anchors(x2, type="second"), ref[anchors(x, type="second", id=TRUE)])
    expect_false(is.unsorted(regions(x2)))
    
    iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
    iset2 <- narrow(iset, start=3)
    expect_identical(class(iset2), class(iset))
    expect_identical(interactions(iset2), x2)

    # Now for CM objects.
    y <- inflate(x, "chrA", "chrB")
    y2 <- narrow(y, start=3)

    expect_identical(class(y), class(y2))
    expect_identical(anchors(y2, type="row"), ref[anchors(y, type="row", id=TRUE)])
    expect_identical(anchors(y2, type="column"), ref[anchors(y, type="column", id=TRUE)])
})

test_that("shifting works for all objects", {
    suppressWarnings(x2 <- shift(x, shift=10))
    ref <- shift(regions(x), shift=10)

    expect_identical(anchors(x2, type="first"), ref[anchors(x, type="first", id=TRUE)])
    expect_identical(anchors(x2, type="second"), ref[anchors(x, type="second", id=TRUE)])
    expect_false(is.unsorted(regions(x2)))
    
    iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
    iset2 <- shift(iset, shift=10)
    expect_identical(class(iset2), class(iset))
    expect_identical(interactions(iset2), x2)

    # Now for CM objects.
    y <- inflate(x, "chrA", "chrB")
    y2 <- shift(y, shift=10)

    expect_identical(class(y), class(y2))
    expect_identical(anchors(y2, type="row"), ref[anchors(y, type="row", id=TRUE)])
    expect_identical(anchors(y2, type="column"), ref[anchors(y, type="column", id=TRUE)])
})

test_that("width calculations work for all objects", {
    x2 <- trim(x)
    expect_identical(width(x2), DataFrame(anchor1=width(anchors(x2, type="first")), anchor2=width(anchors(x2, type="second"))))
    
    iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
    iset2 <- trim(iset)
    expect_identical(width(iset2), width(x2))

    y <- inflate(x, "chrA", "chrB")
    y2 <- trim(y)
    expect_identical(width(y2), list(anchor1=width(anchors(y2, type="row")), anchor2=width(anchors(y2, type="column"))))
})
