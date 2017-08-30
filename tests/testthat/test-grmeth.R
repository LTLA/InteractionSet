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

test_that("resizing works for GI, IS objects", {
    # Common resizing operation.
    new.size <- 30
    suppressWarnings(x2 <- resize(x, fix="center", width=new.size))
    ref <- resize(regions(x), fix="center", width=new.size)

    expect_identical(anchors(x2, type="first"), ref[anchors(x, type="first", id=TRUE)])
    expect_identical(anchors(x2, type="second"), ref[anchors(x, type="second", id=TRUE)])
    expect_false(is.unsorted(regions(x2)))

    # Checking that this works for ISet objects.    
    iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
    iset2 <- resize(iset, fix="center", width=new.size)
    expect_identical(class(iset2), class(iset))
    expect_identical(interactions(iset2), x2)
    
    # Variable resizing operation.
    new.size <- round(runif(Np, 10, 50))
    suppressWarnings(x2 <- resize(x, fix="center", width=new.size))
    expect_identical(anchors(x2, type="first"), resize(anchors(x, type="first"), width=new.size, fix="center"))
    expect_identical(anchors(x2, type="second"), resize(anchors(x, type="second"), width=new.size, fix="center"))
    expect_false(is.unsorted(regions(x2)))
    
    # Highly variable resizing.
    new.size1 <- round(runif(Np, 10, 50))
    new.size2 <- round(runif(Np, 10, 50))
    suppressWarnings(x2 <- resize(x, fix="center", width=list(new.size1, new.size2)))
    expect_identical(anchors(x2, type="first"), resize(anchors(x, type="first"), width=new.size1, fix="center"))
    expect_identical(anchors(x2, type="second"), resize(anchors(x, type="second"), width=new.size2, fix="center"))
    expect_false(is.unsorted(regions(x2)))  

    # Checking that a warning is raised when we recycle.
    expect_warning(x2 <- resize(x, fix="center", width=1:3), "not a multiple")
    expect_warning(x2 <- resize(x, fix="center", width=1:1001), "not a multiple")
})

test_that("resizing works for CM objects", {
    # Now for CM objects.
    new.size <- 30
    y <- inflate(x, "chrA", "chrB")
    y2 <- resize(y, fix="center", width=new.size)
    ref <- resize(regions(x), fix="center", width=new.size)

    expect_identical(class(y), class(y2))
    expect_identical(anchors(y2, type="row"), ref[anchors(y, type="row", id=TRUE)])
    expect_identical(anchors(y2, type="column"), ref[anchors(y, type="column", id=TRUE)])

    # Variable resizing.
    new.size1 <- round(runif(nrow(y), 10, 50))
    new.size2 <- round(runif(ncol(y), 10, 50))
    suppressWarnings(y2 <- resize(y, fix="center", width=list(new.size1, new.size2)))
    expect_identical(anchors(y2, type="row"), resize(anchors(y, type="row"), width=new.size1, fix="center"))
    expect_identical(anchors(y2, type="column"), resize(anchors(y, type="column"), width=new.size2, fix="center"))
    expect_false(is.unsorted(regions(y2)))  
})

# We won't bother to check variable values for the rest, as the underlying code is exactly the same as resize().

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

test_that("flanking works for all objects", {
    suppressWarnings(x2 <- flank(x, width=10))
    ref <- flank(regions(x), width=10)

    expect_identical(anchors(x2, type="first"), ref[anchors(x, type="first", id=TRUE)])
    expect_identical(anchors(x2, type="second"), ref[anchors(x, type="second", id=TRUE)])
    expect_false(is.unsorted(regions(x2)))
    
    iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
    iset2 <- flank(iset, width=10)
    expect_identical(class(iset2), class(iset))
    expect_identical(interactions(iset2), x2)

    # Now for CM objects.
    y <- inflate(x, "chrA", "chrB")
    y2 <- flank(y, width=10)

    expect_identical(class(y), class(y2))
    expect_identical(anchors(y2, type="row"), ref[anchors(y, type="row", id=TRUE)])
    expect_identical(anchors(y2, type="column"), ref[anchors(y, type="column", id=TRUE)])
})


test_that("width calculations work for all objects", {
    expect_identical(width(x), list(first=width(anchors(x, type="first")), second=width(anchors(x, type="second"))))
    
    iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
    expect_identical(width(iset), width(x))

    y <- inflate(x, "chrA", "chrB")
    expect_identical(width(y), list(row=width(anchors(y, type="row")), column=width(anchors(y, type="column"))))
})
