\name{ContactMatrix overlaps}

\alias{overlapsAny,ContactMatrix,Vector-method}
\alias{overlapsAny,ContactMatrix,GInteractions-method}
\alias{overlapsAny,ContactMatrix,InteractionSet-method}

\title{Find overlaps between GRanges and a ContactMatrix}

\description{
Find overlaps for the set of regions representing the rows or columns of a ContactMatrix.
}

\usage{
\S4method{overlapsAny}{ContactMatrix,Vector}(query, subject, maxgap=-1L, minoverlap=0L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=TRUE, ...)

\S4method{overlapsAny}{ContactMatrix,GInteractions}(query, subject, maxgap=-1L, minoverlap=0L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=TRUE, ..., use.region="both")

\S4method{overlapsAny}{ContactMatrix,InteractionSet}(query, subject, maxgap=-1L, minoverlap=0L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=TRUE, ..., use.region="both")
}

\arguments{
\item{query}{A \link{ContactMatrix} object.}

\item{subject}{A \link[GenomicRanges]{GRanges}, \link[GenomicRanges]{GRangesList}, \link{GInteractions} or \link{InteractionSet} object.}

\item{maxgap, minoverlap, type, ..., ignore.strand}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}

\item{use.region}{
    A string specifying how the interacting regions are to be matched to row/column regions, see Details.
}
}

\value{
For \code{overlapsAny,ContactMatrix,GRanges-method}, a named list of two logical vectors is returned specifying the rows of columns of \code{query} overlapped by \code{subject}.

For the other \code{overlapsAny} methods, a logical matrix is returned indicating which entries in \code{query} are overlapped by \code{subject}.
}

\details{
When \code{subject} is a GRanges, overlaps are identified between the row regions of \code{query} and the regions in \code{subject}.
This is repeated for the column regions of \code{query}.
A list of two logical vectors is returned, specifying the rows and columns in \code{query} that are overlapped by any region in \code{subject}.
These vectors can be directly used to subset \code{query}.
Alternatively, they can be used in \code{outer} to generate a logical matrix for masking -- see Examples.

For all other classes of \code{subject}, two-dimensional overlaps are identified.
A logical matrix is returned indicating which entries in the ContactMatrix have overlaps with the specified interactions.
For any given entry, an overlap is only considered if the regions for the corresponding row and column both overlap anchor regions in the \code{subject}.
See \code{?"\link{findOverlaps,GInteractions,GInteractions-method}"} for more details.

If \code{use.region="both"}, overlaps between any row/column region and the first/second interacting region of \code{subject} are considered.
If \code{use.region="same"}, only overlaps between row regions and the first interacting regions, or between column regions and the second interacting regions are considered.
If \code{use.region="reverse"}, only overlaps between row regions and the second interacting regions, or between row regions and the first interacting regions are considered.
}

\examples{
example(ContactMatrix, echo=FALSE)

of.interest <- resize(sample(regions(x), 2), width=1, fix="center")
olap <- overlapsAny(x, of.interest)
olap
x[olap$row,]
x[,olap$column]
x[olap$row,olap$column]

keep <- outer(olap$row, olap$column, "|") # OR mask
temp <- as.matrix(x)
temp[!keep] <- NA

keep <- outer(olap$row, olap$column, "&") # AND mask
temp <- as.matrix(x)
temp[!keep] <- NA

# Two dimensional overlaps.
pairing <- GRangesList(first=regions(x), second=regions(x))
olap <- overlapsAny(x, pairing)
olap
olap <- overlapsAny(sort(x), pairing) # A bit prettier
olap
}

\author{
Aaron Lun
}
