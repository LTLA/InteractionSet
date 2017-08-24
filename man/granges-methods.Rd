\name{GRanges methods}

\alias{trim,GInteractions-method}
\alias{trim,InteractionSet-method}
\alias{trim,ContactMatrix-method}

\alias{shift,GInteractions-method}
\alias{shift,InteractionSet-method}
\alias{shift,ContactMatrix-method}

\alias{flank,GInteractions-method}
\alias{flank,InteractionSet-method}
\alias{flank,ContactMatrix-method}

\alias{narrow,GInteractions-method}
\alias{narrow,InteractionSet-method}
\alias{narrow,ContactMatrix-method}

\alias{resize,GInteractions-method}
\alias{resize,InteractionSet-method}
\alias{resize,ContactMatrix-method}

\alias{width,GInteractions-method}
\alias{width,InteractionSet-method}
\alias{width,ContactMatrix-method}

\title{Equivalents of GRanges methods}

\description{Methods for GInteractions, InteractionSet and ContactMatrix that operate on the internal GenomicRanges.}

\usage{
\S4method{shift}{GInteractions}(x, shift=0L, use.names=TRUE) 
\S4method{shift}{InteractionSet}(x, shift=0L, use.names=TRUE) 
\S4method{shift}{ContactMatrix}(x, shift=0L, use.names=TRUE) 

\S4method{narrow}{GInteractions}(x, start=NA, end=NA, width=NA, use.names=TRUE)
\S4method{narrow}{InteractionSet}(x, start=NA, end=NA, width=NA, use.names=TRUE)
\S4method{narrow}{ContactMatrix}(x, start=NA, end=NA, width=NA, use.names=TRUE)

\S4method{resize}{GInteractions}(x, width, fix="start", use.names=TRUE, ...)
\S4method{resize}{InteractionSet}(x, width, fix="start", use.names=TRUE, ...)
\S4method{resize}{ContactMatrix}(x, width, fix="start", use.names=TRUE, ...)

\S4method{flank}{GInteractions}(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)
\S4method{flank}{InteractionSet}(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)
\S4method{flank}{ContactMatrix}(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)

\S4method{trim}{GInteractions}(x, use.names=FALSE)
\S4method{trim}{InteractionSet}(x, use.names=FALSE)
\S4method{trim}{ContactMatrix}(x, use.names=FALSE)

\S4method{width}{GInteractions}(x)
\S4method{width}{InteractionSet}(x)
\S4method{width}{ContactMatrix}(x)
}

\arguments{
\item{x}{A GInteractions, InteractionSet or ContactMatrix object.}
\item{shift, start, end, width, fix, both}{Further interaction-specific arguments to pass to the relevant GenomicRanges methods, see details.}
\item{use.names, ignore.strand, ...}{Further overall arguments to pass to the relevant GenomicRanges methods.}
}

\details{
The \code{trim}, \code{resize}, \code{narrow} and \code{shift} methods will apply the GenomicRanges methods with the same name to the anchor regions of \code{x}.
Note that this does \emph{not} refer to the \code{regions} slot of \code{x}.
For example, \code{resize(x, width)} will produce an object that is equivalent to running \code{resize} on \code{first(x)} and \code{second(x)} directly.
This is important when the arguments are vectors, as these are assumed to be parallel to the interactions rather than the \code{regions} slot.

After any modifications are applied, resorting may be performed to ensure that the entries in the \code{regions} slot of the output object are ordered.
This means that the order of the ranges in the \code{regions} slot may change between the input and output object.
However, the number and order of the \emph{interactions} will not change.

For GInteractions and InteractionSet objects, the \code{width} method will return a DataFrame with integer fields \code{anchor1} and \code{anchor2}.
Each row contains the width of the two anchor regions corresponding to each interaction.
For ContactMatrix objects, the method will return a list with the same fields (a different format is necessary as the lengths of \code{anchor1} and \code{anchor2} may be different).
}

\value{
Depending on the method, an object of the same class as \code{x}, or information regarding the genomic regions -- see Details.
}

\seealso{
\code{\link{trim,GenomicRanges-method}},
\code{\link{resize,GenomicRanges-method}},
\code{\link{narrow,GenomicRanges-method}},
\code{\link{shift,GenomicRanges-method}},
\code{\link{flank,GenomicRanges-method}},
\code{\link{width,IRanges-method}}
}

\examples{
example(GInteractions, echo=FALSE)

trim(gi)
width(gi)

# Functions are applied along the length of 'gi':
new.sizes <- round(runif(length(gi), 10, 50))
gi2 <- resize(gi, width=new.sizes)
gi2
resize(first(gi), width=new.sizes)
resize(second(gi), width=new.sizes)

# ... not along the length of 'regions(gi)': note the difference!
mod.regions <- suppressWarnings(resize(regions(gi), width=new.sizes))
mod.regions[anchors(gi, type="first", id=TRUE)] # not the same as first(gi2)!
mod.regions[anchors(gi, type="second", id=TRUE)] # not the same second(gi2)!

example(ContactMatrix, echo=FALSE)
trim(x)
width(x)
}

\author{
Aaron Lun, based on suggestions from Liz Ing-Simmons.
}

