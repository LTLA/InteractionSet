\name{GRanges methods}

\alias{trim,GInteractions-method}
\alias{trim,InteractionSet-method}
\alias{trim,ContactMatrix-method}

\alias{shift,GInteractions-method}
\alias{shift,InteractionSet-method}
\alias{shift,ContactMatrix-method}

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
\S4method{trim}{GInteractions}(x, use.names=FALSE)
\S4method{trim}{InteractionSet}(x, use.names=FALSE)
\S4method{trim}{ContactMatrix}(x, use.names=FALSE)

\S4method{resize}{GInteractions}(x, width, fix="start", use.names=TRUE, ...)
\S4method{resize}{InteractionSet}(x, width, fix="start", use.names=TRUE, ...)
\S4method{resize}{ContactMatrix}(x, width, fix="start", use.names=TRUE, ...)

\S4method{narrow}{GInteractions}(x, start=NA, end=NA, width=NA, use.names=TRUE)
\S4method{narrow}{InteractionSet}(x, start=NA, end=NA, width=NA, use.names=TRUE)
\S4method{narrow}{ContactMatrix}(x, start=NA, end=NA, width=NA, use.names=TRUE)

\S4method{shift}{GInteractions}(x, shift=0L, use.names=TRUE) 
\S4method{shift}{InteractionSet}(x, shift=0L, use.names=TRUE) 
\S4method{shift}{ContactMatrix}(x, shift=0L, use.names=TRUE) 

\S4method{width}{GInteractions}(x)
\S4method{width}{InteractionSet}(x)
\S4method{width}{ContactMatrix}(x)
}

\arguments{
\item{x}{A GInteractions, InteractionSet or ContactMatrix object.}
\item{use.names, width, fix, start, end, shift, ...}{Further arguments to pass to the relevant GenomicRanges methods.}
}

\details{
The \code{trim}, \code{resize}, \code{narrow} and \code{shift} methods will apply the GenomicRanges methods with the same name to the entries of \code{regions(x)}.
The modified regions are then potentially resorted to ensure that the entries in the \code{regions} slot of the output object are ordered.
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
\code{\link{width,IRanges-method}}
}

\examples{
example(GInteractions, echo=FALSE)

trim(gi)
width(gi)

example(ContactMatrix, echo=FALSE)
trim(x)
width(x)
}

\author{
Aaron Lun, based on suggestions from Liz Ing-Simmons.
}

