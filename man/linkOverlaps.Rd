\name{linkOverlaps}
\alias{linkOverlaps}

\alias{linkOverlaps,GInteractions,Vector,Vector-method}
\alias{linkOverlaps,InteractionSet,Vector,Vector-method}
\alias{linkOverlaps,GInteractions,Vector,missing-method}
\alias{linkOverlaps,InteractionSet,Vector,missing-method}

\alias{linkOverlaps,GInteractions,Hits,Hits-method}
\alias{linkOverlaps,InteractionSet,Hits,Hits-method}
\alias{linkOverlaps,GInteractions,Hits,missing-method}
\alias{linkOverlaps,InteractionSet,Hits,missing-method}

\title{Link overlapping regions}

\description{
Identify interactions that link two sets of regions by having anchor regions overlapping one entry in each set.
}

\usage{
\S4method{linkOverlaps}{GInteractions,Vector,Vector}(query, subject1, subject2, ..., 
    ignore.strand=TRUE, use.region="both")

\S4method{linkOverlaps}{InteractionSet,Vector,Vector}(query, subject1, subject2, ..., 
    ignore.strand=TRUE, use.region="both")
}

\arguments{
\item{query}{A \link{GInteractions} or \link{InteractionSet} object.}

\item{subject1, subject2}{A \link[S4Vectors]{Vector} object defining a set of genomic intervals, such as a \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GRangesList}.
\code{subject2} may be missing.
Alternatively, both \code{subject1} and \code{subject2} may be Hits objects, see below.}
 
\item{..., ignore.strand}{Additional arguments to be passed to \code{\link[GenomicRanges]{findOverlaps}}.
Note that \code{ignore.strand=TRUE} by default, as genomic interactions are usually unstranded.}

\item{use.region}{
A string specifying which \code{query} regions should be used to overlap which \code{subject}.
Ignored if \code{subject2} is missing.
}
}

\details{
This function identifies all interactions in \code{query} where one anchor overlaps an entry in \code{subject1} and the other anchor overlaps an entry in \code{subject2}.
It is designed to be used to identify regions that are linked by interactions in \code{query}.
For example, one might specify genes as \code{subject1} and enhancers as \code{subject2}, to identify all gene-enhancer contacts present in \code{query}.
This is useful when the exact pairings between \code{subject1} and \code{subject2} are undefined.

The function returns a DataFrame specifying the index of the interaction in \code{query}; the index of the overlapped region in \code{subject1};
    and the index of the overlapped region in \code{subject2}.
If multiple regions in \code{subject1} and/or \code{subject2} are overlapping the anchor regions of a particular interaction,
    all combinations of two overlapping regions (one from each \code{subject*} set) are reported for that interaction.

By default, \code{use.region="both"} such that overlaps will be considered between any first/second interacting region in \code{query} and either \code{subject1} or \code{subject2}.
If \code{use.region="same"}, overlaps will only be considered between the first interacting region in \code{query} and entries in \code{subject1}, 
    and between the second interacting region and \code{subject2}.
The opposite applies with \code{use.region="reverse"}, where the first and second interacting regions are overlapped with \code{subject2} and \code{subject1} respectively.

If \code{subject2} is not specified, links within \code{subject1} are identified instead, i.e., \code{subject2} is set to \code{subject1}.
In such cases, the returned DataFrame is such that the first subject index is always greater than the second subject index, to avoid redundant permutations.
}

\section{Using Hits as input}{
Hits objects can be used for the \code{subject1} and \code{subject2} arguments.
These should be constructed using \code{\link[GenomicRanges]{findOverlaps}} with \code{regions(query)} as the query and the genomic regions of interest as the subject. 
For example, the calls below:
\preformatted{> linkOverlaps(query, subject1) # 1
> linkOverlaps(query, subject1, subject2) # 2}

will yield exactly the same output as:
\preformatted{> olap1 <- findOverlaps(regions(query), subject1)
> linkOverlaps(query, olap1) # 1
> olap2 <- findOverlaps(regions(query), subject2)
> linkOverlaps(query, olap1, olap2) # 2}

This is useful in situations where \code{regions(query)} and the genomic regions in \code{subject1} and \code{subject2} are constant across multiple calls to \code{linkOverlaps}.
In such cases, the overlaps only need to be calculated once, avoiding redundant work within \code{linkOverlaps}.
}

\value{
A DataFrame of integer indices indicating which elements of \code{query} link which elements of \code{subject1} and \code{subject2}.
}

\seealso{
\code{\link{findOverlaps,GInteractions,Vector-method}}
}

\author{
Aaron Lun
}

\examples{
example(GInteractions, echo=FALSE)

all.genes <- GRanges("chrA", IRanges(0:9*10, 1:10*10))
all.enhancers <- GRanges("chrB", IRanges(0:9*10, 1:10*10))

out <- linkOverlaps(gi, all.genes, all.enhancers)
head(out)

out <- linkOverlaps(gi, all.genes)
head(out)

# Same methods apply for InteractionSet objects.

example(InteractionSet, echo=FALSE)
out <- linkOverlaps(iset, all.genes, all.enhancers)
out <- linkOverlaps(iset, all.genes)
}
