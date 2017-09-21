\name{InteractionSetnews}
\title{InteractionSet News}
\encoding{UTF-8}

\section{Version 1.5.7}{\itemize{
\item Removed c() method for InteractionSet, rbind method for GInteractions.

\item Generalized ContactMatrix to allow any type of matrix-like object.

\item Supported inflate() for GInteractions without specifying fill.
Changed default fill for InteractionSet to missing.

\item Separated subsetting and combining documentation into two different pages. 

\item Added convenience wrappers for resize(), narrow() and shift() on the GenomicRanges slot in all objects.

\item Modified anchors() to return a list rather than GRangesList.

\item Modified width() for GInteractions to return a list() rather than a DataFrame().
Also changed names of list element.

\item Added the anchorIds() function for rapid extraction of anchor IDs.

\item Removed the requirement for identical regions in pcompare() and match().
}}

\section{Version 1.4.0}{\itemize{
\item Deprecated anchors<- in favour of anchorIds<-, to avoid confusion about 'value' type.

\item Added first(), second() functions for convenience.

\item Updates to documentation, tests.
}}

\section{Version 1.2.0}{\itemize{
\item
Allowed specification of NULL in row/column arguments to inflate().

\item
Fixes to tests and code in response to updates to BiocGenerics, S4Vectors.

\item 
Added CITATION to the F1000Res article.
}}

\section{Version 1.0.0}{\itemize{
\item
New package InteractionSet, containing base classes for genomic interaction data.
}}

