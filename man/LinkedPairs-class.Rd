\name{LinkedPairs}
\alias{LinkedPairs-class}% the class
\alias{[.LinkedPairs}
\alias{print.LinkedPairs}
\title{
Tables of where syntenic hits link pairs of genes
}
\description{
Syntenic blocks describe where order is shared between two sequences. These blocks are made up of exact match hits. These hits can be overlayed on the locations of sequence features to clearly illustrate where exact sequence similarity is shared between pairs of sequence features.
}
\usage{
\method{print}{LinkedPairs}(x,
      quote = FALSE,
      right = TRUE,
      \dots)
}
\arguments{
  \item{x}{
An object of class \code{LinkedPairs}.
}
  \item{quote}{
Logical indicating whether to print the output surrounded by quotes.
}
  \item{right}{
Logical specifying whether to right align strings.
}
  \item{\dots}{
Other arguments for \code{print}.
}
}
\value{
An object of class ''LinkedPairs''.
}
\details{
Objects of class \code{LinkedPairs} are stored as square matrices of list elements with \code{dimnames} derived from the \code{dimnames} of the object of class ''Synteny'' from which it was created. The diagonal of the matrix is only filled if \code{OutputFormat} ''Comprehensive'' is selected in \code{NucleotideOverlap}, in which case it will be filled with the gene locations supplied to \code{GeneCalls}. The upper triangle is always filled, and contains location information in nucleotide space for all syntenic hits that link features between sequences in the form of an integer matrix with named columns. ''QueryGene'' and ''SubjectGene'' correspond to the integer rownames of the supplied gene calls. ''QueryIndex'' and ''SubjectIndex'' correspond to ''Index1'' and ''Index2'' columns of the source synteny object position. Remaining columns describe the exact positioning and size of extracted hits. The lower triangle is not filled if \code{OutputFormat} ''Sparse'' is selected and contains relative displacement positions for the 'left-most' and 'right-most' hit involved in linking the particular features indicated in the related line up the corresponding position in the upper triangle.

The object serves only as a simple package for input data to the \code{PairSummaries} function, and as such may not be entirely user friendly. However it has been left exposed to the user should they find this data interesting.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}