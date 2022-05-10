\name{DisjointSet}
\alias{DisjointSet}
\title{
Return single linkage clusters from \code{PairSummaries} objects.
}
\description{
Takes in a \code{PairSummaries} object and return a list of identifiers organized into single linkage clusters.
}
\usage{
DisjointSet(Pairs,
            Verbose = FALSE)
}
\arguments{
  \item{Pairs}{
A \code{PairSummaries} object.
}
  \item{Verbose}{
Logical indicating whether to print progress bars and messages. Defaults to \code{FALSE}.
}
}
\details{
Takes in a \code{PairSummaries} object and return a list of identifiers organized into single linkage clusters.
}
\value{
Returns a list of character vectors representing IDs of sequence features, typically genes.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{FindSynteny}}, \code{\link{Synteny-class}}, \code{\link{PairSummaries}}, \code{\link{FindSets}}
}
\examples{
data("Endosymbionts_Pairs03", package = "SynExtend")

Sets <- DisjointSet(Pairs = Endosymbionts_Pairs03,
                    Verbose = TRUE)
}
