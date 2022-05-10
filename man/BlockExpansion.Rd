\name{BlockExpansion}
\alias{BlockExpansion}
\title{
Attempt to expand blocks of paired features in a \code{PairSummaries} object.
}
\description{
Attempt to expand blocks of paired features in a \code{PairSummaries} object.
}
\usage{
BlockExpansion(Pairs,
               GapTolerance = 4L,
               DropSingletons = FALSE,
               Criteria = "PID",
               Floor = 0.5,
               NewPairsOnly = TRUE,
               DBPATH,
               Verbose = FALSE)
}
\arguments{
  \item{Pairs}{
An object of class \code{PairSummaries}.
}
  \item{GapTolerance}{
Integer value indicating the \code{diff} between feature IDs that can be tolerated to view features as part of the same block. Set by default to \code{4L}, implying that a single feature missing in a run of pairs will not cause the block to be split. Setting to \code{3L} would imply that a \code{diff} of 3 between features, or a gap of 2 features, can be viewed as those features being part of the same block.
}
  \item{DropSingletons}{
Ignore solo pairs when planning expansion routes. Set to \code{FALSE} by default.
}
  \item{Criteria}{
Either ``PID'' or ``Score'', indicating which metric to use to keep or reject pairs.
}

  \item{Floor}{
Lower PID limit for keeping a pair that was evaluated during expansion.
}
  \item{NewPairsOnly}{
Logical indicating whether or not to return only the pairs that were kept from all expansion attempts, or to return a \code{PairSummaries} object with the new pairs folded in.
}
  \item{DBPATH}{
A file or connection pointing to the \code{DECIPHER} database supplied to \code{FindSynteny} for the original map construction.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{BlockExpansion} uses a naive expansion algorithm to attempt to fill in gaps in blocks of paired features and to attempt to expand blocks of paired features.
}
\value{
An object of class \code{PairSummaries}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{PairSummaries}}, \code{\link{NucleotideOverlap}}, \code{link{SubSetPairs}}, \code{\link{FindSynteny}}
}
\examples{
DBPATH <- system.file("extdata",
                      "Endosymbionts.sqlite",
                      package = "SynExtend")
                      
data("Endosymbionts_Pairs01", package = "SynExtend")
Pairs02 <- BlockExpansion(Pairs = Endosymbionts_Pairs01,
                          NewPairsOnly = FALSE,
                          DBPATH = DBPATH,
                          Verbose = TRUE)
}
