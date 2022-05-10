\name{FindSets}
\alias{FindSets}
\title{
Find all single linkage clusters in an undirected pairs list.
}
\description{
Take in a pair of vectors representing the columns of an undirected pairs list and return the single linkage clusters.
}
\usage{
FindSets(p1,
         p2,
         Verbose = FALSE)
}
\arguments{
  \item{p1}{
Column 1 of a pairs matrix or list.
}
  \item{p2}{
Column 2 of a pairs matrix or list.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{FindSets} uses a version of the union-find algorithm to collect single linkage clusters from a pairs list. Currently meant to be used inside a wrapper function, but left exposed for user convenience.
}
\value{
A two column matrix with the first column being input nodes, and the second the node representing a single linkage cluster.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{PairSummaries}}
}
\examples{
set.seed(1986)
m <- cbind(as.integer(sample(30, size = 25,
                             replace = TRUE)),
           as.integer(sample(35, size = 25,
                             replace = TRUE)))

Levs <- unique(c(m[, 1],
                 m[, 2]))
m <- cbind("1" = as.integer(factor(x = m[, 1L],
                                   levels = Levs)),
           "2" = as.integer(factor(x = m[, 2L],
                                   levels = Levs)))
z <- FindSets(p1 = m[, 1],
              p2 = m[, 2])
}
