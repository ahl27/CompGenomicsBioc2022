\name{Cophenetic}
\alias{Cophenetic}
\title{
Compute cophenetic distances on dendrogram objects
}
\description{
Calculates the matrix of cophenetic distances represented by a dendrogram object.
}
\usage{
Cophenetic(x)
}
\arguments{
  \item{x}{
A dendrogram object.
}
}
\details{
The cophenetic distance between two observations is defined as the branch length separating them on a dendrogram.  This function differs from the \code{cophenetic} function in that it does not assume the tree is ultrametric and outputs the branch length separating pairs of observations rather than the height of their merger. A dendrogram that better preserves a distance matrix will show higher correlation between the distance matrix and it cophenetic distances.
}
\value{
An object of class 'dist'.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{TreeLine}}
}
\examples{
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
d1 <- DistanceMatrix(dna, type="dist")
dend <- TreeLine(myDistMatrix=d1, method="NJ")
d2 <- Cophenetic(dend)
cor(d1, d2)
}