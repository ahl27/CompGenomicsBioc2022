\name{IdClusters}
\alias{IdClusters}
\title{
Cluster Sequences By Similarity
}
\description{
Groups the sequences into approximate clusters of similarity.
}
\usage{
IdClusters(myXStringSet,
           cutoff = 0,
           processors = 1,
           verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
If \code{method} is \code{"ML"} or \code{reconstruct} is not \code{FALSE}, the \code{DNAStringSet} or \code{RNAStringSet} used in the creation of \code{myDistMatrix}.  If \code{method} is \code{"inexact"}, the \code{DNAStringSet}, \code{RNAStringSet}, or \code{AAStringSet} to cluster.
}
  \item{cutoff}{
A vector with the maximum edge length separating the sequences in the same cluster.  Multiple cutoffs may be provided in ascending or descending order.  (See details section below.)
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{IdClusters} groups the input sequences into approximate clusters using a heuristic algorithm.  First the sequences are ordered by length and the longest sequence becomes the first cluster seed.  If the second sequence is less than \code{cutoff} k-mer distance then it is added to the cluster, otherwise it becomes a new cluster representative.  The remaining sequences are matched to cluster representatives in a similar fashion until all sequences belong to a cluster.  In the majority of cases, this process results in clusters with members separated by less than \code{cutoff} distance.

Multiple cutoffs may be provided in sorted order.  If the \code{cutoff}s are provided in \emph{descending} order then clustering at each new value of \code{cutoff} is continued within the prior \code{cutoff}'s clusters.  In this way clusters at lower values of \code{cutoff} are completely contained within their umbrella clusters at higher values of \code{cutoff}.  This is useful for accelerating the clustering process, because each subsequent group is only clustered within the previous group.  If multiple cutoffs are provided in \emph{ascending} order then clustering at each level of \code{cutoff} is independent of the prior level.
}
\value{
A data.frame is returned with dimensions \eqn{N*M}, where each one of \eqn{N} sequences is assigned to a cluster at the \eqn{M}-level of cutoff.  The row.names of the data.frame correspond to the \emph{names} of \code{myXStingSet}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\examples{
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
aa <- translate(dna)
inexact <- IdClusters(aa, cutoff=seq(0.7, 0.1, -0.1))
head(inexact)
apply(inexact, 2, max) # number of clusters per cutoff
}