\name{TreeLine}
\alias{TreeLine}
\title{
Construct a Phylogenetic Tree
}
\description{
Builds a phylogenetic tree from a set of sequences or distance matrix.
}
\usage{
TreeLine(myXStringSet = NULL,
         myDistMatrix = NULL,
         method = "ML",
         type = "dendrogram",
         model = MODELS,
         cutoff = -Inf,
         showPlot = FALSE,
         collapse = -1,
         reconstruct = FALSE,
         root = 0,
         informationCriterion = "AICc",
         maxGenerations = 20,
         maxTime = Inf,
         quadrature = FALSE,
         costMatrix = 1 - diag(length(x)),
         processors = 1,
         verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
A \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} if \code{method} is \code{"ML"}, \code{method} is \code{"MP"}, or \code{reconstruct} is not \code{FALSE}.
}
  \item{myDistMatrix}{
A symmetric \eqn{N} x \eqn{N} distance matrix with the values of dissimilarity between \eqn{N} sequences or an object of class \code{'dist'}.  If the \code{method} is \code{"ML"} or \code{"MP"} then \code{myDistMatrix} is optional.
}
  \item{method}{
The phylogenetic method to be used.  This should be (an abbreviation of) one of \code{"ML"} (the default), \code{"MP"}, \code{"NJ"}, \code{"complete"}, \code{"single"}, \code{"UPGMA"}, or \code{"WPGMA"}.  (See details section below.)
}
  \item{type}{
Character string indicating the type of output desired.  This should be (an abbreviation of) one of \code{"dendrogram"} (the default), \code{"clusters"}, or \code{"both"}.  (See value section below.)
}
  \item{model}{
One or more of the available \code{MODELS} of evolution.  Only applicable if \code{method} is \code{"ML"} or \code{reconstruct} is not \code{FALSE}.  Automatic model selection will be performed if more than one \code{model} is provided.
}
  \item{cutoff}{
A vector with the maximum edge length separating the sequences in the same cluster.  A negative value (the default) will prevent clustering.  Multiple cutoffs may be provided in ascending or descending order.  (See details section below.)
}
  \item{showPlot}{
Logical specifying whether or not to plot the resulting dendrogram.
}
  \item{collapse}{
Numeric controlling which internal edges of the tree are removed by collapsing their nodes.  If \code{collapse} is zero then nodes at the same height will be collapsed to a single node, resulting in a multifurcating tree.  When \code{collapse} is greater than zero, nodes that are within \code{collapse} difference in height are made into a single node.  A value of \code{collapse} less than zero (the default) will ensure that the \code{dendrogram} is purely bifurcating.  Note that \code{collapse} has no effect on cluster numbers or \code{cutoff}.
}
  \item{reconstruct}{
Logical or numeric determining whether to perform ancestral state reconstruction when \code{myXStringSet} is specified. If \code{TRUE}, maximum likelihood character states are determined at internal nodes of the \code{dendrogram} and provided as the \code{"state"} attribute, unless \code{method} is \code{"MP"}, in which case the most parsimonious state is given at each node.  A numeric value between zero and one (exclusive) can be provided when \code{method} is not \code{MP}, in which case that fraction of the (marginal) likelihood for a state must be greater than the likelihood of all alternative states, otherwise a more ambiguous degeneracy code is used.  Only applicable if \code{type} is \code{"dendrogram"} (the default) or \code{"both"}.
}
  \item{root}{
Integer specifying the index of the outgroup or \code{0} (the default) to midpoint root the dendrogram.
}
  \item{informationCriterion}{
Character string specifying which information criterion to use in automatic model selection.  Must be (an abbreviation of) either \code{"AICc"} or \code{"BIC"}.  The best \code{model} is automatically chosen based on the \code{informationCriterion} calculated from the likelihood and the sample size (defined as the number of variable sites in \code{myXStringSet}).  Only applicable if \code{method} is \code{"ML"}.
}
  \item{maxGenerations}{
Integer indicating the maximum number of generations of optimization to perform.  More generations will potentially better search tree space at the expense of added runtime.  Only applicable if \code{method} is \code{"ML"} or \code{"MP"}.
}
  \item{maxTime}{
Numeric giving the maximum number of hours the algorithm is allowed to run before returning a result.  Once \code{maxTime} is reached, the algorithm will proceed at the next available opportunity.  Therefore, the \code{maxTime} is a soft maximum and does not guarantee results will be returned at exactly the allotted time.  Only applicable if \code{method} is \code{"ML"} or \code{"MP"}.
}
  \item{quadrature}{
Logical determining whether to use the Laguerre quadrature or equal-sized bins when discretizing the rate distribution across sites.  The default is to use equal-sized bins for direct comparison among likelihoods computed by other programs, although the Laguerre quadrature theoretically offers an improvement in accuracy for the same number of rates (Felsenstein, 2001).
}
  \item{costMatrix}{
A function or call written as a function of `x', which will evaluate to a square numeric matrix with one row for each state.  The cost matrix sets the penalties used in Sankoff parsimony.  The default is to apply a cost of \code{1} for state changes and \code{0} otherwise (i.e., equivalent to Fitch parsimony).  Only applicable if \code{method} is \code{"MP"}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.  Note, the number of processors in some steps is automatically selected between one and \code{processors} to optimize performance.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{TreeLine} builds a phylogenetic tree using either \code{myDistMatrix} and/or \code{myXStringSet}.  The output is either a \code{dendrogram} and/or \code{data.frame} containing cluster numbers.

Typical use cases include:
(1) Provide an \code{XStringSet} and obtain a maximum likelihood or maximum parsimony tree, optionally with ancestral states.
(2) Provide a \code{matrix} or \code{dist} object (e.g., from \code{DistanceMatrix}) and obtain a neighbor-joining or ultrametric tree, optionally with clusters of similarity defined by distance \code{cutoff}s.

Multiple \code{method}s of tree building are supported:
(1) Maximum likelihood: \code{ML} iteratively maximizes the likelihood of the tree and \code{model} parameters given aligned sequences (\code{myXStringSet}).  One or more \code{MODELS} of sequence evolution must be specified, of which the best \code{model} is automatically selected.
(2) Maximum parsimony: \code{MP} iterative maximizes the parsimony of a tree given aligned sequences (\code{myXStringSet}) and a \code{costMatrix}.  The default cost matrix is binary, corresponding to Fitch (1971) parsimony.
(3) Neighbor-joining: \code{NJ} uses the Neighbor-Joining method proposed by Saitou and Nei, which creates a reasonable minimum evolution tree from a distance matrix (\code{myDistMatrix}).
(4) Ultrametric: The method \code{complete} assigns clusters using complete-linkage so that sequences in the same cluster are no more than \code{cutoff} distance apart.  The method \code{single} assigns clusters using single-linkage so that sequences in the same cluster are within \code{cutoff} of at least one other sequence in the same cluster.  \code{UPGMA} and \code{WPGMA} assign clusters using average-linkage which is a compromise between the sensitivity of complete-linkage clustering to outliers and the tendency of single-linkage clustering to connect distant relatives that do not appear to be closely related.  \code{UPGMA} produces an unweighted tree, where each leaf contributes equally to the average edge lengths, whereas \code{WPGMA} produces a weighted result.

The returned dendrogram has information stored in its attributes, which can be accessed with the \code{attributes} or \code{attr} functions.  When \code{method} is \code{"ML"} or \code{"MP"}, edges of the tree have branch support values included in their attributes.  These \code{"support"} values represent the fraction of computed trees with that partition separating the sequences to either side, which is a proxy for confidence.  For maximum likelihood trees, the edges also have a \code{"probability"} representing the aBayes support probability (Anisimova, 2011).  If \code{reconstruct} is not \code{FALSE}, each edge of the tree will have a \code{"state"} attribute representing the (marginal) likelihood at that node or most parsimonious state when \code{method} is \code{"MP"}.  Maximum likelihood trees will also provide the likelihoods at each site when \code{reconstruct} is not \code{FALSE}.

When non-negative \code{cutoff}(s) are supplied, \code{TreeLine} will assign clusters based on edge lengths in the tree.  Multiple cutoffs may be provided in sorted order.  If the \code{cutoff}s are provided in \emph{descending} order then clustering at each new value of \code{cutoff} is continued within the prior \code{cutoff}'s clusters.  In this way clusters at lower values of \code{cutoff} are completely contained within their umbrella clusters at higher values of \code{cutoff}.  This is useful for defining taxonomy, where groups need to be hierarchically nested.  If multiple cutoffs are provided in \emph{ascending} order then clustering at each level of \code{cutoff} is independent of the prior level.
}
\value{
If \code{type} is \code{"dendrogram"} (the default), then a tree of class \code{dendrogram} is returned with attributes including pertinent information.
If \code{type} is \code{"clusters"}, then a \code{data.frame} is returned with a column for each cutoff specified.  This data.frame has dimensions \eqn{N*M}, where each one of \eqn{N} sequences is assigned to a cluster at the \eqn{M}-level of cutoff.
If \code{type} is \code{"both"} then a list is returned containing both the \code{"clusters"} and \code{"dendrogram"} outputs.
}
\note{
Note that the cophenetic distance between leaves of the \code{dendrogram} is equal to the sum of branch lengths separating the leaves.  This is the typical phylogenetic interpretation but different than trees produced by \code{hclust} where leaves are merged at a height equal to their cophenetic distance.  Hence, always use \code{Cophenetic} (rather than \code{cophenetic}) to compute correct cophenetic distances.
}
\references{
Anisimova M., Gil, M., Dufayard, J., Dessimoz, C., Gascuel, O. (2011) Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. \emph{Syst Biol.}, \bold{60(5)},685-99.

Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maximum likelihood approach. \emph{Journal of Molecular Evolution}, \bold{17(6)}, 368-376.

Felsenstein J. (2001) Taking variation of evolutionary rates between sites into account in inferring phylogenies. \emph{Journal of molecular evolution}, \bold{53}:447-455.

Fitch, W. M. (1971) Toward defining the course of evolution: minimum change for a specified tree topology. \emph{Systematic Zoology}, \bold{20}:406-416.

Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees. \emph{Molecular Biology and Evolution}, \bold{4(4)}, 406-425.

Sankoff, D. (1975) Minimal mutation trees of sequences. \emph{SIAM Journal of Applied Math}, \bold{28}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{Cophenetic}}, \code{\link{DistanceMatrix}}, \code{\link{MapCharacters}}, \code{\link{MODELS}}, \code{\link{WriteDendrogram}}
}
\examples{
# using the matrix from the original paper by Saitou and Nei
m <- matrix(0,8,8) # only the lower triangle is used
m[2:8,1] <- c(7, 8, 11, 13, 16, 13, 17)
m[3:8,2] <- c(5, 8, 10, 13, 10, 14)
m[4:8,3] <- c(5, 7, 10, 7, 11)
m[5:8,4] <- c(8, 11, 8, 12)
m[6:8,5] <- c(5, 6, 10)
m[7:8,6] <- c(9, 13)
m[8,7] <- c(8)

# returns an object of class "dendrogram"
tree <- TreeLine(myDistMatrix=m, cutoff=10, method="NJ", showPlot=TRUE)

# example of specifying multiple cutoffs
clusters <- TreeLine(myDistMatrix=m, method="UPGMA", type="clusters", cutoff=c(2,6,10,20))
head(clusters)

# example of creating a complete-linkage tree from an alignment
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
dna <- AlignTranslation(dna, verbose=FALSE)
dna # input alignment
d <- DistanceMatrix(dna, type="dist") # returns an object of class 'dist'
complete <- TreeLine(myDistMatrix=d, method="complete", cutoff=0.05, showPlot=TRUE)

\dontrun{
# example of maximum likelihood tree building
tree <- TreeLine(myXStringSet=head(dna, 10), reconstruct=TRUE)
plot(tree)

# displaying ancestral states on each edge
start <- 50
end <- 52
tree <- dendrapply(tree,
	function(x) {
		attr(x, "edgetext") <- substring(attr(x, "state"),
			start,
			end)
		x
	})
plot(tree,
	edgePar=list(p.lwd=1e-5,
		p.col="#00000000",
		t.col=c("#AA3355", "#33FFFF")),
	edge.root=FALSE)
}

# example of maximum parsimony tree building
costs <- matrix(c(0, 2, 1, 2, 2, 0, 2, 1, 1, 2, 0, 2, 2, 1, 2, 0), 4)
dimnames(costs) <- list(DNA_BASES, DNA_BASES)
costs
tree <- TreeLine(myXStringSet=head(dna, 10), method="MP", costMatrix=costs[x, x])

# example of accessing and using the attributes
attributes(tree)
plot(dendrapply(tree,
	function(x) {
		s <- attr(x, "support")
		if (!is.null(s))
			attr(x, "edgetext") <- formatC(as.numeric(s), digits=2, format="f")
		attr(x, "edgePar") <- list(p.col=NA, p.lwd=1e-5, t.col="#CC55AA", t.cex=0.7)
		x
	}),
	horiz=TRUE)
}