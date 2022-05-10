\name{IdTaxa}
\alias{IdTaxa}
\title{
Assign Sequences a Taxonomic Classification
}
\description{
Classifies sequences according to a training set by assigning a confidence to taxonomic labels for each taxonomic level.
}
\usage{
IdTaxa(test,
       trainingSet,
       type = "extended",
       strand = "both",
       threshold = 60,
       bootstraps = 100,
       samples = L^0.47,
       minDescend = 0.98,
       fullLength = 0,
       processors = 1,
       verbose = TRUE)
}
\arguments{
  \item{test}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} of unaligned sequences.
}
  \item{trainingSet}{
An object of class \code{Taxa} and subclass {Train}.
}
  \item{type}{
Character string indicating the type of output desired.  This should be (an abbreviation of) one of \code{"extended"} or \code{"collapsed"}.  (See value section below.)
}
  \item{strand}{
Character string indicating the orientation of the \code{test} sequences relative to the \code{trainingSet}.  This should be (an abbreviation of) one of \code{"both"}, \code{"top"}, or \code{"bottom"}.  The top strand is defined as the input \code{test} sequences being in the same orientation as the \code{trainingSet}, and the bottom strand is its reverse complement orientation.  The default of \code{"both"} will classify using both orientations and choose the result with highest confidence.  Choosing the correct \code{strand} will make classification over 2-fold faster, assuming that all of the reads are in the same orientation. Note that \code{strand} is ignored when \code{test} is an \code{AAStringSet}.
}
  \item{threshold}{
Numeric specifying the confidence at which to truncate the output taxonomic classifications.  Lower values of \code{threshold} will classify deeper into the taxonomic tree at the expense of accuracy, and vise-versa for higher values of \code{threshold}.
}
  \item{bootstraps}{
Integer giving the maximum number of bootstrap replicates to perform for each sequence.  The number of bootstrap replicates is set automatically such that (on average) 99\% of k-mers are sampled in each \code{test} sequence.
}
  \item{samples}{
A function or call written as a function of `L', which will evaluate to a numeric vector the same length as `L'.  Typically of the form ``\code{A + B*L^C}'', where `A', `B', and `C' are constants.
}
  \item{minDescend}{
Numeric giving the minimum fraction of \code{bootstraps} required to descend the tree during the initial tree descend phase of the algorithm.  Higher values are less likely to descend the tree, causing direct comparison against more sequences in the \code{trainingSet}.  Lower values may increase classification speed at the expense of accuracy.  Suggested values are between \code{1.0} and \code{0.9}.
}
  \item{fullLength}{
Numeric specifying the fold-difference in sequence lengths between sequences in \code{test} and \code{trainingSet} that is allowable, or \code{0} (the default) to consider all sequences in \code{trainingSet} regardless of length.  Can be specified as either a single numeric (> 1), or two numerics specifying the upper and lower fold-difference.  If \code{fullLength} is between \code{0} and \code{1} (exclusive), the fold-difference is inferred from the length variability among sequences belonging to each class based on the \code{foldDifference} quantiles.  For example, setting \code{fullLength} to \code{0.99} would use the 1st and 99th percentile of intra-group length variability from the \code{trainingSet}.  In the case of full-length sequences, specifying \code{fullLength} can improve both speed and accuracy by using sequence length as a pre-filter to classification.  Note that \code{fullLength} should only be greater than \code{0} when both the \code{test} and \code{trainingSet} consist of full-length sequences.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Sequences in \code{test} are each assigned a taxonomic classification based on the \code{trainingSet} created with \code{\link{LearnTaxa}}.  Each taxonomic level is given a confidence between 0\% and 100\%, and the taxonomy is truncated where confidence drops below \code{threshold}.  If the taxonomic classification was truncated, the last group is labeled with ``unclassified_'' followed by the final taxon's name.  Note that the reported confidence is not a p-value but does directly relate to a given classification's probability of being wrong.  The default \code{threshold} of \code{60\%} is intended to minimize the rate of incorrect classifications.  Lower values of \code{threshold} (e.g., \code{50\%}) may be preferred to increase the taxonomic depth of classifications.  Values of \code{60\%} or \code{50\%} are recommended for nucleotide sequences and \code{50\%} or \code{40\%} for amino acid sequences.
}
\value{
If \code{type} is \code{"extended"} (the default) then an object of class \code{Taxa} and subclass {Train} is returned.  This is stored as a list with elements corresponding to their respective sequence in \code{test}.  Each list element contains components:
  \item{taxon}{
A character vector containing the taxa to which the sequence was assigned.
}
  \item{confidence}{
A numeric vector giving the corresponding percent confidence for each taxon.
}
  \item{rank}{
If the classifier was trained with a set of \code{rank}s, a character vector containing the rank name of each taxon.
}

If \code{type} is \code{"collapsed"} then a character vector is returned with the taxonomic assignment for each sequence.  This takes the repeating form ``Taxon name [rank, confidence\%]; ...'' if \code{rank}s were supplied during training, or ``Taxon name [confidence\%]; ...'' otherwise.
}
\references{
Murali, A., et al. (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. Microbiome, 6, 140. https://doi.org/10.1186/s40168-018-0521-5

Cooley, N. and Wright, E. (2021). Accurate annotation of protein coding sequences with IDTAXA. NAR Genomics and Bioinformatics, \bold{3(3)}. https://doi.org/10.1093/nargab/lqab080
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{LearnTaxa}}, \code{\link{Taxa-class}}
}
\examples{
data("TrainingSet_16S")

# import test sequences
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

# remove any gaps in the sequences
dna <- RemoveGaps(dna)

# classify the test sequences
ids <- IdTaxa(dna, TrainingSet_16S, strand="top")
ids

# view the results
plot(ids, TrainingSet_16S)
}