\name{AlignSeqs}
\alias{AlignSeqs}
\title{
Align a Set of Unaligned Sequences
}
\description{
Performs profile-to-profile alignment of multiple unaligned sequences following a guide tree.
}
\usage{
AlignSeqs(myXStringSet,
         guideTree = NULL,
         iterations = 2,
         refinements = 1,
         gapOpening = c(-18, -14),
         gapExtension = c(-3, -2),
         useStructures = TRUE,
         structures = NULL,
         FUN = AdjustAlignment,
         levels = c(0.9, 0.7, 0.7, 0.4, 10, 5, 5, 2),
         alphabet = AA_REDUCED[[1]],
         processors = 1,
         verbose = TRUE,
         \dots)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of unaligned sequences.
}
  \item{guideTree}{
Either \code{NULL} or a \code{dendrogram} giving the ordered tree structure in which to align profiles.  If \code{NULL} then a guide tree will be automatically constructed based on the order of shared k-mers.
}
  \item{iterations}{
Number of iteration steps to perform.  During each iteration step the guide tree is regenerated based on the alignment and the sequences are realigned.
}
  \item{refinements}{
Number of refinement steps to perform.  During each refinement step groups of sequences are realigned to rest of the sequences, and the best of these two alignments (before and after realignment) is kept.
}
  \item{gapOpening}{
Single numeric giving the cost for opening a gap in the alignment, or two numbers giving the minimum and maximum costs.  In the latter case the cost will be varied depending upon whether the groups of sequences being aligned are nearly identical or maximally distant.
}
  \item{gapExtension}{
Single numeric giving the cost for extending an open gap in the alignment, or two numbers giving the minimum and maximum costs.  In the latter case the cost will be varied depending upon whether the groups of sequences being aligned are nearly identical or maximally distant.
}
  \item{useStructures}{
Logical indicating whether to use secondary structure predictions during alignment.  If \code{TRUE} (the default), secondary structure probabilities will be automatically calculated for amino acid and RNA sequences if they are not provided (i.e., when \code{structures} is \code{NULL}).
}
  \item{structures}{
Either a list of secondary structure probabilities matching the \code{structureMatrix}, such as that output by \code{PredictHEC} or \code{PredictDBN}, or \code{NULL} to generate the structures automatically.  Only applicable if \code{myXStringSet} is an \code{AAStringSet} or \code{RNAStringSet}.
}
  \item{FUN}{
A function to be applied after each profile-to-profile alignment.  (See details section below.)
}
  \item{levels}{
Numeric with eight elements specifying the levels at which to trigger events.  (See details section below.)
}
  \item{alphabet}{
Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups.  Alphabet reduction helps to find more distant homologies between sequences.  A non-reduced amino acid alphabet can be used by setting \code{alphabet} equal to \code{AA_STANDARD}.  Only applicable if \code{myXStringSet} is an \code{AAStringSet}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{AlignProfiles}}, including \code{perfectMatch}, \code{misMatch}, \code{gapPower}, \code{terminalGap}, \code{restrict}, \code{anchor}, \code{normPower}, \code{standardize}, \code{substitutionMatrix}, and \code{structureMatrix}.
}
}
\details{
The profile-to-profile method aligns a sequence set by merging profiles along a guide tree until all the input sequences are aligned.  This process has three main steps:  (1)  If \code{guideTree=NULL}, an initial single-linkage guide tree is constructed based on a distance matrix of shared k-mers.  Alternatively, a \code{dendrogram} can be provided as the initial \code{guideTree}.  (2)  If \code{iterations} is greater than zero, then a UPGMA guide tree is built based on the initial alignment and the sequences are re-aligned along this tree.  This process repeated \code{iterations} times or until convergence.  (3)  If \code{refinements} is greater than zero, then subsets of the alignment are re-aligned to the remainder of the alignment.  This process generates two alignments, the best of which is chosen based on its sum-of-pairs score.  This refinement process is repeated \code{refinements} times, or until convergence.

The purpose of \code{levels} is to speed-up the alignment process by not running time consuming processes when they are unlikely to change the outcome.  The first four \code{levels} control when \code{refinements} occur and the function \code{FUN} is run on the alignment.  The default \code{levels} specify that these events should happen when above 0.9 (AA; \code{levels[1]}) or 0.7 (DNA/RNA; \code{levels[3]}) average dissimilarity on the initial tree, when above 0.7 (AA; \code{levels[2]}) or 0.4 (DNA/RNA; \code{levels[4]}) average dissimilarity on the iterative tree(s), and after every tenth improvement made during refinement.  The sixth element of levels (\code{levels[6]}) prevents \code{FUN} from being applied at any point to less than 5 sequences.

The \code{FUN} function is always applied before returning the alignment so long as there are at least \code{levels[6]} sequences.  The default \code{FUN} is \code{AdjustAlignment}, but \code{FUN} can be any function that takes in an \code{XStringSet} as its first argument, as well as \code{weights}, \code{processors}, and \code{substitutionMatrix} as optional arguments.  For example, the default \code{FUN} could be altered to not perform any changes by setting it equal to \code{function(x, ...) return(x)}, where \code{x} is an \code{XStringSet}.

Secondary structures are automatically computed for amino acid and RNA sequences unless \code{structures} are provided or \code{useStructures} is \code{FALSE}.  The default \code{structureMatrix} is used unless an alternative is provided.  For RNA sequences, secondary structures are only computed when the total length of the initial guide tree is at least 5 (\code{levels[7]}) or the length of subsequent trees is at least 2 (\code{levels[8]}).
}
\value{
An \code{XStringSet} of aligned sequences.
}
\references{
Wright, E. S. (2015). DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. BMC Bioinformatics, 16, 322. http://doi.org/10.1186/s12859-015-0749-z

Wright, E. S. (2020). RNAconTest: comparing tools for noncoding RNA multiple sequence alignment based on structural consistency. RNA 2020, 26, 531-540.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AdjustAlignment}}, \code{\link{AlignDB}}, \code{\link{AlignProfiles}}, \code{\link{AlignSynteny}}, \code{\link{AlignTranslation}}, \code{\link{ReadDendrogram}}, \code{\link{TreeLine}}, \code{\link{StaggerAlignment}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db, remove="all")
alignedDNA <- AlignSeqs(dna)
BrowseSeqs(alignedDNA, highlight=1)

# use secondary structure with RNA sequences
alignedRNA <- AlignSeqs(RNAStringSet(dna))
BrowseSeqs(alignedRNA, highlight=1)
}
