\name{FindSynteny}
\alias{FindSynteny}
\title{
Finds Synteny in a Sequence Database
}
\description{
Finds syntenic blocks between groups of sequences in a database.
}
\usage{
FindSynteny(dbFile,
            tblName = "Seqs",
            identifier = "",
            useFrames = TRUE,
            alphabet = AA_REDUCED[[1]],
            geneticCode = GENETIC_CODE,
            sepCost = 0,
            gapCost = -0.01,
            shiftCost = 0,
            codingCost = 0,
            maxSep = 2000,
            maxGap = 5000,
            minScore = 30,
            dropScore = -100,
            maskRepeats = TRUE,
            allowOverlap = TRUE,
            storage = 0.5,
            processors = 1,
            verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table where the sequences are located.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If \code{""} then all identifiers are selected.  Repeated identifiers will find synteny between a sequence and itself, while blocking identical positions from matching in both sequences.
}
  \item{useFrames}{
Logical specifying whether to use 6-frame amino acid translations to help find more distant hits.  Using the \code{alphabet} is helpful when the genome is largely composed of coding DNA.  If \code{FALSE} then faster but less sensitive to distant homology.
}
  \item{alphabet}{
Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups.  Alphabet reduction helps to find more distant homologies between sequences.  A non-reduced amino acid alphabet can be used by setting \code{alphabet} equal to \code{AA_STANDARD}.
}
  \item{geneticCode}{
Either a character vector giving the genetic code to use in translation, or a list containing one genetic code for each identifier.  If a list is provided then it must be named by the corresponding identifiers in the database.
}
  \item{sepCost}{
Cost per nucleotide separation between hits to apply when chaining hits into blocks.
}
  \item{gapCost}{
Cost for gaps between hits to apply when chaining hits into blocks.
}
  \item{shiftCost}{
Cost for shifting between different reading frames when chaining reduced amino acid hits into blocks.
}
  \item{codingCost}{
Cost for switching between coding and non-coding hits when chaining hits into blocks.
}
  \item{maxSep}{
Maximal separation (in nucleotides) between hits in the same block.
}
  \item{maxGap}{
The maximum number of gaps between hits in the same block.
}
  \item{minScore}{
The minimum score required for a chain of hits to become a block.  Higher values of \code{minScore} are less likely to yield false positives.
}
  \item{dropScore}{
The change from maximal score required to stop extending blocks.
}
  \item{maskRepeats}{
Logical specifying whether to ``soft'' mask repeats when searching for hits.
}
  \item{allowOverlap}{
Logical specifying whether to permit blocks to overlap on the same sequence.
}
  \item{storage}{
Excess gigabytes available to store objects so that they do not need to be recomputed in later steps.  This should be a number between zero and a (modest) fraction of the available system memory.  Note that more than \code{storage} gigabytes may be required, but will not be stored for later reuse.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Long nucleotide sequences, such as genomes, are often not collinear or may be composed of many smaller segments (e.g., contigs).  \code{FindSynteny} searches for ``hits'' between sequences that can be chained into collinear ``blocks'' of synteny.  Hits are defined as k-mer exact nucleotide matches or k-mer matches in a reduced amino acid alphabet (if \code{useFrames} is \code{TRUE}).  Hits are chained into blocks as long as they are:  (1) within the same sequence, (2) within \code{maxSep} and \code{maxGap} distance, and (3) help maintain the score above \code{minScore}.  Blocks are extended from their first and last hit until their score drops below \code{dropScore} from the maximum that was reached.  This process results in a set of hits and blocks stored in an object of class ``Synteny''.
}
\value{
An object of class ``Synteny''.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\note{
\code{FindSynteny} is intended to be used on sets of sequences with up to ~100 million nucleotides total per identifier.  For this reason, better performance can sometimes be achieved by assigning a unique identifier to each chromosome belonging to a large genome.
}
\seealso{
\code{\link{AlignSynteny}}, \code{\link{Synteny-class}}
}
\examples{
db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
synteny <- FindSynteny(db)
synteny
pairs(synteny) # scatterplot matrix
}