\name{DistanceMatrix}
\alias{DistanceMatrix}
\title{
Calculate the Distances Between Sequences
}
\description{
Calculates a distance matrix for an \code{XStringSet}.  Each element of the distance matrix corresponds to the dissimilarity between two sequences in the \code{XStringSet}.
}
\usage{
DistanceMatrix(myXStringSet,
               type = "matrix",
               includeTerminalGaps = FALSE,
               penalizeGapLetterMatches = TRUE,
               penalizeGapGapMatches = FALSE,
               correction = "none",
               processors = 1,
               verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet}, \code{RNAStringSet}, or \code{AAStringSet} object of aligned sequences.
}
  \item{type}{
Character string indicating the type of output desired.  This should be either \code{"matrix"} or \code{"dist"}.  (See value section below.)
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" or "." characters on each end of the sequence) into the calculation of distance.
}
  \item{penalizeGapLetterMatches}{
Logical specifying whether or not to consider gap-to-letter matches as mismatches.  If \code{FALSE}, then gap-to-letter matches are not included in the total length used to calculate distance.
}
  \item{penalizeGapGapMatches}{
Logical specifying whether or not to consider gap-to-gap matches as mismatches.  If \code{FALSE} (the default), then gap-to-gap matches are not included in the total length used to calculate distance.
}
  \item{correction}{
The substitution model used for distance correction.  This should be (an abbreviation of) either \code{"none"}, \code{"Jukes-Cantor"} (i.e., \code{"JC69"}), or \code{"F81"}.  For \code{"F81"} letter frequencies are derived from \code{myXStringSet}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The uncorrected (\code{correction = "none"}) distance matrix represents the hamming distance between each of the sequences in \code{myXStringSet}.  Ambiguity can be represented using the characters of the \code{IUPAC_CODE_MAP} for \code{DNAStringSet} and \code{RNAStringSet} inputs, or using the \code{AMINO_ACID_CODE} for an \code{AAStringSet} input.  For example, the distance between an 'N' and any other nucleotide base is zero.  The letters B (N or D), J (I or L), Z (Q or E), and X (any letter) are degenerate in the \code{AMINO_ACID_CODE}.

If \code{includeTerminalGaps = FALSE} then terminal gaps ("-" or "." characters) are not included in sequence length.  This can be faster since only the positions common to each pair of sequences are compared.  Sequences with no overlapping region in the alignment are given a value of \code{NA}, unless \code{includeTerminalGaps = TRUE}, in which case distance is 100\%.  Masked characters (\code{"+"}) in either sequence are not considered in distance.

Penalizing gap-to-gap and gap-to-letter mismatches specifies whether to penalize these special mismatch types and include them in the total length when calculating distance.  Both "-" and "." characters are interpreted as gaps.  The default behavior is to calculate distance as the fraction of positions that differ across the region of the alignment shared by both sequences (not including gap-to-gap matches).

Two correction factors are available, \code{"JC69"} and \code{"F81"}, which are described in \code{\link{MODELS}}.  Both transform raw distance (\eqn{d}) by \deqn{-E * \log \left( 1 - d / E \right)}, where \deqn{E = 1 - \sum_{i \in sym} f_i^2} and (\eqn{f}) is the relative frequency of each symbol (\eqn{sym}).  In the \code{"JC69"} model symbols are assumed to have equal frequency, whereas in the \code{"F81"} model the symbol frequencies are empirically derived from the input \code{myXStringSet}.  Note that gaps are treated as an additional symbol when \code{penalizeGapLetterMatches} or \code{penalizeGapGapMatches} is \code{TRUE}.

The elements of the distance matrix can be referenced by \code{dimnames} corresponding to the \code{names} of the \code{XStringSet}.  Additionally, an attribute named "correction" specifying the method of correction used can be accessed using the function \code{attr}.
}
\value{
If \code{type} is \code{"matrix"}, a symmetric matrix where each element is the distance between the sequences referenced by the respective row and column.  The \code{dimnames} of the matrix correspond to the \code{names} of the \code{XStringSet}.

If \code{type} is \code{"dist"}, an object of \code{class} \code{"dist"} that contains one triangle of the distance matrix as a vector.  Since the distance matrix is symmetric, storing only one triangle is more memory efficient.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{TreeLine}}
}
\examples{
# example of using the defaults:
dna <- DNAStringSet(c("ACTG", "ACCG"))
dna
DistanceMatrix(dna)

# changing the output type to "dist":
d <- DistanceMatrix(dna, type="dist")
d
length(d) # minimal memory space required
m <- as.matrix(d)
length(m) # more memory space required

# supplying an AAStringSet
aa <- AAStringSet(c("ASYK", "ATYK", "CTWN"))
aa
DistanceMatrix(aa)

# defaults compare intersection of internal ranges:
dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
dna
d <- DistanceMatrix(dna)
d
# d[1,2] is 1 base in 4 = 0.25

# compare the entire sequence, including gaps:
dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
dna
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE,
                    penalizeGapGapMatches=TRUE)
d
# d[1,2] is now 3 bases in 6 = 0.50

# compare union of internal positions, without terminal gaps:
dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
dna
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE,
                    penalizeGapGapMatches=FALSE)
d
# d[1,2] is now 2 bases in 5 = 0.40

# gap ("-") and unknown (".") characters are interchangeable:
dna <- DNAStringSet(c("ANGCT.", ".ACCT-"))
dna
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE,
                    penalizeGapGapMatches=FALSE)
d
# d[1,2] is still 2 bases in 5 = 0.40
}