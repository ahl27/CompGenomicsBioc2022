\name{AlignTranslation}
\alias{AlignTranslation}
\title{
Align Sequences By Their Amino Acid Translation
}
\description{
Performs alignment of a set of DNA or RNA sequences by aligning their corresponding amino acid sequences.
}
\usage{
AlignTranslation(myXStringSet,
                 sense = "+",
                 direction = "5' to 3'",
                 readingFrame = NA,
                 type = class(myXStringSet),
                 geneticCode = GENETIC_CODE,
                 \dots)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} object of unaligned sequences.
}
  \item{sense}{
Single character specifying sense of the input sequences, either the positive (\code{"+"}) coding strand or negative (\code{"-"}) non-coding strand.
}
  \item{direction}{
Direction of the input sequences, either \code{"5' to 3'"} or \code{"3' to 5'"}.
}
  \item{readingFrame}{
Numeric vector giving a single reading frame for all of the sequences, or an individual reading frame for each sequence in \code{myXStringSet}.  The \code{readingFrame} can be either \code{1}, \code{2}, \code{3} to begin translating on the first, second, and third nucleotide position, or \code{NA} (the default) to guess the reading frame.  (See details section below.)
}
  \item{type}{
Character string indicating the type of output desired.  This should be (an abbreviation of) one of \code{"DNAStringSet"}, \code{"RNAStringSet"}, \code{"AAStringSet"}, or \code{"both"}.  (See value section below.)
}
  \item{geneticCode}{
Either a character vector giving the genetic code in the same format as \code{GENETIC_CODE} (the default), or a list containing one genetic code for each sequence in \code{myXStringSet}.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{AlignSeqs}}, including \code{gapOpening}, \code{gapExtension}, \code{gapPower}, \code{terminalGap}, \code{restrict}, \code{anchor}, \code{normPower}, \code{substitutionMatrix}, \code{structureMatrix}, \code{standardize}, \code{alphabet}, \code{guideTree}, \code{iterations}, \code{refinements}, \code{useStructures}, \code{structures}, \code{FUN}, and \code{levels}.
}
}
\details{
Alignment of proteins is often more accurate than alignment of their coding nucleic acid sequences.  This function aligns the input nucleic acid sequences via aligning their translated amino acid sequences.  First, the input sequences are translated according to the specified \code{sense}, \code{direction}, and \code{readingFrame}.  The resulting amino acid sequences are aligned using \code{AlignSeqs}, and this alignment is (conceptually) reverse translated into the original sequence type, \code{sense}, and \code{direction}.  Not only is alignment of protein sequences generally more accurate, but aligning translations will ensure that the reading frame is maintained in the nucleotide sequences.

If the \code{readingFrame} is \code{NA} (the default) then an attempt is made to guess the reading frame of each sequence based on the number of stop codons in the translated amino acids.  For each sequence, the first reading frame will be chosen (either \code{1}, \code{2}, or \code{3}) without stop codons, except in the final position.  If the number of stop codons is inconclusive for a sequence then the reading frame will default to \code{1}.  The entire length of each sequence is translated in spite of any stop codons identified.  Note that this method is only constructive in circumstances where there is a substantially long coding sequence with at most a single stop codon expected in the final position, and therefore it is preferable to specify the reading frame of each sequence if it is known.
}
\value{
An \code{XStringSet} of the class specified by \code{type}, or a list of two components (nucleotides and amino acids) if \code{type} is \code{"both"}.  Note that incomplete starting and ending codons will be translated into the mask character ("+") if the result includes an \code{AAStringSet}.
}
\references{
Wright, E. S. (2015). DECIPHER: harnessing local sequence context to improve protein multiple sequence alignment. BMC Bioinformatics, 16, 322. http://doi.org/10.1186/s12859-015-0749-z
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AlignDB}}, \code{\link{AlignProfiles}}, \code{\link{AlignSeqs}}, \code{\link{AlignSynteny}}, \code{\link{CorrectFrameshifts}}
}
\examples{
# first three sequences translate to  MFITP*
# and the last sequence translates as MF-TP*
rna <- RNAStringSet(c("AUGUUCAUCACCCCCUAA", "AUGUUCAUAACUCCUUGA",
	"AUGUUCAUUACACCGUAG", "AUGUUUACCCCAUAA"))
RNA <- AlignSeqs(rna, verbose=FALSE)
RNA

RNA <- AlignTranslation(rna, verbose=FALSE)
RNA

AA <- AlignTranslation(rna, type="AAStringSet", verbose=FALSE)
AA

BOTH <- AlignTranslation(rna, type="both", verbose=FALSE)
BOTH

# example of aligning many protein coding sequences:
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
DNA <- AlignTranslation(dna) # align the translation then reverse translate
DNA

# using a mixture of standard and non-standard genetic codes
gC1 <- getGeneticCode(id_or_name2="1", full.search=FALSE, as.data.frame=FALSE)
# Mollicutes use an alternative genetic code
gC2 <- getGeneticCode(id_or_name2="4", full.search=FALSE, as.data.frame=FALSE)
w <- grep("Mycoplasma|Ureaplasma", names(dna))
gC <- vector("list", length(dna))
gC[-w] <- list(gC1)
gC[w] <- list(gC2)
AA <- AlignTranslation(dna, geneticCode=gC, type="AAStringSet")
BrowseSeqs(AA)
}
