\name{BrowseSeqs}
\alias{BrowseSeqs}
\title{
View Sequences in a Web Browser
}
\description{
Opens an html file in a web browser to show the sequences in an \code{XStringSet}.
}
\usage{
BrowseSeqs(myXStringSet,
           htmlFile = paste(tempdir(), "/myXStringSet.html", sep = ""),
           openURL = interactive(),
           colorPatterns = TRUE,
           highlight = NA,
           patterns = c("-", alphabet(myXStringSet, baseOnly=TRUE)),
           colors = substring(rainbow(length(patterns),
                              v=0.8, start=0.9, end=0.7), 1, 7),
           colWidth = Inf,
           \dots)
}
\arguments{
  \item{myXStringSet}{
A \code{XStringSet} object of sequences.
}
  \item{htmlFile}{
Character string giving the location where the html file should be written.
}
  \item{openURL}{
Logical indicating whether the \code{htmlFile} should be opened in a web browser.
}
  \item{colorPatterns}{
Logical specifying whether to color matched \code{patterns}, or an integer vector providing pairs of start and stop boundaries for coloring.
}
  \item{highlight}{
Numeric specifying which sequence in the set to use for comparison or \code{NA} to color all sequences (default).  If \code{highlight} is \code{0} then positions differing from the consensus sequence are highlighted.
}
  \item{patterns}{
Either an \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object, a character vector containing regular expressions, a list of numeric matrices, or \code{NULL}.  (See details section below.)
}
  \item{colors}{
Character vector providing the color for each of the matched \code{patterns}.  Typically a character vector with elements of 7 characters:  ``#'' followed by the red, blue, green values in hexadecimal (after rescaling to 0 ... 255).  Ignored when \code{patterns} is a list of matrices.
}
  \item{colWidth}{
Integer giving the maximum number of nucleotides wide the display can be before starting a new page.  Must be a multiple of \code{20} (e.g., \code{100}), or \code{Inf} (the default) to display all the sequences in one set of rows.
}
  \item{\dots}{
Additional arguments to adjust the appearance of the consensus sequence at the base of the display.  Passed directly to \code{ConsensusSequence} for an \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet}, or to \code{consensusString} for a \code{BStringSet}.
}
}
\details{
\code{BrowseSeqs} converts an \code{XStringSet} into html format for viewing in a web browser.  The sequences are colored in accordance with the \code{patterns} that are provided, or left uncolored if \code{colorPatterns} is \code{FALSE} or \code{patterns} is \code{NULL}.  Character or \code{XStringSet} patterns are matched as regular expressions and colored according to \code{colors}.  If \code{patterns} is a list of matrices, then it must contain one element per sequence.  Each matrix is interpreted as providing the fraction red, blue, and green for each letter in the sequence.  Thus, \code{colors} is ignored when \code{patterns} is a list.  (See examples section below.)

Patterns are not matched across column breaks, so multi-character \code{patterns} should be carefully considered when \code{colWidth} is less than the maximum sequence length.  Patterns are matched sequentially in the order provided, so it is feasible to use nested \code{patterns} such as \code{c("ACCTG", "CC")}.  In this case the ``CC'' could be colored differently inside the previously colored ``ACCTG''.  Note that \code{patterns} overlapping the boundaries of a previously matched pattern will not be matched.  For example, ``ACCTG'' would not be matched if \code{patterns=c("CC", "ACCTG")}.

Some web browsers cannot quickly display a large amount colored text, so it is recommended to use \code{colorPatterns = FALSE} or to \code{highlight} a sequence when viewing a large \code{XStringSet}.  Highlighting will only show all of the characters in the highlighted sequence, and convert all matching positions in the other sequences into dots without \code{color}.  Also, note that some web browsers display small shifts between fixed-width characters that may become noticeable as color offsets between the ends of long sequences.
}
\value{
Creates an html file containing sequence data and (if \code{openURL} is \code{TRUE}) opens it in a web browser for viewing.  The layout has the sequence name on the left, position legend on the top, cumulative number of nucleotides on the right, and consensus sequence on the bottom.

Returns \code{htmlFile} if the html file was written successfully.
}
\references{
ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
Kunzmann P., et al. (2020) "Substitution matrix based color schemes for sequence alignment visualization". BMC Bioinformatics, \bold{21(1):209}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\note{
Some web browsers do not display colored characters with equal widths.  If positions do not align across sequences then try opening the \code{htmlFile} with a different web browser.
}
\seealso{
\code{\link{BrowseDB}}, \code{\link{ConsensusSequence}}
}
\examples{
# load the example DNA sequences
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
dna <- SearchDB(db) # non-coding ribosomal RNA gene sequences

# example of using the defaults with DNA sequences
BrowseSeqs(dna) # view the XStringSet

# color only "ACTG" and "CSC" patterns (where S is C or G)
BrowseSeqs(dna, patterns=DNAStringSet(c("ACTG", "CSC")))

# highlight (i.e., only fully-color) the first sequence
BrowseSeqs(dna, highlight=1) # other sequences are dots where matching

# highlight the consensus sequence at the bottom
BrowseSeqs(dna, highlight=0) # other sequences are dots where matching

# split the wide view into multiple vertical pages (for printing)
BrowseSeqs(dna, colWidth=100, highlight=1)

# specify an alternative color scheme for -, A, C, G, T
BrowseSeqs(dna, colors=c("#1E90FF", "#32CD32", "#9400D3", "black", "#EE3300"))

# only color the positions within certain positional ranges (100-200 & 250-500)
BrowseSeqs(dna, colorPatterns=c(100, 200, 250, 500))

# example of calling attention to letters by coloring gaps black
BrowseSeqs(dna, patterns="-", colors="black")

# color according to base-pairing by supplying the fraction RGB for every position
dbn <- PredictDBN(dna, type="structures") # calculate the secondary structures
# dbn now contains the scores for whether a base is paired (left/right) or unpaired
dbn[[1]][, 1] # the scores for the first position in the first sequence
dbn[[2]][, 10] # the scores for the tenth position in the second sequence
# these positional scores can be used as shades of red, green, and blue:
BrowseSeqs(dna, patterns=dbn) # red = unpaired, green = left-pairing, blue = right
# positions in black are not part of the consensus secondary structure

# color all restriction sites
data(RESTRICTION_ENZYMES) # load dataset containing restriction enzyme sequences
sites <- RESTRICTION_ENZYMES
sites <- gsub("[^A-Z]", "", sites) # remove non-letters
sites <- DNAStringSet(sites) # convert the character vector to a DNAStringSet
rc_sites <- reverseComplement(DNAStringSet(sites))
w <- which(sites != rc_sites) # find non-palindromic restriction sites
sites <- c(sites, rc_sites[w]) # append their reverse complement
sites <- sites[order(nchar(sites))] # match shorter sites first
BrowseSeqs(dna, patterns=sites)

# color bases by quality score
fastq <- system.file("extdata", "s_1_sequence.txt", package="Biostrings")
reads <- readQualityScaledDNAStringSet(fastq, quality.scoring="solexa")
colors <- colorRampPalette(c("red", "yellow", "green"))(42)
colors <- col2rgb(colors)/255
quals <- as(quality(reads), "IntegerList")
quals <- lapply(quals, function(x) colors[, x])
BrowseSeqs(DNAStringSet(reads), patterns=quals) # green = high quality, red = low quality

# load the example protein coding sequences
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

# example of using the defaults with amino acid sequences
aa <- unique(translate(dna)) # the unique amino acid sequences
BrowseSeqs(aa)

# example of highlighting the consensus amino acid sequence
AA <- AlignSeqs(aa)
BrowseSeqs(AA, highlight=0)

# example of highlighting positions that differ from the majority consensus
BrowseSeqs(AA, highlight=0, threshold=0.5)

# specify an alternative color scheme for amino acids (from Kunzmann et al.)
colors <- c(`-`="#000000", `A`="#BDB1E8", `R`="#EFA2C5", `N`="#F6602F",
    `D`="#FD5559", `C`="#12C7FE", `Q`="#DDACB4", `E`="#FEA097", `G`="#F46802",
    `H`="#FCA708", `I`="#369BD9", `L`="#2E95EC", `K`="#CF7690", `M`="#4B8EFE",
    `F`="#76997D", `P`="#FD2AE3", `S`="#A08A9A", `T`="#9A84D5", `W`="#74C80D",
    `Y`="#9BB896", `V`="#89B9F9")
BrowseSeqs(AA, colors=colors, patterns=names(colors))

# example of coloring in a reduced amino acid alphabet
alpha <- AA_REDUCED[[15]]
alpha # clustering of amino acids based on similarity
BrowseSeqs(AA, patterns=c("-", paste("[", alpha, "]", sep="")))

# color amino acids according to their predicted secondary structure
hec <- PredictHEC(AA, type="probabilities") # calculate the secondary structures
# hec now contains the probability that a base is in an alpha-helix or beta-sheet
hec[[3]][, 18] # the 18th position in sequence 3 is likely part of a beta-sheet (E)
# the positional probabilities can be used as shades of red, green, and blue:
BrowseSeqs(AA, patterns=hec) # red = alpha-helix, green = beta-sheet, blue = coil

# color codons according to their corresponding amino acid
DNA <- AlignTranslation(dna) # align the translation then reverse translate
colors <- rainbow(21, v=0.8, start=0.9, end=0.7) # codon colors
m <- match(GENETIC_CODE, unique(GENETIC_CODE)) # corresponding amino acid
codonBounds <- matrix(c(seq(1, width(DNA)[1], 3), # start of codons
	seq(3, width(DNA)[1], 3)), # end of codons
	nrow=2,
	byrow=TRUE)
BrowseSeqs(DNA,
	colorPatterns=codonBounds,
	patterns=c("---", names(GENETIC_CODE)), # codons to color
	colors=c("black", substring(colors[m], 1, 7)))
}