\name{LearnNonCoding}
\alias{LearnNonCoding}
\title{
Learn a Non-Coding RNA Model
}
\description{
Learns a compact representation of patterns representing a set of non-coding RNAs belonging to the same family.
}
\usage{
LearnNonCoding(myXStringSet,
               threshold = 0.3,
               weight = NA,
               maxLoopLength = 500,
               maxPatterns = 20,
               scoreDependence = FALSE,
               structure = NULL,
               processors = 1)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} object of aligned sequence representatives belonging to the same non-coding RNA family.
}
  \item{threshold}{
Numeric specifying the minimum relative frequency of patterns to consider during learning.
}
  \item{weight}{
Either a numeric vector of weights for each sequence, a single number implying equal weights, or \code{NA} (the default) to automatically calculate sequence weights based on \code{myXStringSet}.
}
  \item{maxLoopLength}{
Numeric giving the maximum length of conserved hairpin loops to consider.
}
  \item{maxPatterns}{
A numeric vector of length two specifying the maximum number of motifs and hairpins, respectively, or a single numeric giving the maximum for each.
}
  \item{scoreDependence}{
Logical determining whether to record a log-odds score for dependencies between patterns.  The default (\code{FALSE}) is recommended for most non-coding RNA families.
}
  \item{structure}{
Either a character string providing the consensus secondary structure in dot bracket notation, a matrix of paired positions in the first two columns, or \code{NULL} (the default) to predict the consensus secondary structure with \code{PredictDBN}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
}
\details{
Non-coding RNAs belonging to the same family typically have conserved sequence motifs, secondary structure elements, and k-mer frequencies that can be used to identify members of the family.  \code{LearnNonCoding} identifies these conserved patterns and determines which are best for identifying the non-coding RNA relative to a random sequence background.  Sequence motifs and hairpins are defined relative to their distance from the start or end of the non-coding RNA, allowing the precise and rapid identification of the boundaries of any matches to the non-coding RNA in a genome.
}
\value{
An object of class \code{NonCoding}.
}
\references{
Wright, E. S. (2021). FindNonCoding: rapid and simple detection of non-coding RNAs in genomes. Bioinformatics. https://doi.org/10.1093/bioinformatics/btab708
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{FindNonCoding}}, \code{\link{NonCoding-class}}
}
\examples{
# import a family of non-coding RNAs
fas_path <- system.file("extdata",
	"IhtA.fas",
	package="DECIPHER")
rna <- readRNAStringSet(fas_path)
rna

# align the sequences
RNA <- AlignSeqs(rna)
RNA

y <- LearnNonCoding(RNA)
y
y[["motifs"]]
y[["hairpins"]]
head(y[["kmers"]])
}
