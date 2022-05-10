\name{FindGenes}
\alias{FindGenes}
\title{
Find Genes in a Genome
}
\description{
Predicts the start and stop positions of protein coding genes in a genome.
}
\usage{
FindGenes(myDNAStringSet,
          geneticCode = getGeneticCode("11"),
          minGeneLength = 60,
          includeGenes = NULL,
          allowEdges = TRUE,
          allScores = FALSE,
          showPlot = FALSE,
          processors = 1,
          verbose = TRUE)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} object of unaligned sequences representing a genome.
}
  \item{geneticCode}{
A named character vector defining the translation from codons to amino acids.  Optionally, an \code{"alt_init_codons"} attribute can be used to specify alternative initiation codons.  By default, the bacterial and archael genetic code is used, which has seven possible initiation codons: ATG, GTG, TTG, CTG, ATA, ATT, and ATC.
}
  \item{minGeneLength}{
Integer specifying the minimum length of genes to find in the genome.
}
  \item{includeGenes}{
A \code{Genes} object to include as potential genes or \code{NULL} (the default) to predict all genes \code{de novo}.
}
  \item{allowEdges}{
Logical determining whether to allow genes that run off the edge of the sequences.  If \code{TRUE} (the default), genes can be identified with implied starts or ends outside the boundaries of \code{myDNAStringSet}, although the boundary will be set to the last possible codon position.  This is useful when genome sequences are circular or incomplete.
}
  \item{allScores}{
Logical indicating whether to return information about all possible open reading frames or only the predicted genes (the default).
}
  \item{showPlot}{
Logical determining whether a plot is displayed with the distribution of gene lengths and scores.  (See details section below.)
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to print information about the predictions on each iteration.  (See details section below.)
}
}
\details{
Protein coding genes are identified by learning their characteristic signature directly from the genome, i.e., \emph{ab initio} prediction.  Gene signatures are derived from the content of the open reading frame and surrounding signals that indicate the presence of a gene.  Genes are assumed to not contain introns or frame shifts, making the function best suited for prokaryotic genomes.

If \code{showPlot} is \code{TRUE} then a plot is displayed with four panels.  The upper left panel shows the fitted distribution of background open reading frame lengths.  The upper right panel shows this distribution on top of the fitted distribution of predicted gene lengths.  The lower left panel shows the fitted distribution of scores for the intergenic spacing between genes on the same and opposite genome strands.  The bottom right panel shows the total score of open reading frames and predicted genes by length.

If \code{verbose} is \code{TRUE}, information is shown about the predictions at each iteration of gene finding.  The mean score difference between genes and non-genes is updated at each iteration, unless it is negative, in which case the score is dropped and a \code{"-"} is displayed.  The columns denote the number of iterations (\code{"Iter"}), number of codon scoring models (\code{"Models"}), start codon scores (\code{"Start"}), upstream k-mer motif scores (\code{"Motif"}), mRNA folding scores (\code{"Fold"}), initial codon bias scores (\code{"Init"}), upstream nucleotide bias scores (\code{"UpsNt"}), termination codon bias scores (\code{"Term"}), ribosome binding site scores (\code{"RBS"}), codon autocorrelation scores (\code{"Auto"}), stop codon scores (\code{"Stop"}), and number of predicted genes (\code{"Genes"}).
}
\value{
An object of class \code{Genes}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{ExtractGenes}}, \code{\link{FindNonCoding}}, \code{\link{Genes-class}}, \code{\link{WriteGenes}}
}
\examples{
# import a test genome
fas <- system.file("extdata",
	"Chlamydia_trachomatis_NC_000117.fas.gz",
	package="DECIPHER")
genome <- readDNAStringSet(fas)

z <- FindGenes(genome)
z
genes <- ExtractGenes(z, genome)
genes
proteins <- ExtractGenes(z, genome, type="AAStringSet")
proteins
}
