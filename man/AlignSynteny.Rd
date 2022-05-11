\name{AlignSynteny}
\alias{AlignSynteny}
\title{
Pairwise Aligns Syntenic Blocks
}
\description{
Performs pairwise alignment of all blocks of synteny between sets of sequences.
}
\usage{
AlignSynteny(synteny,
             dbFile,
             tblName = "Seqs",
             identifier = "",
             processors = 1,
             verbose = TRUE,
             ...)
}
\arguments{
  \item{synteny}{
An object of class ``Synteny''.
}
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table where the sequences are located that were used to create the object \code{synteny}.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier, or an integer sequence corresponding to indices of \code{rownames(synteny)}.  If "" (the default), then all identifiers are selected from \code{synteny}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{AlignProfiles}}, including \code{perfectMatch}, \code{misMatch}, \code{gapPower}, \code{terminalGap}, \code{restrict}, \code{normPower}, and \code{substitutionMatrix}.
}
}
\details{
\code{AlignSynteny} will extract all sequence regions belonging to syntenic blocks in \code{synteny}, and perform pairwise alignment with \code{AlignProfiles}.  Hits are used to anchor the alignment such that only the regions between anchors are aligned.
}
\value{
A list with elements for each pair of \code{identifier}s in \code{synteny}.  Each list element contains a \code{DNAStringSetList} one pairwise alignment per syntenic block.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{FindSynteny}}, \code{\link{Synteny-class}}
}
\examples{
db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
synteny <- FindSynteny(db, minScore=50)
DNA <- AlignSynteny(synteny, db)
names(DNA)
DNA[[1]] # the first set of pairwise alignments
DNA[[1]][[1]] # the first block of synteny between H9N2 & H5N1
unlist(DNA[[2]]) # a DNAStringSet of synteny between H9N2 & H2N2
}