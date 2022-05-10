\name{Seqs2DB}
\alias{Seqs2DB}
\title{
Add Sequences from Text File to Database
}
\description{
Adds sequences to a database.
}
\usage{
Seqs2DB(seqs,
        type,
        dbFile,
        identifier,
        tblName = "Seqs",
        chunkSize = 1e7,
        replaceTbl = FALSE,
        fields = c(accession = "ACCESSION", rank = "ORGANISM"),
        processors = 1,
        verbose = TRUE,
        \dots)
}
\arguments{
  \item{seqs}{
A connection object or a character string specifying the file path to the file containing the sequences, an \code{XStringSet} object if \code{type} is \code{XStringSet}, or a \code{QualityScaledXStringSet} object if \code{type} is \code{QualityScaledXStringSet}. Files compressed with \code{gzip}, \code{bzip2}, \code{xz}, or \code{lzma} compression are automatically detected and decompressed during import.  Full URL paths (e.g., "http://" or "ftp://") to uncompressed text files or \code{gzip} compressed text files can also be used.
}
  \item{type}{
The type of the sequences (\code{seqs}) being imported.  This should be (an unambiguous abbreviation of) one of \code{"FASTA"}, \code{"FASTQ"}, \code{"GenBank"}, \code{"XStringSet"}, or \code{"QualityScaledXStringSet"}.
}
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.  If the \code{dbFile} does not exist then a new database is created at this location.
}
  \item{identifier}{
Character string specifying the \code{"id"} to give the imported sequences in the database.
}
  \item{tblName}{
Character string specifying the table in which to add the sequences.
}
  \item{chunkSize}{
Number of characters to read at a time.
}
  \item{replaceTbl}{
Logical indicating whether to overwrite the entire table in the database.  If \code{FALSE} (the default) then the sequences are appended to any already existing in the \code{tblName}.  If \code{TRUE} the entire table is dropped, removing any existing sequences before adding any new sequences.
}
  \item{fields}{
Named character vector providing the fields to import from a \code{"GenBank"} formatted file as text columns in the database (not applicable for other \code{"type"}s).  The default is to import the \code{"ACCESSION"} field as a column named \code{"accession"} and the \code{"ORGANISM"} field as a column named \code{"rank"}.  Other uppercase fields, such as \code{"LOCUS"} or \code{"VERSION"}, can be specified in similar manner.  Note that the \code{"DEFINITION"} field is automatically imported as a column named \code{"description"} in the database.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display each query as it is sent to the database.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{Codec}} for compressing sequence information.
}
}
\details{
Sequences are imported into the database in chunks of lines specified by \code{chunkSize}.  The sequences can then be identified by searching the database for the \code{identifier} provided.  Sequences are added to the database verbatim, so that no sequence information is lost when the sequences are exported from the database.  The sequence (record) names are recorded into a column named ``description'' in the database.
}
\value{
The total number of sequences in the database table is returned after import.
}
\references{
ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\section{Warning}{
If \code{replaceTbl} is \code{TRUE} then any sequences already in the table are overwritten, which is equivalent to dropping the entire table.
}
\seealso{
\code{\link{BrowseDB}}, \code{\link{SearchDB}}, \code{\link{DB2Seqs}}
}
\examples{
gen <- system.file("extdata", "Bacteria_175seqs.gen", package="DECIPHER")
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(gen, "GenBank", dbConn, "Bacteria")
BrowseDB(dbConn)
dna <- SearchDB(dbConn, nameBy="description")
dbDisconnect(dbConn)
}