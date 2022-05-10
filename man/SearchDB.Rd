\name{SearchDB}
\alias{SearchDB}
\title{
Obtain Specific Sequences from a Database
}
\description{
Returns the set of sequences meeting the search criteria.
}
\usage{
SearchDB(dbFile,
         tblName = "Seqs",
         identifier = "",
         type = "XStringSet",
         limit = -1,
         replaceChar = NA,
         nameBy = "row_names",
         orderBy = "row_names",
         countOnly = FALSE,
         removeGaps = "none",
         quality = "Phred",
         clause = "",
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
Optional character string used to narrow the search results to those matching a specific identifier.  If "" (the default) then all identifiers are selected.
}
  \item{type}{
The type of \code{XStringSet} (sequences) to return.  This should be (an unambiguous abbreviation of) one of \code{"XStringSet"}, \code{"DNAStringSet"}, \code{"RNAStringSet"}, \code{"AAStringSet"}, \code{"BStringSet"}, \code{"QualityScaledXStringSet"}, \code{"QualityScaledDNAStringSet"}, \code{"QualityScaledRNAStringSet"}, \code{"QualityScaledAAStringSet"}, or \code{"QualityScaledBStringSet"}.  If \code{type} is \code{"XStringSet"} or \code{"QualityScaledXStringSet"} then an attempt is made to guess the type of sequences based on their composition.
}
  \item{limit}{
Number of results to display.  The default (\code{-1}) does not limit the number of results.
}
  \item{replaceChar}{
Optional character used to replace any characters of the sequence that are not present in the \code{XStringSet}'s alphabet.  Not applicable if \code{type=="BStringSet"}.  The default (\code{NA}) results in an error if an incompatible character exist.  (See details section below.)
}
  \item{nameBy}{
Character string giving the column name for naming the \code{XStringSet}.
}
  \item{orderBy}{
Character string giving the column name for sorting the results.  Defaults to the order of entries in the database.  Optionally can be followed by \code{" ASC"} or \code{" DESC"} to specify ascending (the default) or descending order.
}
  \item{countOnly}{
Logical specifying whether to return only the number of sequences.
}
  \item{removeGaps}{
Determines how gaps ("-" or "." characters) are removed in the sequences.  This should be (an unambiguous abbreviation of) one of \code{"none"}, \code{"all"} or \code{"common"}.
}
  \item{clause}{
An optional character string to append to the query as part of a ``where clause''.
}
  \item{quality}{
The type of quality object to return if \code{type} is a \code{QualityScaledXStringSet}.  This should be (an unambiguous abbreviation of) one of \code{"Phred"}, \code{"Solexa"}, or \code{"Illumina"}.  Note that recent versions of Illumina software provide \code{"Phred"} formatted quality scores.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display queries as they are sent to the database.
}
}
\details{
If \code{type} is \code{"DNAStringSet"} then all U's are converted to T's before creating the \code{DNAStringSet}, and vise-versa if \code{type} is \code{"RNAStringSet"}.  All remaining characters not in the \code{XStringSet}'s alphabet are converted to \code{replaceChar} or removed if \code{replaceChar} is \code{""}.  Note that if \code{replaceChar} is \code{NA} (the default), it will result in an error when an unexpected character is found.  Quality information is interpreted as \code{PredQuality} scores.
}
\value{
An \code{XStringSet} or \code{QualityScaledXStringSet} with the sequences that meet the specified criteria.  The \code{names} of the object correspond to the value in the \code{nameBy} column of the database.
}
\references{
ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{Seqs2DB}}, \code{\link{DB2Seqs}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
# get all sequences in the default table:
dna <- SearchDB(db)
# select a random sequence:
dna <- SearchDB(db, orderBy="random()", limit=1)
# remove gaps from "Mycobacterium" sequences:
dna <- SearchDB(db, identifier="Mycobacterium", removeGaps="all")
# provide a more complex query:
dna <- SearchDB(db, nameBy="description", orderBy="bases", removeGaps="common",
                clause="nonbases is 0")
}