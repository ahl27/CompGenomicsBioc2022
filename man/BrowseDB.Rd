\name{BrowseDB}
\alias{BrowseDB}
\title{
View a Database Table in a Web Browser
}
\description{
Opens an html file in a web browser to show the contents of a table in a database.
}
\usage{
BrowseDB(dbFile,
         htmlFile = paste(tempdir(), "/db.html", sep = ""),
         openURL = interactive(),
         tblName = "Seqs",
         identifier = "",
         limit = -1,
         orderBy = "row_names",
         maxChars = 50,
         clause="")
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{htmlFile}{
Character string giving the location where the html file should be written.
}
  \item{openURL}{
Logical indicating whether the \code{htmlFile} should be opened in a web browser.
}
  \item{tblName}{
Character string specifying the table to view.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.
}
  \item{limit}{
Number of results to display.  The default (\code{-1}) does not limit the number of results.
}
  \item{orderBy}{
Character string giving the column name for sorting the results.  Defaults to the order of entries in the database.  Optionally can be followed by \code{" ASC"} or \code{" DESC"} to specify ascending (the default) or descending order.
}
  \item{maxChars}{
Maximum number of characters to display in each column.
}
  \item{clause}{
An optional character string to append to the query as part of a ``where clause''.
}
}
\value{
Creates an html table containing all the fields of the database table and (if \code{openURL} is \code{TRUE}) opens it in the web browser for viewing.

Returns \code{htmlFile} if the html file was written successfully.
}
\references{
ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\note{
If viewing a table containing sequences, the sequences are purposefully not shown in the output.
}
\seealso{
\code{\link{BrowseSeqs}}
}
\examples{
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
BrowseDB(db)
}