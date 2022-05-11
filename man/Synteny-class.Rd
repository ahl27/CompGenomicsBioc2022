\name{Synteny}
\alias{Synteny-class}% the class
\alias{[.Synteny}
\alias{print.Synteny}
\alias{plot.Synteny}
\alias{pairs.Synteny}
\title{
Synteny blocks and hits
}
\description{
Syntenic blocks are DNA segments composed of conserved hits occurring in the same order on two sequences.  The two sequences are typically chromosomes of different species that are hypothesized to contain homology.  Class \code{"Synteny"} provides objects and functions for storing and viewing syntenic blocks and hits that are shared between sequences.
}
\usage{
\method{pairs}{Synteny}(x,
     bounds = TRUE,
     boxBlocks = FALSE,
     labels = abbreviate(rownames(x), 9),
     gap = 0.5,
     line.main = 3,
     cex.labels = NULL,
     font.labels = 1,
     \dots)

\method{plot}{Synteny}(x,
     colorBy = 1,
     colorRamp = colorRampPalette(c("#FCF9EE", "#FFF272",
                                    "#FFAC28", "#EC5931",
                                    "#EC354D", "#0D0887")),
     barColor = "#CCCCCC",
     barSides = ifelse(nrow(x) < 100, TRUE, FALSE),
     horizontal = TRUE,
     labels = abbreviate(rownames(x), 9),
     cex.labels = NULL,
     width = 0.7,
     scaleBar = TRUE,
     \dots)

\method{print}{Synteny}(x,
      quote = FALSE,
      right = TRUE,
      \dots)
}
\arguments{
  \item{x}{
An object of class \code{Synteny}.
}
  \item{bounds}{
Logical specifying whether to plot sequence boundaries as horizontal or vertical lines.
}
  \item{boxBlocks}{
Logical indicating whether to draw a rectangle around hits belonging to the same block of synteny.
}
  \item{colorBy}{
Numeric giving the index of a reference sequence, or a character string indicating to color by ``neighbor'', ``frequency'', or ``none''.  (See details section below.)
}
  \item{colorRamp}{
A function that will return \code{n} colors when given a number \code{n}.  Examples are \code{rainbow}, \code{heat.colors}, \code{terrain.colors}, \code{cm.colors}, or (the default) \code{colorRampPalette}.
}
  \item{barColor}{
Character string giving the background color of each bar.
}
  \item{barSides}{
Logical indicating whether to draw black lines along the long-sides of each bar.
}
  \item{horizontal}{
Logical indicating whether to plot the sequences horizontally (\code{TRUE}) or vertically (\code{FALSE}).
}
  \item{labels}{
Character vector providing names corresponding to each ``identifier'' for labels on the diagonal.
}
  \item{width}{
Numeric giving the fractional width of each bar between zero and one.
}
  \item{scaleBar}{
Logical controlling whether a scale bar is drawn when \code{colorBy} is ``frequency''.  The scale bar displays the mapping between color and the level of sequence conservation.  Not applicable when \code{colorBy} is a value other than ``frequency''.
}
  \item{gap}{
Distance between subplots, in margin lines.
}
  \item{line.main}{
If \code{main} is specified, \code{line.main} provides the \code{line} argument to \code{mtext}.
}
  \item{cex.labels}{
Magnification of the labels.
}
  \item{font.labels}{
Font of labels on the diagonal.
}
  \item{quote}{
Logical indicating whether to print the output surrounded by quotes.
}
  \item{right}{
Logical specifying whether to right align strings.
}
  \item{\dots}{
Other graphical parameters for \code{pairs} or \code{plot}, including:  \code{main}, \code{cex.main}, \code{font.main}, and \code{oma}.  Other arguments for \code{print}, including \code{print.gap} and \code{max}.
}
}
\details{
Objects of class \code{Synteny} are stored as square matrices of list elements with \code{dimnames} giving the ``identifier'' of the corresponding sequences.  The synteny matrix can be separated into three parts:  along, above, and below the diagonal.  Each list element along the diagonal contains an integer vector with the width of the sequence(s) belonging to that ``identifier''.  List elements above the diagonal (column \emph{j} > row \emph{i}) each contain a \code{matrix} with ``hits'' corresponding to matches between sequences \emph{i} and \emph{j}.  List elements below the diagonal each contain a \code{matrix} with ``blocks'' of synteny between sequences \emph{j} and \emph{i}.

The \code{pairs} method creates a scatterplot matrix from a \code{Synteny} object.  Dot plots above the diagonal show hits between identifier \emph{i} and \emph{j}, where forward hits are colored in black, and hits to the reverse strand of identifier \emph{j} are colored in red.  Plots below the diagonal show blocks of synteny colored by their score, from green (highest scoring) to blue to magenta (lowest scoring).

The \code{plot} method displays a bar view of the sequences in the same order as the input object (\code{x}).  The coloring scheme of each bar is determined by the  \code{colorBy} argument, and the color palette is set by \code{colorRamp}.  When \code{colorBy} is an index, the sequences are colored according to regions of shared homology with the specified reference sequence (by default \code{1}).  If \code{colorBy} is ``neighbor'' then shared syntenic blocks are connected between neighboring sequences.  If \code{colorBy} is ``frequency'' then positions in each sequence are colored based on the degree of conservation with the other sequences.  In each case, regions that have no correspondence in the other sequence(s) are colored \code{barColor}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AlignSynteny}}, \code{\link{FindSynteny}}
}
\examples{
# a small example:
dbConn <- dbConnect(SQLite(), ":memory:")
s1 <- DNAStringSet("ACTAGACCCAGACCGATAAACGGACTGGACAAG")
s3 <- reverseComplement(s1)
s2 <- c(s1, s3)
Seqs2DB(c(c(s1, s2), s3),
        "XStringSet",
        dbConn,
        c("s1", "s2", "s2", "s3"))
syn <- FindSynteny(dbConn, minScore=1)
syn # Note:  > 100% hits because of sequence reuse across blocks
pairs(syn, boxBlocks=TRUE)
plot(syn)
dbDisconnect(dbConn)

# a larger example:
db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
synteny <- FindSynteny(db, minScore=50)
class(synteny) # 'Synteny'
synteny

# accessing parts
i <- 1
j <- 2
synteny[i, i][[1]] # width of sequences in i
synteny[j, j][[1]] # width of sequences in j
head(synteny[i, j][[1]]) # hits between i & j
synteny[j, i][[1]] # blocks between i & j

# plotting
pairs(synteny) # dot plots
pairs(synteny, boxBlocks=TRUE) # boxes around blocks

plot(synteny) # bar view colored by position in genome 1
plot(synteny, barColor="#268FD6") # emphasize missing regions
plot(synteny, "frequency") # most regions are shared by all
plot(synteny, "frequency", colorRamp=rainbow) # change the colors
plot(synteny, "neighbor") # connect neighbors
}