\name{GetProtWebData.ProtWeb}
\alias{GetProtWebData.ProtWeb}
\alias{GetProtWebData}
\title{
Extract information from a ProtWeb object
}
\description{
ProtWeb objects are outputted from \code{\link{predict.ProtWeaver}}. 

This function extracts the underlying data from the object.
}
\usage{
\method{GetProtWebData}{ProtWeb}(x, AsDf=FALSE, ...)
}
\arguments{
  \item{x}{
    A ProtWeb object
  }
  \item{AsDf}{
    Should data be printed as a pairwise entry? If \code{TRUE}, returns a matrix 
    with three columns, where the first two columns define the pair of genes/proteins 
    and the third column defines the prediction/score. If \code{FALSE}, returns an 
    adjacency matrix encoding the same information (but may be sparse depending on
    how many predictions were made in the original \code{predict} call.)
  }
  \item{...}{
    Additional parameters for consistency with generic.
  }
}

\details{
\code{\link{predict.ProtWeaver}} returns a \code{ProtWeb} object, which bundles some methods
to make formatting and printing of results slightly nicer. This method extracts data form the \code{ProtWeb} object.

If \code{AsDf=TRUE}, the return data is a \code{3xN} \code{\link[base]{data.frame}}, with columns 3 showing
the prediction for the pair of genes/proteins specified in columns 1 and 2. This format 
is when predictions are made on small number of pairs (meaning the resulting adjacency 
matrix is sparse).

If \code{AsDf=FALSE}, the return data is a \code{NxN} adjacency matrix, with entry
\code{i,j} containing the prediction for genes \code{i} and \code{j}. 
}
\value{
Either a \code{\link[base]{data.frame}} or a matrix.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{predict.ProtWeaver}}
}
\examples{
##############
## Prediction with built-in model and data 
###############

exData <- get(data("ExampleStreptomycesData"))

# Subset isn't necessary but is faster for a working example
pw <- ProtWeaver(exData$Genes[1:10])

protweb <- predict(pw, method='Jaccard')

# print out results as an adjacency matrix
GetProtWebData(protweb)

# print out results as a pairwise data.frame
GetProtWebData(protweb, AsDf=TRUE)
}