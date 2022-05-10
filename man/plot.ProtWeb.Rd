\name{plot.ProtWeb}
\alias{plot.ProtWeb}
\title{
Plot predictions in a ProtWeb object
}
\description{
ProtWeb objects are outputted from \code{\link{predict.ProtWeaver}}. 

This function plots the predictions in the object using a force-directed embedding
of connections in the adjacency matrix.

\emph{This function is still a work in progress.}
}
\usage{
\method{plot}{ProtWeb}(x, NumSims=10, 
      Gravity=0.05, Coulomb=0.1, Connection=5, 
      MoveRate=0.25, Cutoff=0.2, ColorPalette=topo.colors, 
      Verbose=TRUE, ...)
}
\arguments{
  \item{x}{
    A ProtWeb object. See \code{\link{predict.ProtWeaver}}
  }
  \item{NumSims}{
    Number of iterations to run the model for.
  }
  \item{Gravity}{
    Strength of Gravity force. See 'Details'.
  }
  \item{Coulomb}{
    Strength of Coulomb force. See 'Details'.
  }
  \item{Connection}{
    Strength of Connective force. See 'Details'.
  }
  \item{MoveRate}{
    Controls how far each point moves in each iteration.
  }
  \item{Cutoff}{
    Cutoff value; if \code{abs(val) < Cutoff}, that Connection is shrunk to zero.
  }
  \item{ColorPalette}{
    Color palette for graphing. Valid inputs are any palette available in
    \code{palette.pals()}. See \code{\link[grDevices]{palette}} for more info.
  }
  \item{Verbose}{
    Logical indicating whether to print progress bars and messages. Defaults to \code{TRUE}.
  }
  \item{...}{
    Additional parameters for consistency with generic.
  }
}

\details{
  This function plots the \code{ProtWeb} object using a force-directed embedding.
  This embedding has three force components:
  \itemize{
  \item Gravity Force: Attractive force pulling nodes towards \code{(0,0)}
  \item Coulomb Force: Repulsive force pushing close nodes away from each other
  \item Connective Force: Tries to push node connections to equal corresponding values in the adjacency matrix
  }
  
  The parameters in the function are sufficient to get an embedding, though users
  are welcome to try to tune them for a better visualization. This function is meant
  to aid with visualization of the adjacency matrix, not for concrete analyses of clusters.
  
  The function included in this release is early stage.
  Next release cycle will update this function with an updated version of this 
  algorithm to improve plotting, visualization, and runtime. 
}
\value{
No return value; creates a plot in the graphics window.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\seealso{
\code{\link{predict.ProtWeaver}}

\code{\link{GetProtWebData.ProtWeb}}
}
\examples{
exData <- get(data("ExampleStreptomycesData"))
pw <- ProtWeaver(exData$Genes)

# Subset isn't necessary but is faster for a working example
# Same w/ method='Jaccard'
protweb <- predict(pw, 'Jaccard', subset=1:50)

plot(protweb)
}