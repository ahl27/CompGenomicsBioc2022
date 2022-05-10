\name{predict.ProtWeaver}
\alias{predict.ProtWeaver}
\title{
Make predictions with ProtWeaver objects
}
\description{
This S3 method predicts a functional association network from a ProtWeaver object.
This returns an object of type \code{ProtWeb}, which is essentially an adjacency
matrix with some extra S3 methods to make printing cleaner.
}
\usage{
\method{predict}{ProtWeaver}(object, Method='Ensemble',
         Subset=NULL, NumCores=1,
         MySpeciesTree=NULL, PretrainedModel=NULL,
         RawZScores=FALSE, NoPrediction=FALSE,
         ReturnRawData=FALSE, Verbose=TRUE, ...)
}
\arguments{
  \item{object}{
    A ProtWeaver object
  }
  \item{Method}{
    Method to use for prediction. See 'Details'.
  }
  \item{Subset}{
    Subset of data to predict on. This can either be a vector or a \code{2xN} matrix.
    
    
    If a vector, prediction proceeds for all possible pairs of elements specified in the vector 
    (either by name, for character vector, or by index, for numeric vector). For example, 
    \code{subset=1:3} will predict for pairs \code{(1,2), (1,3), (2,3)}. 
    
    If a matrix, subset is interpreted as a matrix of pairs, where each row of the matrix 
    specifies a pair to evaluate. These can also be specifed by name (character) or by index 
    (numeric). 
    
    \code{subset=cbind(c(1,1,2), c(2,3,3))} produces equivalent functionality to 
    \code{subset=1:3}.
  }
  \item{NumCores}{
    Number of cores to use for methods that support multithreaded execution. 
    Currently only supported for methods \code{'ProfDCA'} and \code{'Ensemble'}. Setting
    to a negative value will use one less than the value of \code{detectCores()}, 
    or one core if the number of available cores cannot be determined. See Note
    for more information. This parameter has no effect on Windows due to reliance on
    forking via \code{mclapply}.
  }
  \item{MySpeciesTree}{
    Phylogenetic tree of all genomes in the dataset. Required for \code{Method='Behdenna'},
    and can improve predictions for other methods. \code{'Behdenna'} requires a 
    rooted, bifurcating tree (other values of \code{Method} can handle arbitrary trees).
  }
  \item{PretrainedModel}{
    A pretrained model for use with ensemble predictions. If unspecified when 
    \code{Method='Ensemble'}, the program will use built-in models 
    (see \code{\link{BuiltInEnsembles}}). See the examples for how to train an
    ensemble method to pass to \code{PretrainedModel}. 
    
    Has no effect if \code{Method != 'Ensemble'}.
  }
  \item{RawZScores}{
    For methods that return z-scores, should raw scores be returned? If FALSE, instead
    returns normalized absolute value of predictions. These tend to be better predictions.
    
    Currently, only \code{Method='Behdenna'} uses this parameter.
  }
  \item{NoPrediction}{
    For \code{Method='Ensemble'}, should data be returned prior to making predictions?
    
    If \code{TRUE}, this will instead return a dataframe (\link[base]{data.frame})
    with predictions from each algorithm for each pair. This dataframe is typically
    used to train an ensemble model. 
    
    If \code{FALSE}, ProtWeaver will return predictions for each pair (using user model
    if provided or a built-in otherwise).
  }
  \item{ReturnRawData}{
      Internal parameter used for ensemble predictions. If TRUE, returns predictions
      without formatting them into a ProtWeb object. Users should specify 
      \code{NoPrediction=TRUE} rather than use this parameter (see Details).
  }
  \item{Verbose}{
    Logical indicating whether to print progress bars and messages. Defaults to \code{TRUE}.
  }
  \item{...}{
    Additional parameters for other predictors and consistency with generic.
  }
}

\details{
\code{predict.ProtWeaver} wraps several methods to create an easy interface for multiple prediction types. The following values of \code{Method} are currently supported: 

\itemize{
  \item \code{'Jaccard'}: Jaccard distance of PA profiles
  \item \code{'Hamming'}: Hamming distance of PA profiles
  \item \code{'MutualInformation'}: MI of PA profiles
  \item \code{'ProfDCA'}: Direct Coupling Analysis of PA profiles
  \item \code{'Behdenna'}: Analysis of Gain/Loss events following Behdenna et al. (2016) 
  \item \code{'Coloc'}: Co-localization analysis
  \item \code{'MirrorTree'}: MirrorTree
  \item \code{'ContextTree'}: ContextTree
}

(PA = Presence/Absence)

This returns a \code{ProtWeb} object, an S3 class that makes formatting and printing of 
results slightly nicer. Data can be extracted form the \code{ProtWeb} object with:

\code{GetProtWebData(ProtWebObject, AsDf=c(T,F))}

Different methods require different types of input. The constructor 
\code{\link{ProtWeaver}} will notify the user which methods are
runnable with the given data. Note that method \code{Behdenna} requires a species 
tree, which must be bifurcating. Method \code{Ensemble} automatically selects the 
methods that can be run with the given input data. 

See \code{\link{ProtWeaver}} for more information on input data types.
}
\value{
Returns a ProtWeb object. See \code{\link[=GetProtWebData.ProtWeb]{GetProtWebData}} for more info.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}

\references{
Behdenna, A., et al., \emph{Testing for Independence between Evolutionary Processes.} Systematic Biology, 2016. \bold{65}(5): p. 812-823.

Franceschini, A., et al., \emph{SVD-phy: improved prediction of protein functional associations through singular value decomposition of phylogenetic profiles.} Bioinformatics, 2016. \bold{32}(7): p. 1085-1087.

Fukunaga, T. and W. Iwasaki, \emph{Inverse Potts model improves accuracy of phylogenetic profiling.} Bioinformatics, 2022.

Lokhov, A.Y., et al., \emph{Optimal structure and parameter learning of Ising models.} Science advances, 2018. \bold{4}(3): p. e1700791.

Pazos, F. and A. Valencia, \emph{Similarity of phylogenetic trees as indicator of protein–protein interaction.} Protein Engineering, Design and Selection, 2001. \bold{14}(9): p. 609-614.

Pazos, F., et al., \emph{Assessing protein co-evolution in the context of the tree of life assists in the prediction of the interactome.} J Mol Biol, 2005. \bold{352}(4): p. 1002-15.

Sadreyev, I.R., et al., \emph{PhyloGene server for identification and visualization of co-evolving proteins using normalized phylogenetic profiles.} Nucleic Acids Research, 2015. \bold{43}(W1): p. W154-W159.

Sato, T., et al., \emph{The inference of protein-protein interactions by co-evolutionary analysis is improved by excluding the information about the phylogenetic relationships.} Bioinformatics, 2005. \bold{21}(17): p. 3482-9.

Sato, T., et al., \emph{Partial correlation coefficient between distance matrices as a new indicator of protein-protein interactions.} Bioinformatics, 2006. \bold{22}(20): p. 2488-92.

Weigt, M., et al., \emph{Identification of direct residue contacts in protein-protein interaction by message passing.} Proceedings of the National Academy of Sciences, 2009. \bold{106}(1): p. 67-72.

}

\note{
  Note that the pairwise associations are stored in a matrix, meaning that if the
  ProtWeaver object contains 100 entries, the output ProtWeb object contains a
  \code{100x100} matrix. Users should be advised that predicting too many pairs
  can lead to vector memory exhaustion errors. On my machine, ProtWeaver supports 
  predictions on up to around \eqn{53,000^2} pairs (meaning all possible pairs between
  53,000 proteins/genes). Next release will add more memory efficient storage to
  further increase this limit, though the memory capacity can never be better than
  \eqn{O(P^2)}, with \code{P} the number of proteins/genes.
  
  \code{NumCores} only uses 1 less core than is detected, or 1 core if \code{detectCores()}
  cannot detect the number of available cores. This is because of a recurring issue
  on my machine where the R session takes all available cores and is then locked
  out of forking processes, with the only solution to restart the entire R session.
  This may be an issue specific to ARM Macs, but out of an abundance of caution
  I've made the default setting to be slightly slower but guarantee completion
  rather than risk bricking a machine.

  More models will be implemented in the future. Planned models for next release
  include:
  \itemize{
    \item Random Forests for Ensemble predictions
    \item XGBoost for Ensemble predictions
    \item Normalized Phylogenetic Profiles
    \item SVDPhy
    \item DCA at the residue level (Weigt et al. 2009)
  }
  
  Feel free to contact me regarding other models you would like to see added.
}

\seealso{
\code{\link{ProtWeaver}}

\code{\link[=GetProtWebData.ProtWeb]{GetProtWebData}}
}
\examples{
##############
## Prediction with built-in model and data 
###############

exData <- get(data("ExampleStreptomycesData"))
pw <- ProtWeaver(exData$Genes[1:50])

# Subset isn't necessary but is faster for a working example
protweb1 <- predict(pw, Subset=1:10, MySpeciesTree=exData$Tree)

# print out results as an adjacency matrix
GetProtWebData(protweb1)

###############
## Training own ensemble model
###############

datavals <- predict(pw, NoPrediction=TRUE)
                  
actual_values <- sample(c(0,1), nrow(datavals), replace=TRUE)
# This example just picks random numbers
# ***Do not do this for your own models***

# Make sure the actual values correspond to the right pairs! 
datavals[,'y'] <- actual_values                  
myModel <- glm(y~., datavals[,-c(1,2)], family='binomial')

testProtWeaverObject <- ProtWeaver(exData$Genes[51:60])
protweb2 <- predict(testProtWeaverObject, 
                     PretrainedModel=myModel)
                     
# Print result as a 3xN matrix of pairwise scores
GetProtWebData(protweb2, AsDf=TRUE)
}