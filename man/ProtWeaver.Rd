\name{ProtWeaver}
\alias{ProtWeaver}
\alias{ProtWeaver-class}
\title{
ProtWeaver: Predicting Protein Functional Association Networks
}
\description{
ProtWeaver is an S3 class with methods for predicting functional association using
protein or gene data. ProtWeaver implements several methods utilized in the literature, with many more planned for future implementation. For details on predictions, see \code{\link{predict.ProtWeaver}}.
}
\usage{
ProtWeaver(ListOfData, NoWarn=FALSE)
}
\arguments{
  \item{ListOfData}{
A list of gene data, where each entry corresponds to information on a particular
gene. List must contain either dendrograms or vectors, and cannot contain a mixture. 
If list is composed of dendrograms, each dendrogram is a gene tree for the corresponding entry. If list is composed of vectors, vectors should be numeric or character vectors denoting the genomes containing that gene. 
}
\item{NoWarn}{
Several algorithms depend on having certain data. When a \code{ProtWeaver} object is initialized, it automatically selects which algorithms can be used given the input data. By default, ProtWeaver will notify the user of algorithms that cannot be used with warnings. Setting \code{NoWarn=TRUE} will suppress these messages.
}
}
\details{
ProtWeaver expects input data to be a list. All entries must be one of the following:
\enumerate{
  \item{
    \code{ListOfData[[i]] = c('ID#1', 'ID#2', ..., 'ID#k')}
  }
  \item{
  \code{ListOfData[[i]] = c('i1_d1_p1', 'i2_d2_p2', ..., 'ik_dk_pk') }
  }
  \item{
  \code{ListOfData[[i]] = dendrogram(...)}
  }
}
In (1), each ID#i corresponds to the unique identifier for genome #i. For entry #j in the list, the presence of 'ID#i' means genome #i has an ortholog for gene/protein #j. 

Case (2) is the same as (1), just with the formatting of names slightly different. Each entry is of the form \code{i_d_p}, where \code{i} is the unique identifier for the genome, \code{d} is which chromosome the ortholog is located, and \code{p} is what position the ortholog appears in on that chromosome. \code{p} must be a \code{numeric}, while the other entries can be any value.

Case (3) expects gene trees for each gene, with labeled leaves corresponding to each source genome. If \code{ListOfData} is in this format, taking \code{labels(ListOfData[[i]])} should produce a character vector that matches the format of either (2) or (1).

\emph{See the Examples section for illustrative examples.}

ProtWeaver requires input of scenario (3) to use MirrorTree or ContextTree, and requires input of scenario (2) (or (3) with leaves labeled according to (2)) for co-localization analyses.

Note that ALL entries must belong to the same category--a combination of character vectors and dendrograms is not allowed. 

Prediction of a functional association network is done using \code{predict(ProtWeaverObject)}. See \code{\link{predict.ProtWeaver}} for more information.
}
\value{
Returns a ProtWeaver object.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}
\seealso{
\code{\link{predict.ProtWeaver}},
\code{\link{ExampleStreptomycesData}},
\code{\link{BuiltInEnsembles}}
}
\examples{
# I'm using gene to mean either a gene or protein

## Imagine we have the following 4 genomes:
## (each letter denotes a distinct gene)
##    Genome 1: a b c d
##    Genome 2: d c e
##    Genome 3: b a e 
##    Genome 4: a e

## We have 5 total genes: (a,b,c,d,e)
##    a is present in genomes 1, 3, 4
##    b is present in genomes 1, 3
##    c is present in genomes 1, 2
##    d is present in genomes 1, 2
##    e is present in genomes 2, 3, 4

## Constructing a ProtWeaver object according to (1):
l <- list()
l[['a']] <- c('1', '3', '4') 
l[['b']] <- c('1', '3') 
l[['c']] <- c('1', '2') 
l[['d']] <- c('1', '2') 
l[['e']] <- c('2', '3', '4') 

## Each value of the list corresponds to a gene
## The associated vector shows which genomes have that gene
pwCase1 <- ProtWeaver(l)

## Constructing a ProtWeaver object according to (2):
##  Here we need to add in the chromosome and the position
##  As we only have one chromosome, 
##  we can just set that to 1 for all.
##  Position can be identified with knowledge, or with
##  FindGenes(...) from DECIPHER.

## In this toy case, genomes are small so it's simple.
l <- list()
l[['a']] <- c('1_1_1', '3_1_2', '4_1_1') 
l[['b']] <- c('1_1_2', '3_1_1') 
l[['c']] <- c('1_1_3', '2_1_2') 
l[['d']] <- c('1_1_4', '2_1_1') 
l[['e']] <- c('2_1_3', '3_1_3', '4_1_2') 

pwCase2 <- ProtWeaver(l)


## For Case 3, we just need dendrogram objects for each
# l[['a']] <- dendrogram(...)
# l[['b']] <- dendrogram(...)
# l[['c']] <- dendrogram(...)
# l[['d']] <- dendrogram(...)
# l[['e']] <- dendrogram(...)

## Leaf labels for these will be the same as the 
##  entries in Case 1.
}
