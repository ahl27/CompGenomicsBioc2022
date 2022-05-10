# Workshop Overview

The `SynExtend` and `DECIPHER` packages for R incorporate a wealth of easy to use functions for comparative genomics analyses. This workshop will introduce users to these packages by walking them through a complete workflow of identifying co-evolving genes from a dataset of genome sequences. 

This workflow comprises several steps, each of which are detailed below. More detail will be added later; this is currently a work in progress.

## Loading Genome Data with `DECIPHER`

The first step in analyzing genomics data is loading the data itself. Here we will download sequencing data from NCBI as a `.fasta`, load it into R, then visualize and align the sequences.

## Gene Calling and Annotation with `DECIPHER`

Once we have the data, we will identify genetic regions with `FindGenes()` from `DECIPHER`, then annotate them with `IdTaxa()`.

## Annotation of COGs with `SynExtend`

Now that we have annotated gene calls, we can combine them into clusters of orthologous genes (COGs) with `DisjointSets()` from `SynExtend`.

## Constructing Gene Trees with `DECIPHER`

We will then build phylogenetic trees from each COG using the new `TreeLine()` function from `DECIPHER`.

## Identifying Co-evolving Gene Collectives with `SynExtend`

With COGs and gene trees in hand, we can identify gene clusters under shared evolutionary pressure with the `ProtWeaver` class in `SynExtend`. Co-evolutionary signal implies functional association, so the resulting pairwise associations are useful for finding functionally associated genes/proteins.

## Conclusion

At the conclusion of this workshop, users will be able to perform the following tasks in R:
* Visualize sequence data
* Identify and annotate genes from sequence data
* Identify COGs from a set of gene calls
* Build phylogenies at the species and gene level
* Predict COGs under shared evolutionary pressure

## Useful Links
* [DECIPHER](http://bioconductor.org/packages/release/bioc/html/DECIPHER.html)
* [SynExtend](http://bioconductor.org/packages/release/bioc/html/SynExtend.html)
* [Related Tutorials](http://www2.decipher.codes/Tutorials.html)
* [Our Lab!](https://www.wrightlabscience.com/p/index.html)
