# Welcome!

The `SynExtend` and `DECIPHER` packages for R incorporate a wealth of easy to use functions for comparative genomics analyses. This interactive tutorial series will introduce users to these packages by walking them through a complete workflow of identifying co-evolving genes from a dataset of genome sequences. This webpage was created for presentation at [Bioconductor 2022](https://bioc2022.bioconductor.org/), but the content will be freely available forever.

I've summarized on this page all the skills you can expect to learn by working through the tutorials on this site. See the [Overview](https://www.ahl27.com/CompGenomicsBioc2022/articles/CompGenomicsBioc2022.html) page for more information, along with code examples you can run for yourself! 

(Note: I'm still working on the tutorials, but they will be posted by mid-June! Code is not guaranteed to work until the tutorials are finished--check the [Changelog](https://www.ahl27.com/CompGenomicsBioc2022/news/index.html) for the latest updates).

### [Loading Genome Data with `DECIPHER`](https://www.ahl27.com/CompGenomicsBioc2022/articles/LoadingGenomeData.html)

The first step in analyzing genomics data is loading the data itself. Here we will download sequencing data from NCBI as a `.fasta`, load it into R, then perform some basic operations with the data. Users will learn to efficiently work with large scale genomics data, including visualization and alignment of sequencing data. 

[<sup>Function Reference</sup>](https://www.ahl27.com/CompGenomicsBioc2022/reference/index.html#loading-genome-data)

### [Gene Calling and Annotation with `DECIPHER`](https://www.ahl27.com/CompGenomicsBioc2022/articles/GeneCallingAnnotation.html)

A natural next step is identifying what elements comprise each genome in our dataset. Users will learn to programmatically identify coding and non-coding regions of genomes, and annotate them with predicted KEGG orthology groups using `IDTAXA`. 

[<sup>Function Reference</sup>](https://www.ahl27.com/CompGenomicsBioc2022/reference/index.html#gene-calling-and-annotation)

### [Annotation of COGs with `SynExtend`](https://www.ahl27.com/CompGenomicsBioc2022/articles/FindingCOGs.html)

Annotated genetic regions can be mapped across organisms into clusters of orthologous genes (COGs). Users will learn how to identify COGs at scale using the data generated in the previous step. 

[<sup>Function Reference</sup>](https://www.ahl27.com/CompGenomicsBioc2022/reference/index.html#constructing-cogs)

### [Constructing Gene Trees with `DECIPHER`](https://www.ahl27.com/CompGenomicsBioc2022/articles/ConstructingPhylogenies.html)

Each COG comprises sets of conserved orthologs across species. These data, combined with sequencing data for each ortholog, allow us to reconstruct the evolutionary history of each COG. Users will learn how to construct, visualize, and save phylogenetic trees from sets of genomes using the `TreeLine()` function. 

[<sup>Function Reference</sup>](https://www.ahl27.com/CompGenomicsBioc2022/reference/index.html#constructing-phylogenies)

### [Identifying Co-evolving Gene Collectives with `SynExtend`](https://www.ahl27.com/CompGenomicsBioc2022/articles/CoevolutionNetworks.html)

With these data, we can analyze patterns in evolutionary signal across COGs. Co-evolutionary signal between genes implies functional association, so finding COGs under shared selective pressure aids us in uncovering the mechanisms of intracellular pathways. Users will learn to use the `ProtWeaver` class to tease out subtle evidence of correlated evolutionary pressure in order to create co-evolutionary networks. 

[<sup>Function Reference</sup>](https://www.ahl27.com/CompGenomicsBioc2022/reference/index.html#finding-co-evolving-gene-collectives)


### Conclusion

At the conclusion of this workshop, users will be able to perform the following tasks in R:

* Visualize sequence data
* Identify and annotate genes from sequence data
* Identify COGs from a set of gene calls
* Build phylogenies at the species and gene level
* Predict COGs under shared evolutionary pressure

### Useful Links
* [DECIPHER](http://bioconductor.org/packages/release/bioc/html/DECIPHER.html)
* [SynExtend](http://bioconductor.org/packages/release/bioc/html/SynExtend.html)
* [Related Tutorials](http://www2.decipher.codes/Tutorials.html)
* [Our Lab!](https://www.wrightlabscience.com/p/index.html)


&nbsp;

&nbsp;

&nbsp;

&nbsp;



[![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
