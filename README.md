# ScSpAlleleExperiment

Defines a S4 class that is based on SingleCellExperiment, but in addition      
to the usual gene layer can also store data for immune genes such as HLAs, Igs and KIRs at allele and functional level.

## Biological background and motivation

Immune molecules such as B and T cell receptors, human leukocyte antigens (HLAs) or killer Ig-like
receptors (KIRs) are encoded in the genetically most diverse loci of the human genome. Many of
these immune genes are hyperpolymorphic, showing high allelic diversity across human populations.
In addition, typical immune molecules are polygenic, which means that multiple functionally similar
genes encode the same protein subunit.

We have developed a workflow, that allows quantification and interactive exploration of 
donor-specific alleles of different immune genes. The workflow is composed of three steps: 
1. Typing of donor-specific alleles, 
2. Quantification of these alleles in single-cell transcriptomic data, 
3. exploration of single-cell transcriptomic data using `SingleCellAlleleExperiment (SCAE)` class. The software for steps 1 and 2 is currently not yet publicly available.

The aim of the here presented `SingleCellAlleleExperiment (SCAE)` 
class is to provide an addition to the `SingleCellExperiment` object that
can simultaneously contain quantification of alleles, genes, and groups of functionally similar genes
and thus allows data analysis across these immunologically relevant different layers of annotation. 

## The `SingleCellAlleleExperiment (SCAE)` class

The `SingleCellAlleleExperiment (SCAE)` class is a container for storing and handling allele-aware quantification data for immune genes. The SCAE class is derived
from the `r Biocpkg("SingleCellExperiment")``(SCE)` class and uses the same overall object architecture. However, multiple data layers are integrated into the object during object generation.
Data from a novel allele-aware quantification method contains information about alleles for immune genes. During object generation, the allele information is aggregated into two additional data layers
via an ontology based MHC design principle and appended to the initial raw data using a lookup table. Thus the final SCAE object contains non-immune genes, alleles for immune genes, an immune gene layer and a functional class layer (refer to Figure 1). 

For example, the counts of the alleles `A*01:01:01:01` and `A*02:01:01:01` that are present in the raw input data will be summed up and transformed into the `HLA-A` immune gene layer. Next, all counts for the present HLA-class I immune genes will be summed up and transformed into the `HLA-class I` functional class layer. The information necessary to perform these transformations is saved in a lookup table retrieved from the IPD-IMGT/HLA database. 

The implemented object follows similar conventions like the SCE class, where rows should represent features (genes, transcripts) and columns should represent cells. Established single cell packages like `r Biocpkg("scater")` and `r Biocpkg("scran")` can be used with the SCAE object to perform downstream analysis on immune gene expression. 

This allows new insights on functional as well as allele level to uncover the high diversity of immune genes.

![alt text here](./inst/extdata/scae_advanced.png)
**Figure 1:** Scheme of SingleCellAlleleExperiment object structure with lookup table.