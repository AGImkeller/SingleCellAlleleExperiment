# SingleCellAlleleExperiment

Defines a S4 class that is based on `SingleCellExperiment`. In addition to the usual gene layer `SingleCellAlleleExperiment` can also store data for immune genes such as HLAs, Immunoglobulins and KIRs at the allele level and at the level of functionally similar groups of immune genes.

## Biological background and motivation

Immune molecules such as B and T cell receptors, human leukocyte antigens (HLAs) or killer Ig-like
receptors (KIRs) are encoded in the genetically most diverse loci of the human genome. Many of
these immune genes are hyperpolymorphic, showing high allelic diversity across human populations.
In addition, typical immune molecules are polygenic, which means that multiple functionally similar
genes encode the same protein subunit.

We have developed a workflow, that allows quantification of expression and interactive exploration of 
donor-specific alleles of different immune genes. The workflow is composed of three steps: 
1. typing of donor-specific alleles, 
2. quantification of these alleles in single-cell transcriptomic data, 
3. exploration of single-cell transcriptomic data using `SingleCellAlleleExperiment (SCAE)` class. 

The software for steps 1 and 2 can be found in the form of a *Snakemake* workflow named *[scIGD](https://github.com/AGImkeller/scIGD)*. For detailed instructions on utilizing this workflow, please refer to its documentation.

The present repository is dedicated to the implementation of the `SingleCellAlleleExperiment (SCAE)` class. The aim of the here presented class is to provide an addition to the `SingleCellExperiment (SCE)` object that can simultaneously contain quantification of alleles, genes, and groups of functionally similar genes and thus allows data analysis across these immunologically relevant different layers of annotation. 

## The `SingleCellAlleleExperiment (SCAE)` class

The `SingleCellAlleleExperiment (SCAE)` class is a container for storing and handling allele-aware quantification data for immune genes. The SCAE class is derived from the `SingleCellExperiment (SCE)` class and uses the same overall object architecture. However, multiple data layers are integrated into the object during object generation. During object generation, the allele information is aggregated into two additional data layers via an ontology based design principle and appended to the initial raw data using a lookup table. Thus the final `SCAE` object contains quantification of all classical genes ("non-immune genes") and additionally a multi-layer representation of a set of genes of interest (e.g. "immune genes" such as HLAs, Igs, KIRs). For these genes of interest, the quantification is stored on the levels of alleles, genes, and functionally similar groups of genes (refer to Figure 1). 

For example, the counts of the alleles `A*01:01:01:01` and `A*02:01:01:01` that are present in the raw input data will be combined into the `HLA-A` immune gene layer. Next, all counts for the present HLA-class I immune genes will be combined into the `HLA-class I` functional class layer. The information necessary to perform these transformations is saved in a lookup table retrieved from the IPD-IMGT/HLA database and other immunogenetic databases. 

The implemented object follows similar conventions like the `SCE` class, where rows should represent features (genes, transcripts) and columns should represent cells. Established single cell packages like `r Biocpkg("scater")` and `r Biocpkg("scran")` can be used with the newly implemented `SCAE` object to perform downstream analysis on immune gene expression. This allows data exploration on functional as well as allele level.

![alt text here](./inst/extdata/scae_advanced.png)

**Figure 1:** Scheme of SingleCellAlleleExperiment object structure with lookup table.

## Installation

`SingleCellAlleleExperiment` and its data package `scaeData` are available in Bioconductor and can be installed as follows:

```markdown
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("scaeData")
BiocManager::install("SingleCellAlleleExperiment")
```

Alternatively, they can be installed from GitHub using the [devtools](https://github.com/r-lib/devtools) package:

```markdown
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("AGImkeller/scaeData", build_vignettes = TRUE)
devtools::install_github("AGImkeller/SingleCellAlleleExperiment", build_vignettes = TRUE)
```

## Citation

To be added..

## Authors 

- [Jonas Schuck](https://github.com/Jonas-Schuck), [Ahmad Al Ajami](https://github.com/ahmadalajami), [Federico Marini](https://github.com/federicomarini), [Katharina Imkeller](https://github.com/imkeller)