
<p align="center">
<img src="man/figures/logo.png" height="150" />
</p>

**XploR** is an R package for robust, large-scale copy number and
allelic imbalance analysis from whole exome sequencing (WES) data.  
It provides tools for segmentation, purity/ploidy estimation, allelic
imbalance analysis, ISCN annotation, interactive plotting, and more.

------------------------------------------------------------------------

## Features

- Exome-wide copy number segmentation and allelic imbalance detection
- Purity and ploidy estimation with model selection
- BAF and coverage smoothing, binning, and quality control
- ISCN and gene annotation of CNV segments
- Interactive and publication-ready visualization
- Batch and command-line workflows for high-throughput analysis

------------------------------------------------------------------------

## Installation

Install the latest version from GitHub using
[devtools](https://github.com/r-lib/devtools):

``` r
install.packages("devtools")
devtools::install_github("lingxi1002/XploR")
```

## Example

This is a basic example to run the program:

``` r
library(XploR)

#Example: Run Allelic inbalance bases segmentation (see vignette for details)
result <- runAIsegmentation(
 seg = "sample.seg", 
  cov = "sample.counts",
  ai = "sample.tumor.baf.bedgraph.gz",
  gender = "female",
  out_dir = "results/",
  prefix = "Sample1"
)

#Example: Run CNV calling (see vignette for details)
runModelLikelihood(
    seg = opt$seg,
    out_dir = opt$out_dir,
    prefix = opt$prefix,
    gender = opt$gender,
    lambda = 1,
    gamma = 1,
    epsilon = 0.01,
    modelminprobes = 20,
    modelminAIsize = 5000000,
    minsf = 0.4,
    callcov = 0.3,
    thread = 4,
    diploidweight = 0.5
)


#Annotate segments with ISCN and gene information
AnnotateSegments(
  input = "results/Sample1_GATK_DRAGEN_merge.tsv",
  out_dir = "results/",
  prefix = "Sample1",
  cytoband = "data/cytoBand.txt",
  whitelist_edge = "data/whitelist.txt",
  gene = "data/gene_anno.txt"
)

# Plot results
plot_cnv_profile(
  seg = "results/Sample1_CNV_annotation.tsv",
  coverage = "sample.denoisedCR.tsv",
  ballele = "sample.tumor.ballele.counts.gz",
  whitelist = "data/whitelist.txt",
  gender = "female",
  out_dir = "results/",
  prefix = "Sample1"
)
```
