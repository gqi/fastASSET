# fastASSET: Fast ASSET using pre-screening

fastASSET is an R package for joint genetic association analysis across a large number of traits. It is an accelerated version of the *as*sociation analysis based on sub*set*s (ASSET) method (Bhattacharjee et al, AJHG 2012). The input is GWAS summary statistics for one genetic variant across multiple traits, and the output is a p-value for global association and a list of associated traits. fastASSET accelerates the computation by restricting subset search among the traits with suggestive level of associations, defined by a liberal p-value threshold. 

See the original `ASSET` R package for more information: <https://bioconductor.org/packages/release/bioc/html/ASSET.html>

### Installation

There are two ways to install `fastASSET`:

1. `fastASSET` has been incorporated into the latest version of the main `ASSET` package. First install `devtools` by 

```
install.packages("devtools")
```

and then install the latest version of `ASSET` by

```
devtools::install_github("sbstatgen/ASSET")
```

To load the package and view instructions, type

```
library(ASSET)
?fast_asset
```

2. For the latest developmental version of `fastASSET`, install `devtools` and `ASSET` as in option 1, then install

```
devtools::install_github("gqi/fastASSET")
```

Call the function by `fastASSET::fast_asset()` to access the latest version.

### Reference

If you use this package or other custom scripts from this repository, please cite:

Qi, G., Chhetri, S. B., Ray, D., Dutta, D., Battle, A., Bhattacharjee, S.\*, & Chatterjee, N.\* (2022). Genome-Wide Large-Scale Multi-Trait Analysis Characterizes Global Patterns of Pleiotropy and Unique Trait-Specific Variants. bioRxiv.

*Corresponding authors.

See folders `simulations` and `analysis` for other custom scripts used in the paper.	
