# fastASSET: Fast ASSET using pre-screening

This R package implements a fast version of ASSET for joint genetic association analysis across a large number of traits. It is developed based on ASSET, a subset-based approach for multi-trait association testing (Bhattacharjee et al, AJHG 2012). It provides a p-value of global association and a list of traits associated with each variant. fastASSET accelerates the computation by restricting subset search among the traits with suggestive level of associations. The traits are selected by a liberal p-value threshold. The input is GWAS summary statistics for one SNP across multiple traits.

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

Qi G, Chhetri SB, Ray D, Dutta D, Battle A, Bhattacharjee S*, Chatterjee N*. Genome-wide multi-trait analysis across 116 studies identifies complex patterns of pleiotropy and unique trait-specific loci. (2022).

*Corresponding authors.

Other scripts for simulations studies and data analysis used in this paper are in folders `simulations` and `analysis`.	
