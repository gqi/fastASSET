# fastASSET: Fast ASSET using pre-screening

This package implements a fast version of ASSET for analysis of a large number of traits. It accelerates the computation by restricting subset search among the traits with suggestive level of associations. The traits are selected by a liberal p-value threshold. The input is GWAS summary statistics for one SNP across multiple traits.

### Installation

There are two ways to install fast ASSET.

1. Fast ASSET has been incorporated into the latest version of the main ASSET package. First install `devtools` by 

```
install.packages("devtools")
```

and then install the latest version of ASSET by

```
devtools::install_github("sbstatgen/ASSET")
```

Type `?fast_asset` for help file and example.

2. For the latest developmental version of fast ASSET, install `devtools` and `ASSET` as in option 1, then install

```
devtools::install_github("gqi/fastASSET")
```
	
