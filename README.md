# fastASSET: Fast ASSET using pre-screening

fastASSET is an R package for joint genetic association analysis across a large number of traits. It is an accelerated version of the *as*sociation analysis based on sub*set*s (ASSET) method (Bhattacharjee et al, AJHG 2012). The input is GWAS summary statistics for one genetic variant across multiple traits, and the output is a p-value for global association and a list of associated traits. fastASSET accelerates the computation by restricting subset search among the traits with suggestive level of associations, defined by a liberal p-value threshold. 

See the original `ASSET` R package for more information: <https://bioconductor.org/packages/release/bioc/html/ASSET.html>

### Installation

First install `devtools` and `ASSET`

```
install.packages("devtools")
devtools::install_github("sbstatgen/ASSET")
```

then install `fastASSET`

```
devtools::install_github("gqi/fastASSET")
```

Type `?fastASSET::fast_asset()` to view the documentation.


### Step-by-step tutorial

#### 1. Load example dataset.
```
library(fastASSET)
data("example_rs6678982", package="fastASSET")
```

Example dataset `example_rs6678982` has been provided as part of the fastASSET package. It includes data of a single SNP required to run fastASSET :

* `SNP`: SNP name
* `traits`: Vector of trait names.
* `betahat`: Effect of SNP on traits, obtained from GWAS summary statistics.
* `SE`: Standard error of `betahat`.
* `Neff`: Vector of effective sample size of GWAS. For continuous traits, the effective sample size is the total sample size; for binary traits, the effective sample size is `Ncase*Ncontrol/(Ncase+Ncontrol)`.
* `ldscintmat`: Matrix of bivariate LD score regression intercepts. It estimates the correlation of z-statistics across traits under the global null hypothesis.

We have provided `ldscintmat` with this example. In real data analysis, `ldscintmat` can be obtained by running bivariate [LD score regression](https://github.com/bulik/ldsc) for each pair of traits. See the [LDSC tutorial](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) for running bivaraite LD score regression.

#### 2. Create correlated trait blocks using hierarchical clustering

This step partitions the traits into smaller blocks correlated traits. Traits from independent blocks are treated as independent in subsequent analysis.
```
block <- create_blocks(ldscintmat)
```

#### 3. Conduct fastASSET analysis

```
test <- fast_asset(snp=SNP, traits.lab=traits, beta.hat=betahat, sigma.hat=SE,
Neff=Neff, cor=ldscintmat, block=block, scr_pthr=0.05)
```



### Reference

If you use this package or other custom scripts from this repository, please cite:

Qi, G., Chhetri, S. B., Ray, D., Dutta, D., Battle, A., Bhattacharjee, S.\*, & Chatterjee, N.\* (2022). Genome-Wide Large-Scale Multi-Trait Analysis Characterizes Global Patterns of Pleiotropy and Unique Trait-Specific Variants. bioRxiv. In press for Nature Communications.

*Corresponding authors.

See folders `simulations` and `analysis` for other custom scripts used in the paper.	
