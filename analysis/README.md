# Real data analysis

`01_ASSET_GRASP_analysis.R`: Run fastASSET analysis on GRASP and GWAS consortium data. Instead of using the `fast_asset` wrapper function, this script is an older version that uses indvidual component functions of `fast_asset`. Main component functions are `scr.transform.block` and `h.traits`.

`01_process_ASSET_results.R`: Process fastASSET results and remove SNPs with convergence issues.

`02_plot_gcorr_ldscint_hclust.R`: Plot LDSC intercept matrix (`ldscintmat`) and hierarchical clustering results.

`03_brisbane_plot.R`: Create brisbane plot (Figure 3).

`04_LD_clumping.R`: Perform LD clumping based on fastASSET global association p-values.

`05_pleio_snps_excluded.R`: Analyze the SNPs excluded from MR analysis due to large number of traits passing pre-screening.

`06_download_UKB.R`: Download UK Biobank summary statistics for validation of estimated pleiotropy.

`07_UKB_search.R`: Search UK Biobank summary statistics and extract results for significant SNPs of fastASSET analysis.

`08_eqtl.R`: Extract multi-tissue eQTL association statistics from GTEx v8 for significant SNPs of fastASSET analysis.

`09_combine_clumpedSNPs_info.R`: Combine variant annotation for top SNPs into a single table.

`10_numtraits_figure.R`: Plot the distribution of estimated pleiotropy and relationship with variant annotations (Figure 4).

`11_match_traitspec_pleio_SNPs.R`: Match each trait-specific SNP to a highly pleiotropic SNP (associaed with >15 traits) that is also associated with the trait. If there are multiple highly pleiotropic SNPs, choose the one with the most significant association with the trait.

`12_eqtl_dag_matched_pleio.R`: Create eQTL connection plot (Figure 6) for trait-specific SNPs and matched highly pleiotropic SNPs.

`13_phewas_lineplt_matched.R`: Create the plot of p-values for top 10 associated traits for trait-specific and highly pleiotropic SNPs (Figure 5).

`14_chromstate_IDEAS_matched.R`: Identify cell types where selected SNPs are in active chromatin states based on Roadmap Epigenomics data (Figure 7).
