# Download Neale lab UK Biobank GWAS summary statistics for validation of pleiotropy
rm(list=ls())
library(data.table)
library(dplyr)

# Loading list of Lead SNPs for 2,293 significant loci identified by ASSET
load("../support_files/snps_clumped_all.rda")
# Lead SNPs of locus that are not present in Neale lab UKB GWAS, and their proxy SNPs
load("../support_files/asset_noukb.rda")

snps_clumped_all$BP_proxy = snps_clumped_all$BP
snps_clumped_all$BP_proxy[match(asset_noukb$SNP,snps_clumped_all$SNP)] = asset_noukb$BP_proxy
snps_clumped_all = snps_clumped_all %>% filter(!is.na(BP_proxy)) %>%
    select(CHR,BP,SNP,numtraits,BP_proxy) %>%
    mutate(variant_id_proxy = paste0(CHR,":",BP_proxy))

nealeukb <- fread("../support_files/Neale_UKB_traitlist.csv")
# Remove non-phenotype files, gender and age
# Choose GWASs that pool both sexes, and use inverse normal transformed phenotypes
nealeukb <- nealeukb %>% filter(!(`Phenotype Code`%in%c("N/A","is_female","age")) &
                                    Sex=="both_sexes" & !grepl("_raw",`Phenotype Code`))
grp <- as.integer(commandArgs(trailingOnly=TRUE)) # 1:107; 4280 datasets in total, 40 datasets per job
nealeukb.sub <- nealeukb[(40*(grp-1)+1):(40*grp),]

## Download and read UKB data
for (i in 1:nrow(nealeukb.sub)){
    pheno <- nealeukb.sub$`Phenotype Code`[i]
    filename = gsub("\\.bgz","",nealeukb.sub$File[i])
    # paste0(pheno,".gwas.imputed_v3.both_sexes.tsv")
    print(paste((40*(grp-1)+i),pheno))
    
    system(nealeukb.sub$`wget command`[i])
    system(paste0("mv ",filename,".bgz ", filename,".gz"))
    system(paste0("gunzip ",filename,".gz"))
    
    ukb_sumstats = fread(filename, nThread = 1)
    ukb_sumstats <- ukb_sumstats %>% mutate(variant_id = gsub(":[A-Z]:[A-Z]","",variant)) %>%
        select(variant, variant_id, n_complete_samples, beta, se)
    colnames(ukb_sumstats)[c(3,4,5)] <- paste0(colnames(ukb_sumstats)[c(3,4,5)],".",pheno)
    if (i==1){
        sumstats <- snps_clumped_all %>% left_join(ukb_sumstats, by=c("variant_id_proxy"="variant_id"))
    } else{
        sumstats <- sumstats %>% 
            left_join(select(ukb_sumstats,-variant), by=c("variant_id_proxy"="variant_id"))
    }
    
    save(sumstats, file = paste0("ukb_sumstats_grp",grp,".rda"))
    system(paste0("rm ",filename))
}
save(sumstats, file = paste0("ukb_sumstats_grp",grp,".rda"))
