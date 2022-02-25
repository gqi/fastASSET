# Combine information on significant SNPs and variant annotations into one file
rm(list=ls())
library(dplyr)
library(data.table)

load("snps_clumped_all.rda")
load("UKB_sumstats/ukb_search.rda")
load("eqtl/eqtlmat.rda")
load("PATH_TO_DATA/numtraits_avail.rda")
load("PATH_TO_LD_SCORE/l2.rda")
snps_clumped_all$numtraits_avail <- 
    numtraits_avail[match(snps_clumped_all$SNP,names(numtraits_avail))]
snps_clumped_all$ldscore <- l2$L2[match(snps_clumped_all$SNP,l2$SNP)]

hg38 = fread("eqtl/hglft_genome_6c385_8068a0.bed", header = F)
temp = strsplit(hg38$V1, split="-")
snps_clumped_all$hg38 = gsub(":","_",sapply(temp, function(x) x[1]))

eqtlmat[is.na(eqtlmat)] = 0
if (identical(snps_clumped_all$hg38,rownames(eqtlmat))){
    snps_clumped_all$numtissues = rowSums(eqtlmat>0)
    snps_clumped_all$numegenes <- sapply(eqtlegene, function(x) length(unique(unlist(x))))
}

if (identical(snps_clumped_all$SNP,ukb.search$SNP)){
    snps_clumped_combn <- cbind(snps_clumped_all,ukb.search[,-1]) %>%
        select(-c(pheno_pos,pheno_neg)) %>%
        mutate(numtraits_grp = cut(numtraits, breaks=c(0.5,1.5,5.5,10.5,15.5,Inf), 
                                   labels=c("1","2-5","6-10","11-15",">15")))
}
save(snps_clumped_combn, file = "snps_clumped_combn.rda")
