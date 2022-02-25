# Match trait-specific SNPs to highly pleiotropic SNPS (>15) that has the strongest effect for the trait
# This enables direct comparison between the functional mechanisms of trait-specific and highly pleiotropic SNPs
## Match for eQTL plot - Figure 5
rm(list=ls())
library(dplyr)
library(data.table)
library(readxl)
library(ggplot2)
library(ggthemes)
library(gridExtra)

load("../support_files/snps_clumped_combn.rda")
sumstats1t <- snps_clumped_combn %>% filter(numtraits==1)
sumstats.mt15t <- snps_clumped_combn %>% filter(numtraits>15)

traitvec <- c("alzheimer","BrCa","PrCa","CD","BMDheel","DBP","Height18","mono")
sumstats1t.eqtl <- sumstats1t %>% filter(numtissues>0) %>%
    mutate(pheno=unlist(pheno))

## Match trait-specific SNPs to the pleiotropic SNP with the strongest trait association
load("PATH_TO_FILE/snps_clumped_sumstats.rda")
match.specpleio <- data.frame(matrix(NA,nrow=0,ncol=5,dimnames=list(NULL,c("pheno.spec","SNP","numtraits","numtissues","numegenes"))))
for (i in 1:length(traitvec)){
    trait <- traitvec[i]
    z.trait.all <- snps_clumped_sumstats[[paste0("beta_",trait)]]/
        snps_clumped_sumstats[[paste0("se_",trait)]]
    temp <- sumstats.mt15t %>% filter(numtissues>0) %>%
        filter(sapply(pheno, function(x) trait%in%x) ) %>%  
        mutate(z.trait=z.trait.all[match(SNP,snps_clumped_sumstats$SNP)]) %>%
        slice_max(abs(z.trait), n=1) %>% 
        select(SNP,numtraits,numtraits.ukb,numtissues,numegenes) %>%
        mutate(pheno.spec=trait) %>% relocate(pheno.spec)
    
    match.specpleio <- rbind(match.specpleio, temp)
}
save(match.specpleio, file="match_specpleio.rda")

## Match SNPs for Figure 4
rm(list=ls())
library(dplyr)
library(data.table)
library(readxl)

# Load lists of trait-specific and pleiotropic SNPs
trait.list.all <- read_excel("../support_files/supptab1_final_data_list_w_DISPLAYNAME.xlsx", sheet=1)

pheno.order <- c("alzheimer","BrCa","PrCa","CD","BMDheel","DBP","Height18",
                 "intelligence","male_baldness","menarche","mono")
loci1t <- read.csv("loci_1trait.csv", header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(pheno = factor(pheno, levels=pheno.order)) %>%
    arrange(pheno) %>%
    mutate(pheno_display = trait.list.all$Display_name[match(pheno,trait.list.all$trait)])

load("../support_files/snps_clumped_combn.rda")
loci.mt15t <- snps_clumped_combn %>% filter(numtraits>15)

load("PATH_TO_FILE/snps_clumped_sumstats.rda")
match.specpleio <- data.frame(matrix(NA,nrow=0,ncol=3,dimnames=list(NULL,c("pheno.spec","SNPpleio","SNPspec"))))
for (i in 1:nrow(loci1t)){
    trait <- as.character(loci1t$pheno[i]) # 
    z.trait.all <- snps_clumped_sumstats[[paste0("beta_",trait)]]/
        snps_clumped_sumstats[[paste0("se_",trait)]]
    temp <- loci.mt15t %>%
        filter(sapply(pheno, function(x) trait%in%x) & !(SNP%in%match.specpleio$SNPpleio)) %>%  
        mutate(z.trait=z.trait.all[match(SNP,snps_clumped_sumstats$SNP)],
               SNPspec=loci1t$SNP[i]) %>%
        slice_max(abs(z.trait), n=1) %>% # 
        select(SNPpleio=SNP,SNPspec) %>% mutate() %>%
        mutate(pheno.spec=trait) %>% relocate(pheno.spec)
    
    match.specpleio <- rbind(match.specpleio, temp)
}

loci1t$SNP_pleio_match <- match.specpleio$SNPpleio
save(loci1t, file="loci1t_match_specpleio.rda")
