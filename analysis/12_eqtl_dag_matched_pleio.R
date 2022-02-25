##eQTL connection plot (Figure 5)

# Extract eQTL results for trait-specific and highly pleiotropic SNPs
rm(list=ls())
library(dplyr)
library(data.table)

load("PATH_TO_FILES/GTEx_v8_glist.rda")
load("../support_files/snps_clumped_combn.rda")
load("eqtlegene.rda")

# Trait-specific SNPs
loci1t <- snps_clumped_combn %>% filter(numtraits==1 & numtissues>0) %>% 
    mutate(pheno = sapply(pheno, function(x) x[[1]])) %>% 
    arrange(pheno) %>%
    select(-c(traits.ukb,traits.ukb.name))
eqtl1t <- eqtlegene[loci1t$hg38]
names(eqtl1t) <- loci1t$SNP

dag1t <- matrix(NA, nrow=0, ncol=5)
for (snp in names(eqtl1t)){
    print(snp)
    for (ts in names(eqtl1t[[snp]])){
        genes <- eqtl1t[[snp]][[ts]]
        temp <- cbind(rep(snp,nrow(genes)), genes, rep(ts,nrow(genes)))
        dag1t <- rbind(dag1t, temp)
    }
}
colnames(dag1t) <- c("SNP",colnames(genes),"tissue")
dag1t <- data.frame(dag1t, stringsAsFactors = FALSE)

dag1t$gene <- glist$gene_name[match(dag1t$gene_id,glist$gene_id)]
dag1t$SNP <- as.character(dag1t$SNP)
dag1t$tissue <- as.character(dag1t$tissue)
save(dag1t, file = "dag1t.rda")

# Highly pleiotropic SNPs
loci_mt15t <- snps_clumped_combn %>% filter(numtraits>=15 & numtissues>0) %>% 
    mutate(pheno = sapply(pheno, function(x) x[[1]])) %>% 
    arrange(pheno) %>%
    select(-c(traits.ukb,traits.ukb.name))
eqtl_mt15t <- eqtlegene[loci_mt15t$hg38]
names(eqtl_mt15t) <- loci_mt15t$SNP

dag_mt15t <- matrix(NA, nrow=0, ncol=5)
for (snp in names(eqtl_mt15t)){
    print(snp)
    for (ts in names(eqtl_mt15t[[snp]])){
        genes <- eqtl_mt15t[[snp]][[ts]]
        temp <- cbind(rep(snp,nrow(genes)), genes, rep(ts,nrow(genes)))
        dag_mt15t <- rbind(dag_mt15t, temp)
    }
}
colnames(dag_mt15t) <- c("SNP",colnames(genes),"tissue")
dag_mt15t <- data.frame(dag_mt15t, stringsAsFactors = FALSE)

dag_mt15t$gene <- glist$gene_name[match(dag_mt15t$gene_id,glist$gene_id)]
dag_mt15t$SNP <- as.character(dag_mt15t$SNP)
dag_mt15t$tissue <- as.character(dag_mt15t$tissue)
save(dag_mt15t, file = "dag_mt15t.rda")


## Create plot
rm(list=ls())
library(dplyr)
library(data.table)
library(readxl)
library(plotrix)
library(RColorBrewer)
source("eqtl_connection_plot.R")

load("../loci1t_match_specpleio.rda")
top_pleio <- loci1t %>% group_by(pheno_display) %>% summarize(top_pleio_match=SNP_pleio_match[1])
loci1t$SNP_pleio_match <- top_pleio$top_pleio_match[match(loci1t$pheno_display,top_pleio$pheno_display)]

traitvec <- c("alzheimer","BrCa","PrCa","CD","BMDheel","DBP","Height18","mono")
sumstats1t.eqtl <- loci1t %>% filter(numtissues>0) %>%
    mutate(pheno=unlist(pheno), pheno=factor(pheno,levels=traitvec)) %>%
    arrange(pheno)

# Match trait-specific SNPs to pleiotropic SNPs - select the pleiotropic SNP with the strongest association with the trait for specific SNP
load("../../sumstats_clumped_pleio/snps_clumped_sumstats.rda")

# Plot disease traits
load("gtex_colcode.rda") # GTEx color code
coding.genes <- fread("protein-coding_gene.txt")
datalist <- read_excel("../support_files/supptab1_final_data_list_w_DISPLAYNAME.xlsx", sheet=1)

traitvec <- c("alzheimer","BrCa","PrCa","CD")
load("dag1t.rda")
loci1t <- sumstats1t.eqtl %>% filter(pheno %in% traitvec)
dag1t <- dag1t %>% filter(SNP%in%loci1t$SNP) %>%
    mutate(pheno = loci1t$pheno[match(SNP,loci1t$SNP)]) %>%
    filter(gene%in%coding.genes$symbol) %>%
    arrange(pheno) %>%
    mutate(pheno = datalist$Display_name[match(as.character(pheno),datalist$trait)],
           colgrp=pheno)

load("dag_mt15t.rda")
dag_mt15t <- dag_mt15t %>% 
    filter((gene%in%coding.genes$symbol) & (SNP%in%loci1t$SNP_pleio_match)) %>%
    mutate(colgrp=loci1t$pheno[match(SNP,loci1t$SNP_pleio_match)]) %>% # 
    mutate(colgrp=datalist$Display_name[match(colgrp,datalist$trait)]) %>%
    arrange(factor(SNP,levels=unique(loci1t$SNP_pleio_match)))
write.csv(dag_mt15t, file="dag_mt15t_coding.csv")

dag_mt15t$pheno <- NA
for (i in 1:nrow(dag_mt15t)){
    snp <- dag_mt15t$SNP[i]
    dag_mt15t$pheno[i] <- paste0("Matched to ",dag_mt15t$colgrp[i],"\n",
                                 paste(loci1t$SNP[loci1t$SNP_pleio_match==snp],collapse="\n"))
}

dag <- rbind(dag1t,dag_mt15t)

eqtl_connection_plot(dag, gtex.colcode = gtex.colcode, filename="eqtl_matched_disease.png",
                     height = 7000, width=5000, hflen.gene=7, hfht = 0.7)


## Continuous traits
traitvec <- c("BMDheel","DBP","Height18","mono") 
load("dag1t.rda")
loci1t <- sumstats1t.eqtl %>% filter(pheno %in% traitvec)
dag1t <- dag1t %>% filter(SNP%in%loci1t$SNP) %>%
    mutate(pheno = loci1t$pheno[match(SNP,loci1t$SNP)]) %>%
    filter(gene%in%coding.genes$symbol) %>%
    arrange(pheno) %>%
    mutate(pheno = datalist$Display_name[match(as.character(pheno),datalist$trait)],
           colgrp=pheno)

load("dag_mt15t.rda")
dag_mt15t <- dag_mt15t %>% 
    filter((gene%in%coding.genes$symbol) & (SNP%in%loci1t$SNP_pleio_match)) %>%
    mutate(pheno=SNP, colgrp=loci1t$pheno[match(SNP,loci1t$SNP_pleio_match)]) %>% # 
    mutate(colgrp=datalist$Display_name[match(colgrp,datalist$trait)]) %>%
    arrange(factor(SNP,levels=unique(loci1t$SNP_pleio_match)))

dag_mt15t$pheno <- NA
for (i in 1:nrow(dag_mt15t)){
    snp <- dag_mt15t$SNP[i]
    dag_mt15t$pheno[i] <- paste0("Matched to ",dag_mt15t$colgrp[i],"\n",
                                 paste(loci1t$SNP[loci1t$SNP_pleio_match==snp],collapse="\n"))
}

dag <- rbind(dag1t,dag_mt15t) 

eqtl_connection_plot(dag, gtex.colcode = gtex.colcode, filename="eqtl_matched_continuous.png",
                     height = 6500, width=5000, hflen.gene=7, hfht = 0.7)
