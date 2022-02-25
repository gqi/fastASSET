# eQTL lookup of significant SNPs using GTEX v8
## Lift 1000G GRCh37 coordinates to GTEX v8 GRCh38 coordinates
rm(list=ls())
library(dplyr)
load("../support_files/snps_clumped_all.rda")
snps_clumped_all %>% 
    mutate(position=paste0("chr",CHR,":",BP,"-",BP)) %>% 
    select(position) %>% 
    write.table(file="snps_clumped_all_position.txt", row.names=F, col.names=F, quote=F)
# Go to UCSC liftover

## eQTL lookup - original only eqtlmat and egenes
rm(list=ls())
library(data.table)
library(dplyr)
load("../support_files/snps_clumped_all.rda")
hg38 = fread("hglft_genome_6c385_8068a0.bed", header = F)
temp = strsplit(hg38$V1, split="-")
snps_clumped_all$hg38 = gsub(":","_",sapply(temp, function(x) x[1]))

tissues = list.files("PATH_TO_GTEX/GTEx_Analysis_v8_eQTL", pattern=".v8.egenes.txt")
tissues = gsub(".v8.egenes.txt", "", tissues)
eqtlmat = matrix(NA, nrow = nrow(snps_clumped_all), ncol = length(tissues))
rownames(eqtlmat) = snps_clumped_all$hg38
colnames(eqtlmat) = tissues

eqtlegene <- vector("list",length=nrow(snps_clumped_all))
names(eqtlegene) <-  snps_clumped_all$hg38
for (i in 1:length(eqtlegene)){
    eqtlegene[[i]] <- list()
}

for (ts in tissues){
    print(ts)
    dt = fread(paste0("PATH_TO_GTEX/GTEx_Analysis_v8_eQTL/",ts,".v8.signif_variant_gene_pairs.txt"))
    dt = dt %>% mutate(hg38 = gsub("_[A-Z]_[A-Z]_b38","",variant_id)) %>%
        filter(hg38%in%snps_clumped_all$hg38)
    temp = table(dt$hg38)
    eqtlmat[names(temp),ts] = temp
    
    for (i in 1:nrow(snps_clumped_all)){
        temp <- dt %>% filter(hg38==snps_clumped_all$hg38[i])
        if (nrow(temp)>0){
            eqtlegene[[i]][[ts]] <- temp$gene_id
        }
    }
}
save(eqtlmat,eqtlegene, file = "eqtlmat.rda")

## eQTL lookup - with p-values and effect sizes
rm(list=ls())
library(data.table)
library(dplyr)
load("../support_files/snps_clumped_all.rda")
hg38 = fread("hglft_genome_6c385_8068a0.bed", header = F)
temp = strsplit(hg38$V1, split="-")
snps_clumped_all$hg38 = gsub(":","_",sapply(temp, function(x) x[1]))

tissues = list.files("PATH_TO_GTEX/GTEx_Analysis_v8_eQTL", pattern=".v8.egenes.txt")
tissues = gsub(".v8.egenes.txt", "", tissues)
eqtlmat = matrix(NA, nrow = nrow(snps_clumped_all), ncol = length(tissues))
rownames(eqtlmat) = snps_clumped_all$hg38
colnames(eqtlmat) = tissues

eqtlegene <- vector("list",length=nrow(snps_clumped_all))
names(eqtlegene) <-  snps_clumped_all$hg38
for (i in 1:length(eqtlegene)){
    eqtlegene[[i]] <- list()
}

for (ts in tissues){
    print(ts)
    dt = fread(paste0("PATH_TO_GTEX/GTEx_Analysis_v8_eQTL/",ts,".v8.signif_variant_gene_pairs.txt"))
    dt = dt %>% mutate(hg38 = gsub("_[A-Z]_[A-Z]_b38","",variant_id)) %>%
        filter(hg38%in%snps_clumped_all$hg38)
    
    for (i in 1:nrow(snps_clumped_all)){
        temp <- dt %>% filter(hg38==snps_clumped_all$hg38[i])
        if (nrow(temp)>0){
            eqtlegene[[i]][[ts]] <- temp %>% select(gene_id, pval_nominal, slope)
        }
    }
}
save(eqtlegene, file = "eqtlegene.rda")

