## Validating pleiotropy using UK Biobank data

# Merge UKB sumstats
rm(list=ls())
library(data.table)
library(dplyr)

for (i in 1:107){
    if (i%%10==0) print(i)
    load(paste0("ukb_sumstats_grp",i,".rda"))
    beta <- as.matrix(sumstats[,seq(9,ncol(sumstats),by=3)])
    se <- as.matrix(sumstats[,seq(10,ncol(sumstats),by=3)])
    temp <- 2*pnorm(-abs(beta/se))
    colnames(temp) <- gsub("^beta","pval",colnames(beta))
    if (i==1){
        sumstats.all <- cbind(sumstats[,1:7],temp)
    } else{
        if (all.equal(sumstats.all$SNP,sumstats$SNP)){
            sumstats.all <- cbind(sumstats.all,temp)
        }
    }
}
save(sumstats.all, file = "sumstats_all.rda")

for (i in 1:107){
    if (i%%10==0) print(i)
    filename <- paste0("ukb_sumstats_grp",i,".rda")
    load(filename)
    if (ncol(sumstats)<127){
        print(filename)
        print(dim(sumstats))
    }
}


# Search UKB
rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)

nealeukb <- fread("../support_files/Neale_UKB_traitlist.csv")
# Remove non-phenotype files, gender and age
nealeukb <- nealeukb %>% filter(!(`Phenotype Code`%in%c("N/A","is_female","age")) &
                                    Sex=="both_sexes" & !grepl("_raw",`Phenotype Code`))
phenorep <- table(nealeukb$`Phenotype Code`)
phenorep <- names(phenorep)[phenorep>1]
# If the trait is repeated, remove the second version, keep the original version
nealeukb <- nealeukb %>% mutate(repeated = `Phenotype Code`%in%phenorep & grepl("\\.v2\\.",File))

load("../support_files/snps_clumped_all.rda")
load("sumstats_all.rda")
all.equal(nealeukb$`Phenotype Code`, gsub("pval.","",colnames(sumstats.all)[-(1:7)]))

sumstats.all <- sumstats.all[,-c(7+which(nealeukb$repeated))]
temp <- as.matrix(sumstats.all[,8:ncol(sumstats.all)])
ukb.traits <- apply(temp, 1, function(x) {y=x[!is.na(x)]; gsub("pval.","",colnames(temp)[p.adjust(y,method="BH")<0.05])})

sumstats.all$traits.ukb <- ukb.traits
sumstats.all$numtraits.ukb <- sapply(ukb.traits,length)

snps_clumped_all <- snps_clumped_all %>% 
    left_join(sumstats.all[,c("SNP","numtraits.ukb","traits.ukb")],by="SNP") %>%
    select(SNP,numtraits.ukb,traits.ukb)

nealeukb <- nealeukb %>% filter(!repeated)
ukb.search <- snps_clumped_all %>% 
    mutate(traits.ukb.name = lapply(traits.ukb, function(x) nealeukb$`Phenotype Description`[match(x,nealeukb$`Phenotype Code`)]))
save(ukb.search, file = "ukb_search.rda")

