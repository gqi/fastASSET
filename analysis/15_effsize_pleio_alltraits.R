## Create the effect size ranking vs pleiotropy plot (Supplementary Figure 10)

## Extract data
rm(list=ls())
library(data.table)
library(dplyr)

load("../res_all.rda")
res_all <- res_all %>% filter(error_warning!=6)

for (i in 1:18){
    print(i)
    load(paste0("PATH_TO_FILE/data_split",i,".rda"))
    if (i==1){
        zstat_all <- data_split$zstat
        sampsize_all <- data_split$sampsize
    } else{
        zstat_all <- rbind(zstat_all, data_split$zstat)
        sampsize_all <- rbind(sampsize_all, data_split$sampsize)
    }
}

# Only keep variants that reach genome-wide significance with at least one trait
ind <- rowSums(abs(zstat_all)>qnorm(2.5e-8,lower.tail=FALSE),na.rm=T)>=1 &
    rownames(zstat_all) %in% res_all$SNP
zstat_sig <- zstat_all[ind,]
sampsize_sig <- sampsize_all[ind,]
save(zstat_sig, sampsize_sig, file="data_sig.rda")


### LD clumping
rm(list=ls())
library(data.table)
library(dplyr)
source("clump.R")
load("PATH_TO_1000G/bim.rda")
load("data_sig.rda")

temp <- as.integer(commandArgs(trailingOnly = TRUE)) # 1:12

pthr <- 5e-8
for (trait in colnames(zstat_sig)[(10*(temp-1)+1):min(10*temp,ncol(zstat_sig))]){
    zstat <- zstat_sig[,trait]
    zstat <- zstat[!is.na(zstat)]
    
    ind <- match(names(zstat), bim$SNP)
    assoc <- data.frame(SNP = names(zstat), CHR = bim$CHR[ind],
                        BP = bim$BP[ind], P = 2*pnorm(-abs(zstat)))
    snps.clumped <- clump(assoc, pthr = pthr, r2thr = 0.1, bpthr = 5e5,
                          plink.path = "PATH_TO_PLINK",
                          ref.geno = paste0("PATH_TO_1000G/1000G.EUR.QC.",1:22),
                          ld.out = paste0("out_",trait,"_"))
    
    save(snps.clumped, file = paste0("snps_clumped_",trait,"_pthr",pthr,".rda"))
}

## Combine datasets
rm(list=ls())
library(data.table)
library(dplyr)
load("data_sig.rda")

pthrvec <- c(1e-4,1e-6,5e-8)
numtraits.hardthr <- matrix(NA, nrow=nrow(zstat_sig), ncol=length(pthrvec))
rownames(numtraits.hardthr) <- rownames(zstat_sig)
colnames(numtraits.hardthr) <- paste("pthr",pthrvec)
for (i in 1:length(pthrvec)){
    print(i)
    numtraits.hardthr[,i] <- rowSums(abs(zstat_sig)>qnorm(pthrvec[i]/2, lower.tail=F), na.rm=T)
}

beta <- zstat_sig/sqrt(sampsize_sig)

cor.beta.numtraits <- matrix(NA, nrow=ncol(zstat_sig), ncol=length(pthrvec))
rownames(cor.beta.numtraits) <- colnames(zstat_sig)
colnames(cor.beta.numtraits) <- paste("pthr",pthrvec)
pval.beta.numtraits <- cor.beta.numtraits

beta_numtraits_clump <- data.frame()
for (i in 1:ncol(zstat_sig)){
    if (i%%10==0) print(i)
    trait <- colnames(zstat_sig)[i]
    load(paste0("snps_clumped_",trait,"_pthr5e-08.rda"))
    
    if (nrow(snps.clumped)>=30){
        snps.clumped$beta <- beta[match(snps.clumped$SNP,rownames(beta)),trait]
        snps.clumped$beta_rank <- rank(abs(snps.clumped$beta))
        snps.clumped$beta_revrank <- rank(-abs(snps.clumped$beta))
        snps.clumped <- cbind(snps.clumped, numtraits.hardthr[match(snps.clumped$SNP,rownames(numtraits.hardthr)),])
        snps.clumped$trait <- trait
        beta_numtraits_clump <- rbind(beta_numtraits_clump, snps.clumped)
    }
}

load("../res_all.rda")
res_all$numtraits <- sapply(res_all$pheno_pos,length)+sapply(res_all$pheno_neg,length)
ind <- match(beta_numtraits_clump$SNP,res_all$SNP)
beta_numtraits_clump$numtraits <- res_all$numtraits[ind]
beta_numtraits_clump$pval_asset <- res_all$pval[ind]
save(beta_numtraits_clump, file = "beta_numtraits_clump.rda")


## Plots
rm(list=ls())
library(data.table)
library(readxl)
library(dplyr)
library(ggplot2)
library(gridExtra)

load("PATH_TO_FILE/beta_numtraits_clump.rda")
trait.list.all <- read_excel("../support_files/supptab1_final_data_list_w_DISPLAYNAME.xlsx", sheet=1)

beta_numtraits_clump <- beta_numtraits_clump %>% filter(beta_revrank<=100)
ind <- match(beta_numtraits_clump$trait,trait.list.all$trait)
beta_numtraits_clump$Domain <- trait.list.all$Domain[ind]
beta_numtraits_clump$trait_display <- trait.list.all$Display_name[ind]
beta_numtraits_clump$N <- trait.list.all$N[ind]
beta_numtraits_clump %>% filter(grepl("/",N)) %>% select(trait) %>% distinct()

beta_numtraits_clump_disease <- beta_numtraits_clump %>% filter(grepl("/",N))
temp <- beta_numtraits_clump_disease %>% select(trait_display,Domain) %>% distinct()
beta_numtraits_clump_disease$trait_display <- 
    factor(beta_numtraits_clump_disease$trait_display, temp$trait_display)

beta_numtraits_clump_biomarker <- beta_numtraits_clump %>% 
    filter(Domain=="Biomarker" & !grepl("/",N)) %>%
    filter(as.integer(N)>=1e5)
temp <- beta_numtraits_clump_biomarker %>% select(trait_display,Domain) %>% distinct()
beta_numtraits_clump_biomarker$trait_display <- 
    factor(beta_numtraits_clump_biomarker$trait_display, temp$trait_display)

beta_numtraits_clump_anthssl <- beta_numtraits_clump %>% 
    filter((Domain=="Social science and lifestyle"|Domain=="Anthropometric") & !grepl("/",N)) %>%
    filter(as.integer(N)>=1e5)
temp <- beta_numtraits_clump_anthssl %>% select(trait_display,Domain) %>% distinct()
beta_numtraits_clump_anthssl$trait_display <- 
    factor(beta_numtraits_clump_anthssl$trait_display, temp$trait_display)

gg1 <- ggplot(beta_numtraits_clump_disease,
              aes(x=beta_revrank,y=`pthr 1e-04`)) +
    geom_point(size=0.7) +
    geom_smooth(method="loess") +
    labs(x="Rank of effect size (largest to smallest)", y="Number of traits p<1e-4") +
    theme_bw() +
    facet_wrap(~trait_display, nrow=1, scale="free")
gg2 <- ggplot(beta_numtraits_clump_biomarker,
              aes(x=beta_revrank,y=`pthr 1e-04`)) +
    geom_point(size=0.7) +
    geom_smooth(method="loess") +
    labs(x="Rank of effect size (largest to smallest)", y="Number of traits p<1e-4") +
    theme_bw() +
    facet_wrap(~trait_display, nrow=1, scale="free")
gg3 <- ggplot(beta_numtraits_clump_anthssl,
              aes(x=beta_revrank,y=`pthr 1e-04`)) +
    geom_point(size=0.7) +
    geom_smooth(method="loess") +
    labs(x="Rank of effect size (largest to smallest)", y="Number of traits p<1e-4") +
    theme_bw() +
    facet_wrap(~trait_display, nrow=1, scale="free")

gg4 <- ggplot(beta_numtraits_clump %>% filter(beta_revrank<=30), 
       aes(x=reorder(trait_display,numtraits,FUN=median),y=numtraits,fill=Domain)) +
    geom_boxplot() +
    labs(x = "Trait", y = "Number of traits by ASSET (top 30 hits)") +
    theme(axis.text.x = element_text(angle = 90, size=10, hjust=1, vjust=0.5),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_text(vjust=0.75))

png(filename = "effsize_pleio_all.png",
    width = 3500, height = 3500, res = 300)
grid.arrange(gg1,gg2,gg3,gg4,
             ncol = 1, nrow = 5, 
             layout_matrix = matrix(c(1,2,3,5,5),ncol=1))
dev.off()
