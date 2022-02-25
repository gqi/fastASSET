# Plot top 10 associateion for disease-traits specific SNPs and matched pleiotropic SNPs
# This script was used to create Figure 4
# Binary disease traits
rm(list=ls())
library(readxl)
library(dplyr)
library(reshape2)
library(RColorBrewer)

load("PATH_TO_FILE/loci1t_match_specpleio.rda")
top_pleio <- loci1t %>% group_by(pheno_display) %>% summarize(top_pleio_match=SNP_pleio_match[1])
loci1t$SNP_pleio_match <- top_pleio$top_pleio_match[match(loci1t$pheno_display,top_pleio$pheno_display)]

spec <- c("rs6733839","rs7667257","rs12653946","rs2260801","rs4631223","rs4958425")
loci1t <- loci1t %>% filter(SNP%in%spec) %>% mutate(snppheno=paste(SNP,pheno_display))
loci1t.pleio <- loci1t %>% select(pheno_display, SNP_pleio_match) %>% distinct()

trait.list.all <- read_excel("../support_files/supptab1_final_data_list_w_DISPLAYNAME.xlsx", sheet=1)

load("PATH_TO_FILE/snps_clumped_sumstats.rda")
snps_clumped_sumstats <- snps_clumped_sumstats[match(c(loci1t$SNP,unique(loci1t$SNP_pleio_match)),
                                                     snps_clumped_sumstats$SNP),]
zstat <- as.matrix(snps_clumped_sumstats[,10:125])/
    as.matrix(snps_clumped_sumstats[,126:241])
pval <- 2*pnorm(-abs(zstat))
rownames(pval) <- c(loci1t$SNP,unique(loci1t$SNP_pleio_match))
colnames(pval) <- trait.list.all$Display_name[match(gsub("beta_","",colnames(pval)),trait.list.all$trait)]

##
pval.lnplt <- melt(pval) %>% filter(!is.na(value)) %>%
    mutate(value = -log10(value)) %>%
    group_by(Var1) %>% 
    mutate(valrank=rank(-value)) %>%
    ungroup() %>%
    filter(valrank<=10) %>%
    rename(SNP=Var1, trait=Var2)

snpvec <- levels(pval.lnplt$SNP)

png(filename="phewas_lineplot_matched_disease.png",
    width=500*4, height=500*3, res=300)
par(mfrow=c(3,4), mar=c(5,2.7,2,1))

for (i in 1:nrow(loci1t)){
    if (i==4) plot.new()
    
    snp <- loci1t$SNP[i]
    df <- pval.lnplt %>%
        filter(SNP==snp) %>% arrange(valrank) 
    plot(df$valrank, df$value, xaxt = "n", pch=19, ylim=c(0,max(df$value)),
         xlab="", ylab="", main=loci1t$snppheno[i], cex.main=0.7, cex.axis=0.7,
         bty="n", cex=0.8)
    lines(df$valrank, df$value)
    dist <- df$value[1]-df$value[2]
    segments(x0=2.5,y0=df$value[2]+0.01*dist, x1=2.5,y1=df$value[2]+0.4*dist)
    segments(x0=2.5,y0=df$value[1]-0.01*dist, x1=2.5,y1=df$value[1]-0.4*dist)
    
    segments(x0=2.3,y0=df$value[2]+0.01*dist, x1=2.7,y1=df$value[2]+0.01*dist)
    segments(x0=2.3,y0=df$value[1]-0.01*dist, x1=2.7,y1=df$value[1]-0.01*dist)
    text(x=2.5,y=mean(df$value[1:2]), labels=round(dist,1), cex=0.75, adj=0.5)
    
    axis(1, at=1:10, labels=df$trait, cex.axis=0.5, las=2)
    title(ylab="-log10(P)", line=2, cex.lab=0.7)
}

plot.new()

for (i in 1:nrow(loci1t.pleio)){
    snp <- loci1t.pleio$SNP_pleio_match[i]
    df <- pval.lnplt %>%
        filter(SNP==snp) %>% arrange(valrank)

    plot(df$valrank, df$value, xaxt = "n", pch=19, ylim=c(0,max(df$value)),
         xlab="", ylab="", main=paste0(snp,"\n(Matched to ",loci1t.pleio$pheno_display[i],")"), 
         cex.main=0.7, cex.axis=0.7, bty="n", cex=0.8)
    lines(df$valrank, df$value)
    axis(1, at=1:10, labels=df$trait, cex.axis=0.5, las=2)
    title(ylab="-log10(P)", line=2, cex.lab=0.7)
}
dev.off()

# Continuous traits
rm(list=ls())
library(readxl)
library(dplyr)
library(reshape2)
library(RColorBrewer)

# List of trait-specific SNPs
spec <- c("rs6733839","rs7667257","rs12653946","rs2260801","rs4631223","rs4958425")
load("ASSET_final_analysis/results/20_01_31_ASSET_GRASP_cornumr_F/clump/loci1t_match_specpleio.rda")
top_pleio <- loci1t %>% group_by(pheno_display) %>% summarize(top_pleio_match=SNP_pleio_match[1])
loci1t$SNP_pleio_match <- top_pleio$top_pleio_match[match(loci1t$pheno_display,top_pleio$pheno_display)]

loci1t <- loci1t %>% filter(!(SNP%in%spec)) %>% mutate(snppheno=paste(SNP,pheno_display))
loci1t.pleio <- loci1t %>% select(pheno_display, SNP_pleio_match) %>% distinct()
trait.list.all <- read_excel("../support_files/supptab1_final_data_list_w_DISPLAYNAME.xlsx", sheet=1)

load("PATH_TO_FILE/snps_clumped_sumstats.rda")
snps_clumped_sumstats <- snps_clumped_sumstats[match(c(loci1t$SNP,unique(loci1t$SNP_pleio_match)),
                                                     snps_clumped_sumstats$SNP),]
zstat <- as.matrix(snps_clumped_sumstats[,10:125])/
    as.matrix(snps_clumped_sumstats[,126:241])
pval <- 2*pnorm(-abs(zstat))
rownames(pval) <- c(loci1t$SNP,unique(loci1t$SNP_pleio_match))
colnames(pval) <- trait.list.all$Display_name[match(gsub("beta_","",colnames(pval)),trait.list.all$trait)]

##
pval.lnplt <- melt(pval) %>% filter(!is.na(value)) %>%
    mutate(value = -log10(value)) %>%
    group_by(Var1) %>% 
    mutate(valrank=rank(-value)) %>%
    ungroup() %>%
    filter(valrank<=10) %>%
    rename(SNP=Var1, trait=Var2)

snpvec <- levels(pval.lnplt$SNP)

png(filename="phewas_lineplot_matched_continuous.png",
    width=500*5, height=500*5, res=300)
par(mfrow=c(5,5), mar=c(5,2.7,2,1))

for (i in 1:nrow(loci1t)){
    snp <- loci1t$SNP[i]
    df <- pval.lnplt %>%
        filter(SNP==snp) %>% arrange(valrank) %>%
        mutate(trait=factor(trait,levels=trait))
    plot(df$valrank, df$value, xaxt = "n", pch=19, ylim=c(0,max(df$value)),
         xlab="", ylab="", main=loci1t$snppheno[i], cex.main=0.7, cex.axis=0.7,
         bty="n", cex=0.8)
    lines(df$valrank, df$value)
    dist <- df$value[1]-df$value[2]
    segments(x0=2.5,y0=df$value[2]+0.01*dist, x1=2.5,y1=df$value[2]+0.4*dist)
    segments(x0=2.5,y0=df$value[1]-0.01*dist, x1=2.5,y1=df$value[1]-0.4*dist)
    
    segments(x0=2.3,y0=df$value[2]+0.01*dist, x1=2.7,y1=df$value[2]+0.01*dist)
    segments(x0=2.3,y0=df$value[1]-0.01*dist, x1=2.7,y1=df$value[1]-0.01*dist)
    text(x=2.5,y=mean(df$value[1:2]), labels=round(dist,1), cex=0.75, adj=0.5)
    
    axis(1, at=1:10, labels=df$trait, cex.axis=0.5, las=2)
    title(ylab="-log10(P)", line=2, cex.lab=0.7)
}

for (i in 1:nrow(loci1t.pleio)){
    snp <- loci1t.pleio$SNP_pleio_match[i]
    df <- pval.lnplt %>%
        filter(SNP==snp) %>% arrange(valrank) 
    
    plot(df$valrank, df$value, xaxt = "n", pch=19, ylim=c(0,max(df$value)),
         xlab="", ylab="", main=paste0(snp,"\n(Matched to ",loci1t.pleio$pheno_display[i],")"), 
         cex.main=0.7, cex.axis=0.7, bty="n", cex=0.8)
    lines(df$valrank, df$value)
    axis(1, at=1:10, labels=df$trait, cex.axis=0.5, las=2)
    title(ylab="-log10(P)", line=2, cex.lab=0.7)
}

dev.off()
