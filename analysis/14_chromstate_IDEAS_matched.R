# Chromatin states results - Figure 6
# Run in clump
rm(list=ls())
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggthemes)
library(gridExtra)

ideas.bed <- fread("PATH_TO_FILE/ideas_loci_snps_all.bed")
colnames(ideas.bed) <- c("CHR","BP","BPplus1","SNP","pval","pheno","numtraits","CHR_ideas","BP_start_ideas","BP_end_ideas","state", "eid")
chromstate <- ideas.bed %>% filter(state%in%c("10_TssA","8_TssAFlnk","2_TxWk","14_TssWk","4_Enh","6_EnhG","17_EnhGA","5_Tx"))

# Add group labels
roadmap.annot <- read_excel("PATH_TO_FILE/Roadmap.metadata.qc.jul2013.xlsx",
                            sheet = 1)[-c(1,2),] %>% 
    select(eid = `Epigenome ID (EID)`, GROUP, COLOR)
roadmap.color <- roadmap.annot %>% select(GROUP,COLOR) %>% distinct()
roadmap.colorvec <- roadmap.color$COLOR
names(roadmap.colorvec) <- roadmap.color$GROUP
chromstate$group <- roadmap.annot$GROUP[match(chromstate$eid,roadmap.annot$eid)]

chromstate <- chromstate %>% filter(group!="ENCODE2012")

# Load lists of trait-specific and pleiotropic SNPs
trait.list.all <- read_excel("../support_files/supptab1_final_data_list_w_DISPLAYNAME.xlsx", sheet=1)

pheno.order <- c("alzheimer","BrCa","PrCa","CD","BMDheel","DBP","Height18",
                 "intelligence","male_baldness","menarche","mono")
load("loci1t_match_specpleio.rda")
top_pleio <- loci1t %>% group_by(pheno_display) %>% summarize(top_pleio_match=SNP_pleio_match[1])
loci1t$SNP_pleio_match <- top_pleio$top_pleio_match[match(loci1t$pheno_display,top_pleio$pheno_display)]

loci1t <- loci1t %>%
    filter(SNP%in%chromstate$SNP) %>%
    mutate(pheno = factor(pheno, levels=pheno.order)) %>%
    arrange(pheno) %>%
    mutate(snppheno = paste0(SNP,"\n",pheno_display)) %>% 
    mutate(snppheno = factor(snppheno,levels=snppheno))

chromstate.1t <- chromstate %>% filter(SNP%in%loci1t$SNP)
chromstate.1t$snppheno <- loci1t$snppheno[match(chromstate.1t$SNP,loci1t$SNP)]
chromstate.mt15t <- chromstate %>% filter(SNP%in%loci1t$SNP_pleio_match)
chromstate.mt15t$snppheno <- factor(chromstate.mt15t$SNP,levels=unique(loci1t$SNP_pleio_match))
levels(chromstate.mt15t$snppheno) <- paste0(levels(chromstate.mt15t$snppheno),
                                            "\n(",top_pleio$pheno_display[match(levels(chromstate.mt15t$snppheno),top_pleio$top_pleio_match)],")")

gg1 <- ggplot(chromstate.1t, aes(x=snppheno,fill=group)) +
    geom_bar() +
    theme(axis.text.x = element_text(size=7, angle = 90, vjust = 0.5, hjust=1),
          panel.background = element_blank(), 
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    ylim(0,120) +
    scale_fill_manual(values=roadmap.colorvec, name="Cell type") +
    labs(x = "SNP and trait", y = "Number of cell types in Roadmap", title = "Trait specific")

gg2 <- ggplot(chromstate.mt15t, aes(x=snppheno,fill=group)) +
    geom_bar() +
    theme(axis.text.x = element_text(size=7, angle = 90, vjust = 0.5, hjust=1),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    ylim(0,120) +
    scale_fill_manual(values=roadmap.colorvec[names(roadmap.colorvec)!="ENCODE2012"], name="Cell type") +
    labs(x = "SNP", y = "Number of cell types in Roadmap", title = "Highly pleiotropic")

png("chromstate_1ormt15t_matched_IDEAS.png",
    height = 2000, width = 3500, res = 300)
grid.arrange(gg1,gg2,nrow=1)
dev.off()



