## Brisbane plot - Figure 3
# This type of figure is originally created by Yengo, Loic, et al. "A Saturated Map of Common Genetic Variants Associated with Human Height from 5.4 Million Individuals of Diverse Ancestries." bioRxiv (2022).

# Step 1: repeat LD clumping
rm(list=ls())
library(data.table)
library(dplyr)
bim = fread("PATH_TO_1000G_DATA/1000G.EUR.QC.1.bim")
for (i in 2:22){
    print(i)
    bim = rbind(bim, fread(paste0("PATH_TO_1000G_DATA/1000G.EUR.QC.",i,".bim")))
}
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")

outpath <- "."
KGpath <- "PATH_TO_1000G_DATA"
plinkpath <- "PATH_TO_PLINK_SOFTWARE"
pthr = 5e-8 
r2thr = 0.1
kbpthr = 500

load("../res_all.rda") # & pval>0
# This step automatically removes SNPs with error_warning==6
sumstats = res_all %>% filter(pval<pthr) %>% left_join(bim,by="SNP") %>%
    select(SNP, CHR, BP, pval)

tempdata = sumstats %>% select(SNP, CHR, BP)
tempdata$NMISS = NA
tempdata$BETA = NA 
tempdata$SE = NA 

# Find SNPs that are significant for any of the traits, then partition by individual traits, genpcs
tempdata$P = sumstats$pval

tempdata$R2 = NA
tempdata$T = NA

tempdata = tempdata %>% select(CHR = chr, SNP = rsid, BP = bp, NMISS, BETA, SE, R2, T, P) 

# Call PLINK from R to conduct preliminary LD clumping
for (chrnum in 1:22){
    if (sum(tempdata$CHR==chrnum)>0){
        tempdata %>% filter(CHR==chrnum) %>% write.table(paste0(outpath,"/chr",chrnum,".qassoc"), row.names = FALSE, quote = FALSE, sep = "\t")
        plinkcode = paste(paste0(plinkpath,"/plink"),
                          "--bfile", paste0(KGpath,"/1000G.EUR.QC.",chrnum),
                          "--clump", paste0(outpath,"/chr",chrnum,".qassoc"),
                          "--clump-p1", pthr,
                          "--clump-p2", pthr,
                          "--clump-r2", r2thr, # 
                          "--clump-kb", 1000,
                          "--out", paste0(outpath,"/chr",chrnum,"_clump"))
        system(plinkcode)
    }
}

snps_clumped_all = NULL
# Read the .clumped file into R
for (chrnum in 1:22){
    filename = paste0(outpath,"/chr",chrnum,"_clump.clumped")
    if (file.exists(filename)){
        temp = read.table(filename, header = TRUE, stringsAsFactors = FALSE)
        snps_clumped_all = rbind(snps_clumped_all, temp)
    }
}

snps_clumped_all = snps_clumped_all %>% left_join(res_all,by="SNP") %>%
    select(CHR, BP, SNP, pval, error_warning, pheno_pos, pheno_neg)
snps_clumped_all <- snps_clumped_all %>% 
    mutate(numtraits=sapply(pheno_pos,length)+sapply(pheno_neg,length))
save(snps_clumped_all, file="snps_clumped_no_bpthr.rda")

# Run on local computer in folder brisbane_plot
rm(list=ls())
library(data.table)
library(tidyverse)
library(qqman)
library(gridExtra)

# snps_clumped_combn.rda includes top SNPs for all significant loci (p<5e-8) after LD clumping at r2<0.1 and >500kb apart
load("../support_files/snps_clumped_combn.rda")
# snps_clumped_no_bpthr.rda includes independent (r2<0.01) SNPs with p-value<5e-8, but the SNPs can be close to each other in physical distance.
load("../support_files/snps_clumped_no_bpthr.rda")

snps_clumped_combn$numhits100k <- snps_clumped_combn$numtraits_avg <- 
    snps_clumped_combn$numsnps1trait <- 
    snps_clumped_combn$numsnps1_3traits <- snps_clumped_combn$numsnps4_10traits <-
    numsnps_mt11traits <- NA
for (i in 1:nrow(snps_clumped_combn)){
    if (i%%500==0) print(i)
    temp <- snps_clumped_all %>% filter(CHR==snps_clumped_combn$CHR[i] &
                                            abs(BP-snps_clumped_combn$BP[i])<1e5) %>%
        select(CHR,BP,SNP,numtraits)
    
    snps_clumped_combn$numhits100k[i] <- nrow(temp)
    snps_clumped_combn$numtraits_avg[i] <- mean(temp$numtraits[temp$SNP!=snps_clumped_combn$SNP[i]])
    snps_clumped_combn$numsnps1trait[i] <- mean(temp$numtraits==1)
    snps_clumped_combn$numsnps2_5traits[i] <- mean(temp$numtraits>=2 & temp$numtraits<=5)
    snps_clumped_combn$numsnps6_10traits[i] <- mean(temp$numtraits>=6 & temp$numtraits<=10)
    snps_clumped_combn$numsnps_mt11traits[i] <- mean(temp$numtraits>=11)
}

ggplot(snps_clumped_combn, aes(x=numtraits, y=numtraits_avg)) +
    geom_point() +
    geom_smooth()

bim = fread("PATH_TO_1000G_DATA/1000G.EUR.QC.1.bim")
for (i in 2:22){
    print(i)
    bim = rbind(bim, fread(paste0("PATH_TO_1000G_DATA/1000G.EUR.QC.",i,".bim")))
}
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")

chr_info = bim %>% group_by(CHR) %>% summarise(chr_len=as.numeric(max(BP))) %>%
    mutate(chrstart_manhattan = cumsum(c(0,chr_len[1:21])) + mean(chr_len)*0.5*(0:21), # Place a gap between each pair of chromosomes
           mid = chrstart_manhattan+chr_len/2)

snps_clumped_combn = snps_clumped_combn %>% mutate(chrpos = chr_info$chrstart_manhattan[CHR]+BP,
                             numtraits_color = cut(numtraits, breaks=c(0.5,1.5,5.5,10.5,15.5,Inf), 
                                         labels=c("1","2-5","6-10","11-15",">15")),
                             numtraits_avg_color = cut(numtraits_avg, breaks=c(0.5,1.5,5.5,10.5,15.5,Inf), 
                                                   labels=c("1","2-5","6-10","11-15",">15")))

# Make plot
gg1 <- ggplot(snps_clumped_combn, aes(x=chrpos, y=numhits100k)) + # , size = 0.2
    geom_point(aes(color=numtraits_color), shape=15) +
    scale_color_manual(values=c("1"="red","2-5"="grey50","6-10"="lightskyblue",
                                "11-15"="blue",">15"="darkgoldenrod1"),
                       name="Number of traits\nassociated with\nlead SNP",
                       guide=guide_legend(override.aes = list(size = 3))) +
    # custom X axis:
    scale_x_continuous(label = chr_info$CHR, breaks= chr_info$mid , position = "bottom") +
    # Custom the theme:
    theme_bw() +
    theme(# text = element_text(size=18),
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) + 
    labs(x = "CHR", y = "Number of independent\nsignals within 100kb", title="(a)")

gg2 <- ggplot(snps_clumped_combn, aes(x=chrpos, y=numhits100k)) + # , size = 0.2
    geom_point(aes(color=100*numsnps1trait), shape=15) +
    scale_color_gradient(low="grey80",high="red3") +
    # custom X axis:
    scale_x_continuous(label = chr_info$CHR, breaks= chr_info$mid , position = "bottom") +
    # Custom the theme:
    theme_bw() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) + 
    labs(x = "CHR", y = "Number of independent\nsignals within 100kb", 
         title="(b)", color="% 1 trait")

gg3 <- ggplot(snps_clumped_combn, aes(x=chrpos, y=numhits100k)) + # , size = 0.2
    geom_point(aes(color=100*numsnps2_5traits), shape=15) +
    scale_color_gradient(low="grey80",high="saddlebrown") +
    # custom X axis:
    scale_x_continuous(label = chr_info$CHR, breaks= chr_info$mid , position = "bottom") +
    # Custom the theme:
    theme_bw() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) + 
    labs(x = "CHR", y = "Number of independent\nsignals within 100kb", 
         title="(c)", color="% 2-5 traits")

gg4 <- ggplot(snps_clumped_combn, aes(x=chrpos, y=numhits100k)) + # , size = 0.2
    geom_point(aes(color=100*numsnps6_10traits), shape=15) +
    scale_color_gradient(low="grey80",high="blue") +
    # custom X axis:
    scale_x_continuous(label = chr_info$CHR, breaks= chr_info$mid , position = "bottom") +
    # Custom the theme:
    theme_bw() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) + 
    labs(x = "CHR", y = "Number of independent\nsignals within 100kb", 
         title="(d)", color="% 6-10 traits")

gg5 <- ggplot(snps_clumped_combn, aes(x=chrpos, y=numhits100k)) + # , size = 0.2
    geom_point(aes(color=100*numsnps_mt11traits), shape=15) +
    scale_color_gradient(low="grey80",high="navy") +
    # custom X axis:
    scale_x_continuous(label = chr_info$CHR, breaks= chr_info$mid , position = "bottom") +
    # Custom the theme:
    theme_bw() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) + 
    labs(x = "CHR", y = "Number of independent\nsignals within 100kb", 
         title="(e)", color="% >10 traits")

png(filename = "brisbane.png", width = 3200, height = 830*5, res = 300, type = "cairo")
grid.arrange(gg1,gg2,gg3,gg4,gg5,ncol=1)
dev.off()

write.csv(snps_clumped_combn %>% select(SNP,CHR,BP,numhits100k,numsnps1trait,numsnps2_5traits,numsnps6_10traits,numsnps_mt11traits), 
          file="snps_clumped_combn_w_brisbane.csv", row.names = FALSE)
