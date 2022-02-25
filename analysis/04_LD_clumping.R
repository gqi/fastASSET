## LD clumping to obtain SNPs of significant loci
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
        # temp$SP2 = gsub("\\(1\\)","",temp$SP2)
        if (nrow(temp)>1 & kbpthr>0){
            temp = temp[order(temp$P,decreasing=FALSE),]
            temp_pruned = NULL
            while (nrow(temp)>0){
                temp_pruned = rbind(temp_pruned,temp[1,])
                ind = abs(temp$BP-temp$BP[1])<=kbpthr*1000
                # Remove SNPs within kbpthr*1000 distance from the index SNP
                temp = temp[!ind,]
            }
            snps_clumped_all = rbind(snps_clumped_all, temp_pruned)
        } else{
            temp$SP2 = gsub("NONE","",temp$SP2) # Remove NONE's
            snps_clumped_all = rbind(snps_clumped_all, temp)
        }
    }
}

snps_clumped_all = snps_clumped_all %>% left_join(res_all,by="SNP") %>%
    select(CHR, BP, SNP, pval, error_warning, pheno_pos, pheno_neg)
snps_clumped_all$pheno = sapply(1:nrow(snps_clumped_all), 
                                function(x) c(snps_clumped_all$pheno_pos[[x]],snps_clumped_all$pheno_neg[[x]]))
snps_clumped_all$numtraits = sapply(snps_clumped_all$pheno,length)

save(snps_clumped_all, file = "snps_clumped_all.rda")
