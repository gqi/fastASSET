## Check and LD clump SNPs that were not analyzed by ASSET (too many traits passed pre-screening)
# Number of traits that passed ASSET pre-screening threshold
rm(list=ls())
numtraits_thr = vector(length = 0) # number of traits that pass the pre-screening treshold
for (i in 1:18){
    nbatch = ifelse(i<18,10,2)
    for (j in 1:nbatch){
        load(paste0("../res_split",i,"_",j,".rda"))
        numtraits_thr = c(numtraits_thr, sapply(res$beta_scr,length))
        rm(res)
    }
    print(c(i,length(numtraits_thr)))
}
save(numtraits_thr, file = "numtraits_thr.rda")

# Check pleiotropic SNPs: run in folder check_pleio_snps
rm(list=ls())
library(data.table)
library(dplyr)

bim = fread("PATH_TO_1000G_DATA/1000G.EUR.QC.1.bim")
for (i in 2:22){
    print(i)
    bim = rbind(bim, fread(paste0("PATH_TO_1000G_DATA/1000G.EUR.QC.",i,".bim")))
}
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")

load("../res_all.rda")
sig = res_all %>% filter(pval<5e-8) %>%
    left_join(bim,by="SNP") %>% select(CHR, BP, SNP, pval, error_warning)
pleio = res_all %>% filter(error_warning==6) %>%
    left_join(bim,by="SNP") %>% select(CHR, BP, SNP, pval, error_warning)
rm(res_all,bim)

for (i in 1:22){
    print(i)
    snplist = rbind(sig,pleio) %>% filter(CHR==i) %>% select(SNP)
    write.table(snplist, file = paste0("snplist_chr",i,".txt"), quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
}

plinkpath = "PATH_TO_PLINK_SOFTWARE"
KGpath = "PATH_TO_1000G_DATA"
outpath = "."
for (chrnum in 1:22){
    plinkcode = paste(paste0(plinkpath,"/plink"),
                      "--bfile", paste0(KGpath,"/1000G.EUR.QC.",chrnum),
                      "--extract", paste0(outpath,"/snplist_chr",chrnum,".txt"),
                      "--ld-window 99999",
                      "--ld-window-r2 0.1",
                      "--ld-window-kb 1000", 
                      "--r2",
                      "--out", paste0(outpath,"/chr_LD",chrnum))
    system(plinkcode)
}

pleio$R2max = NA
for (i in 1:22){
    print(i)
    pleio_temp = pleio %>% filter(CHR==i)
    ld = fread(paste0("chr_LD",i,".ld")) %>% filter(((SNP_A%in%pleio$SNP & SNP_B%in%sig$SNP)|(SNP_B%in%pleio$SNP & SNP_A%in%sig$SNP)))
    for (j in 1:nrow(pleio_temp)){
        snp = pleio_temp$SNP[j]
        ld_temp = ld %>% filter(SNP_A==snp|SNP_B==snp)
        pleio$R2max[pleio$SNP==snp] = max(ld_temp$R2)
        if (j%%100==0) print(j)
    }
}
pleio$dist = NA
for (i in 1:nrow(pleio)){
    pleio$dist[i] = min(abs(sig$BP[sig$CHR==pleio$CHR[i]]-pleio$BP[i]))
}
save(pleio, file = "pleio.rda")

## Select SNPs that are not analyzed by ASSET and not tagged by any other SNP with p<5e-8 
# Clump them to obtain independent loci
rm(list=ls())
library(dplyr)
load("pleio.rda")
load("numtraits_thr.rda")
pleio = pleio %>% mutate(R2max = ifelse(is.infinite(R2max),0.1,R2max)) %>%
    filter(R2max<0.5)
pleio$numtraits_thr = numtraits_thr[match(pleio$SNP,names(numtraits_thr))]

tempdata = pleio %>% select(SNP, CHR, BP)
tempdata$NMISS = NA
tempdata$BETA = NA 
tempdata$SE = NA 
tempdata$P = 1/(pleio$numtraits_thr+runif(nrow(pleio),0,0.1))
tempdata$R2 = NA
tempdata$T = NA

tempdata = tempdata %>% select(CHR = chr, SNP = rsid, BP = bp, NMISS, BETA, SE, R2, T, P) 

outpath <- "."
plinkpath = "PATH_TO_PLINK_SOFTWARE"
KGpath = "PATH_TO_1000G_DATA"
pthr = 0.1 
r2thr = 0.1
kbpthr = 500
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
pleio_untagged = pleio
pleio_untagged$lead_snp = FALSE
pleio_untagged$lead_snp[match(snps_clumped_all$SNP,pleio_untagged$SNP)] = TRUE

pleio_untagged = pleio_untagged %>% 
    filter(lead_snp) %>%
    select(CHR, BP, SNP, pval, error_warning)
pleio_untagged$numtraits_thr <- numtraits_thr[match(pleio_untagged$SNP,names(numtraits_thr))]
save(pleio_untagged, file = "pleio_untagged.rda")

