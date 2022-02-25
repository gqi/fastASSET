#' assoc is a data.frame with columns CHR, BP, SNP, P. Other columns are optional.
#' plink.path should end with /
clump = function(assoc_min, pthr = 5e-8, r2thr = 0.1, bpthr = 1e6, plink.path = "", 
                 ref.geno = paste0("FILE_NAME",1:22),
                 ld.window = 99999, ld.window.kb = 10000, threads = 1,
                 ld.out = "out"){
    # SNP list
    assoc_min = assoc_min %>% filter(P<pthr)
    for (i in 1:22){
        print(i)
        snplist.name = paste0("snplist_",ld.out,i,".txt")
        assoc_min %>% filter(CHR==i) %>% select(SNP) %>%
            write.table(file = snplist.name, row.names = F, col.names = F, quote = F)
        plink.code = paste(paste0(plink.path,"plink"), 
                           "--bfile", ref.geno[i],
                           "--extract", snplist.name,
                           "--r2",
                           "--ld-window", ld.window,
                           "--ld-window-kb", ld.window.kb,
                           "--ld-window-r2", r2thr,
                           "--threads", threads,
                           "--out", paste0(ld.out,i))
        system(plink.code)
    }
    
    assoc_clumped = data.frame()
    for (chrnum in 1:22){
        print(chrnum)
        assoc = assoc_min %>% filter(CHR==chrnum) %>% arrange(P)
        if (nrow(assoc)<=1){
            assoc_clumped = assoc_clumped = rbind(assoc_clumped, assoc)
        } else{
            ld = fread(paste0(ld.out,chrnum,".ld"))
            while(nrow(assoc)>0){
                assoc_clumped = rbind(assoc_clumped, assoc[1,])
                assoc = assoc %>% filter(abs(BP-BP[1])>bpthr & 
                                             !(SNP %in% c(ld$SNP_A[ld$SNP_B==assoc_clumped$ID[nrow(assoc_clumped)]],
                                                          ld$SNP_B[ld$SNP_A==assoc_clumped$ID[nrow(assoc_clumped)]]))) 
            }
            
        }
    }
    system(paste0("rm ",ld.out,"*"))
    system(paste0("rm snplist_",ld.out,"*"))
    return(assoc_clumped)
}


