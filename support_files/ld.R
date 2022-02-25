#' assoc is a data.frame with columns CHR, BP, SNP, P. Other columns are optional.
#' plink.path should end with /
ld = function(snplist, r2thr = 0.1, plink.path = "/home-2/gqi1@jhu.edu/data/gqi1/plink/", 
                 ref.geno = paste0("/home-2/gqi1@jhu.edu/data/gqi1/1000G_EUR_Phase3_plink/1000G.EUR.QC.",1:22),
                 chr, ld.window = 99999, ld.window.kb = 10000, threads = 1,
                 ld.out = "LDinfo"){
    # SNP list
    write.table(snplist, file = paste0(ld.out,"_snplist.txt"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    plink.code = paste(paste0(plink.path,"plink"), 
                       "--bfile", ref.geno[chr],
                       "--extract", paste0(ld.out,"_snplist.txt"),
                       "--r2",
                       "--ld-window", ld.window,
                       "--ld-window-kb", ld.window.kb,
                       "--ld-window-r2", r2thr,
                       "--threads", threads,
                       "--out", ld.out)
    system(plink.code)
    
    ldres = fread(paste0(ld.out,".ld"), nThread = 1)
    
    system(paste0("rm ",paste0(ld.out,"*")))
    return(ldres)
}


