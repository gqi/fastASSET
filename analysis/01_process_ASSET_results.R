## Process and filter results from fast ASSET analysis

## Assemble results
rm(list=ls())
library(dplyr)
library(data.table)
library(ASSET)

# Trait list
traitdf = fread("../support_files/datalist_lowh2higcorr_rm.csv")
# Load LDSC intercept matrix: correlation of z stats
load("../support_files/ldscintmat.rda")
ldscintmat = ldscintmat[traitdf$trait,traitdf$trait]
traits = colnames(ldscintmat)
ldscintmat = solve(diag(sqrt(diag(ldscintmat)))) %*% ldscintmat %*% solve(diag(sqrt(diag(ldscintmat))))
colnames(ldscintmat) = rownames(ldscintmat) = traits

pval = vector(length = 0)
error_warning = vector(length = 0)
pheno_pos = list()
pheno_neg = list()
ast1s_info = matrix(NA, nrow = 0, ncol = 2)
colnames(ast1s_info) = c("pval.1","pval.2")
z.meta = matrix(NA, nrow = 0, ncol = 2)
colnames(z.meta) = c("z1.meta","z2.meta")

calc.meta.zstat = function(x){
    ntraits = length(res$pheno_pos[[x]])
    if (ntraits>0){
        z1 = traits.meta(rep(TRUE,length(res$pheno_pos[[x]])), 1, res$beta_scr[[x]][res$pheno_pos[[x]]], 
                         res$sigma_scr[[x]][res$pheno_pos[[x]]], 
                         ncase=2/res$sigma_scr[[x]][res$pheno_pos[[x]]]^2, ncntl=2/res$sigma_scr[[x]][res$pheno_pos[[x]]]^2, 
                         rmat=ldscintmat[res$pheno_pos[[x]],res$pheno_pos[[x]]], side=2, cor.numr=F, wt.sigma=F)$z
    } else{
        z1 = 0
    }
    
    ntraits = length(res$pheno_neg[[x]])
    if (ntraits>0){
        z2 = traits.meta(rep(TRUE,length(res$pheno_neg[[x]])), 1, res$beta_scr[[x]][res$pheno_neg[[x]]], 
                         res$sigma_scr[[x]][res$pheno_neg[[x]]], 
                         ncase=2/res$sigma_scr[[x]][res$pheno_neg[[x]]]^2, ncntl=2/res$sigma_scr[[x]][res$pheno_neg[[x]]]^2, 
                         rmat=ldscintmat[res$pheno_neg[[x]],res$pheno_neg[[x]]], side=2, cor.numr=F, wt.sigma=F)$z
    } else{
        z2 = 0
    }
    
    c(z1,z2)
}

for (i in 1:18){
    nbatch = ifelse(i<18,10,2)
    for (j in 1:nbatch){
        load(paste0("res_split",i,"_",j,".rda"))
        pval = c(pval, res$pval_ast)
        error_warning = c(error_warning, res$error_warning)
        pheno_pos = c(pheno_pos, res$pheno_pos)
        pheno_neg = c(pheno_neg, res$pheno_neg)
        ast1s_info = rbind(ast1s_info, res$ast1s_info[,c("pval.1","pval.2")])
        
        z.meta = rbind(z.meta, t(sapply(1:length(res$pval_ast), calc.meta.zstat)))
        rm(res)
    }
    print(c(i,length(pval)))
}
res_all = tibble(SNP = names(pval), pval = pval, error_warning = error_warning,
                 pheno_pos = pheno_pos, pheno_neg = pheno_neg)
res_all = cbind(res_all, ast1s_info, z.meta)

# SNP filtering: 
# 1) Remove SNPs for which ASSET did not converge
# 2) Remove SNPs for which sumstats for <50 traits are available, indicated by (pval.1<pval.1.meta | pval.2<pval.2.meta)

# Number of available traits for each SNP
load("PATH_TO_DATA/numtraits_avail.rda")
res_all$numtraits_avail <- numtraits_avail[match(res_all$SNP,names(numtraits_avail))]
rm(numtraits_avail)
sum(res_all$numtraits_avail<50)
res_all = res_all %>% filter(numtraits_avail>=50) %>%
    mutate(pval.1.meta = 2*pnorm(z1.meta,lower.tail=FALSE), pval.2.meta = 2*pnorm(z2.meta))
snplist_all = res_all$SNP # Save full list of SNPs
res_ast_error = res_all %>% filter(error_warning==0 & (pval.1<pval.1.meta | pval.2<pval.2.meta))
dim(res_ast_error)
res_all = res_all %>% filter(!(SNP%in%res_ast_error$SNP))

save(res_all, file = "res_all.rda")
save(snplist_all, file = "snplist_all.rda")
