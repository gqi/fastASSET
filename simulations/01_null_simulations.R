## Null simulations to evaluate type I error of fastASSET, metaUSAT and metaMANOVA
rm(list=ls())
library(MASS)
library(data.table)
library(ASSET)
source("../support_files/ASSET_block_orth.R")
source("PATH_TO_METAUSAT_SCRIPT/metaUSAT_v1.17.R")

# Trait list
traitdf = fread("../support_files/datalist_lowh2higcorr_rm.csv")
# Load LDSC intercept matrix: correlation of z stats
load("../support_files/ldscintmat.rda")
ldscintmat = ldscintmat[traitdf$trait,traitdf$trait]
traits = colnames(ldscintmat)
ldscintmat = solve(diag(sqrt(diag(ldscintmat)))) %*% ldscintmat %*% solve(diag(sqrt(diag(ldscintmat))))
colnames(ldscintmat) = rownames(ldscintmat) = traits

# Hierarchical clustering
corrdist = as.dist(1-abs(ldscintmat))
hc = hclust(corrdist)
htree = cutree(hc, h=0.8) # Criterion: corr>0.2
block = as.integer(names(table(htree))[table(htree)>=2])
block = lapply(block, function(x) names(htree)[htree==x])

n = 20000
nS = 200000
level=5e-8

temp = as.integer(commandArgs(trailingOnly = TRUE)) # 100 jobs
set.seed(675453*temp)
betahat = mvrnorm(nS, mu = rep(0,nrow(ldscintmat)), Sigma = ldscintmat)/sqrt(n) # error term
Neff = rep(n,nrow(ldscintmat))
names(Neff) = colnames(betahat)

# Create list/data.frame to store results
res_alltraits <- matrix(NA,nrow = nS, ncol = 5, dimnames=list(NULL,c("pval","pval.1","pval.2","z1.meta","z2.meta")))
res <- data.frame(pval.metausat=rep(NA, nS), errmsg.metausat=rep(NA, nS), 
                  pval.metamanova=rep(NA, nS))
for (i in 1:nS){
    print(i)
    # Apply fast ASSET. This uses an older versions of the codes (ASSET_block_orth.R).
    # It can now be run by calling wrapper function fast_asset().
    dt_scr = scr.transform.block(betahat = betahat[i,], SE = 1/sqrt(Neff), 
                                 SNP = paste0("SNP",i), colnames(ldscintmat), cor = ldscintmat, block, scr_pthr = 0.05)
    betahat_orth = dt_scr$betahat_orth
    sigma_orth = dt_scr$sigma_orth
    traits_i=names(betahat_orth)
    
    if (length(traits_i)>=2){
        res.asset = h.traits(snp.vars=paste0("SNP",i), traits.lab=traits_i,
                             beta.hat=betahat_orth, sigma.hat=sigma_orth, ncase=2/sigma_orth^2,                                    
                             ncntl=2/sigma_orth^2, cor=ldscintmat[traits_i,traits_i], 
                             side=2, search=2, cor.numr = F)
        res_alltraits[i,"pval"] = res.asset$Subset.2sided$pval
        res_alltraits[i,"pval.1"] = res.asset$Subset.2sided$pval.1
        res_alltraits[i,"pval.2"] = res.asset$Subset.2sided$pval.2
        
        tt = colnames(res.asset$Subset.2sided$pheno.1)[res.asset$Subset.2sided$pheno.1]
        if (length(tt)>0){
            res_alltraits[i,"z1.meta"] = ASSET:::traits.meta(rep(TRUE,length(tt)), 1, betahat_orth[tt], sigma_orth[tt], ncase=2/sigma_orth[tt]^2, ncntl=2/sigma_orth[tt]^2, rmat=ldscintmat[tt,tt], side=2, cor.numr=F, wt.sigma=F)$z
        } else{
            res_alltraits[i,"z1.meta"] = 0
        }
        
        tt = colnames(res.asset$Subset.2sided$pheno.2)[res.asset$Subset.2sided$pheno.2]
        if (length(tt)>0){
            res_alltraits[i,"z2.meta"] = ASSET:::traits.meta(rep(TRUE,length(tt)), 1, betahat_orth[tt], sigma_orth[tt], ncase=2/sigma_orth[tt]^2, ncntl=2/sigma_orth[tt]^2, rmat=ldscintmat[tt,tt], side=2, cor.numr=F, wt.sigma=F)$z
        } else{
            res_alltraits[i,"z2.meta"] = 0
        }
        
    } else if (length(traits_i)==1){
        res_alltraits[i,"pval"] = as.numeric(2*pnorm(-abs(betahat_orth/sigma_orth)))
    } else{
        res_alltraits[i,"pval"] = 1 
    }
    
    # Apply metaUSAT
    res.metausat = metausat(betahat[i,]*sqrt(n), ldscintmat)
    res$pval.metausat[i] = res.metausat$p.metausat
    res$errmsg.metausat[i] = res.metausat$error.msg
    
    # Apply metaMANOVA
    res$pval.metamanova[i] = pchisq(n*diag(t(betahat[i,])%*%solve(ldscintmat,betahat[i,])), 
                                df=nrow(ldscintmat), lower.tail = FALSE)
    # Save results
    if (i%%500==0) {
        print(i)
        save(res_alltraits,res, file = paste0("pval",temp,".rda"))
    }
}

## Filter and summarize results
rm(list=ls())
library(dplyr)
res_alltraits_all = matrix(NA, nrow = 0, ncol = 5)
for (i in 1:100){
    print(i)
    load(paste0("pval",i,".rda"))
    res_alltraits_all = rbind(res_alltraits_all,res_alltraits)
}
dd = data.frame(res_alltraits_all)
dd = dd %>% filter(!is.na(pval))
dd = dd %>% mutate(pval.1.meta = 2*pnorm(z1.meta,lower.tail=FALSE), 
                   pval.2.meta = 2*pnorm(z2.meta),
                   pval.1 = ifelse(is.na(pval.1),1,pval.1), 
                   pval.2 = ifelse(is.na(pval.2),1,pval.2))

# If the either of the one-sided raw meta-analysis p-value is smaller than the fast ASSET p-value,
# it indicates convergence issue of the fast ASSET p-value approximation. Remove those SNPs.
ddflt = dd %>% mutate(pval.1.meta = ifelse(is.na(pval.1.meta),1,pval.1.meta),
                      pval.2.meta = ifelse(is.na(pval.2.meta),1,pval.2.meta)) %>% 
    filter(pval.1>=pval.1.meta & pval.2>=pval.2.meta)

sapply(c(1e-4,1e-5,1e-6,1e-7), function(x) mean(ddflt$pval<x))
