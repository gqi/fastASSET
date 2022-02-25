## Non-null simulations to evaluate power
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
nS = 2000
level=5e-8

effvec <- c(0.02,0.03)

temp0 = as.integer(commandArgs(trailingOnly = TRUE)) # 1:2;  1:50
eff <- effvec[temp0[1]]
temp <- temp0[2]
set.seed(675453*temp)

betahat = mvrnorm(nS, mu = rep(0,nrow(ldscintmat)), Sigma = ldscintmat)/sqrt(n) # error term
ind <- sample(100, 10)
betahat[,ind] = betahat[,ind] + eff

Neff = rep(n,nrow(ldscintmat))
names(Neff) = colnames(betahat)

res_alltraits = matrix(NA,nrow = nS, ncol = 5)
colnames(res_alltraits) = c("pval","pval.1","pval.2","z1.meta","z2.meta")
res <- data.frame(pval.metausat=rep(NA, nS), errmsg.metausat=rep(NA, nS), 
                  pval.metamanova=rep(NA, nS))

for (i in 1:nS){
    print(i)
    
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
    
    res.metausat = metausat(betahat[i,]*sqrt(n), ldscintmat)
    res$pval.metausat[i] = res.metausat$p.metausat
    res$errmsg.metausat[i] = res.metausat$error.msg
    
    res$pval.metamanova[i] = pchisq(n*diag(t(betahat[i,])%*%solve(ldscintmat,betahat[i,])), 
                                    df=nrow(ldscintmat), lower.tail = FALSE)
    
    
    if (i%%200==0) {
        print(i)
        save(res_alltraits, res, file = paste0("res_eff",eff,"_",temp,".rda"))
    }
}


## Assemble results
rm(list=ls())
library(dplyr)

powermat <- matrix(NA, nrow=0, ncol=3, dimnames=list(NULL,c("ASSET","USAT","MANOVA")))
for (eff in c(0.02,0.03)){
    res_alltraits_all = matrix(NA, nrow = 0, ncol = 5)
    
    for (i in 1:50){
        filename <- paste0("res_eff",eff,"_",i,".rda")
        if (file.exists(filename)){
            load(filename)
            if (i==1){
                res.all <- res
            } else{
                res.all <- rbind(res.all, res)
            }
            
            res_alltraits_all = rbind(res_alltraits_all,res_alltraits)
            
            rm(res, res_alltraits)
        }
    }

    dd = data.frame(res_alltraits_all)
    dd = dd %>% filter(!is.na(pval))
    dd = dd %>% mutate(pval.1.meta = 2*pnorm(z1.meta,lower.tail=FALSE), 
                       pval.2.meta = 2*pnorm(z2.meta),
                       pval.1 = ifelse(is.na(pval.1),1,pval.1), 
                       pval.2 = ifelse(is.na(pval.2),1,pval.2))
    
    ddflt = dd %>% mutate(pval.1.meta = ifelse(is.na(pval.1.meta),1,pval.1.meta),
                          pval.2.meta = ifelse(is.na(pval.2.meta),1,pval.2.meta)) %>% 
        filter(pval.1>=pval.1.meta & pval.2>=pval.2.meta)
    
    powermat <- rbind(powermat,
                    cbind(sapply(c(1e-4,1e-5,1e-6,1e-7), function(x) mean(ddflt$pval<x)),
                    sapply(c(1e-4,1e-5,1e-6,1e-7), function(x) mean(res.all$pval.metausat<x,na.rm=T)),
                    sapply(c(1e-4,1e-5,1e-6,1e-7), function(x) mean(res.all$pval.metamanova<x,na.rm=T))))
    rm(res.all,res_alltraits_all,dd,ddflt)
}
write.csv(round(powermat,2), file="powermat.csv")

