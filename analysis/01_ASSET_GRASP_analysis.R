## Use fast ASSET to analysis summary statistics fo 116 traits (mostly from GRASP repository)
rm(list = ls())
library(data.table)
library(ASSET)
source("../support_files/ASSET_block_orth.R")

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

temp = as.integer(commandArgs(trailingOnly = TRUE)) # Indices of small batchs of summary stats
print(temp)
split = temp[1] # 1:18
grp = temp[2] # 1:10
load(paste0("PATH_TO_DATA/data_split",split,".rda"))

snpind = (50000*(grp-1)+1):min(50000*grp,nrow(data_split$zstat))
nsnp = length(snpind) # Number of SNPs
data_split$zstat = data_split$zstat[snpind,] # Remove blood cell traits - this group is too big and needs to be dealt with separately
data_split$sampsize = data_split$sampsize[snpind,]
rm(snpind)

tempvec = rep(NA,nsnp)
names(tempvec) = rownames(data_split$zstat)
pval_ast = error_warning = tempvec
ast1s_info = matrix(NA, nrow = nsnp, ncol = 8) # Information of two one-sided searches
colnames(ast1s_info) = c("pval.1","beta.1","sd.1","sd.1.meta","pval.2","beta.2","sd.2","sd.2.meta")
pheno_pos = pheno_neg = beta_scr = sigma_scr = vector("list", length=nsnp) # pheno_scr are the traits that passed pre-screening
names(pheno_pos) = names(pheno_neg) = names(beta_scr) = rownames(data_split$zstat)

## Two-sided ASSET
side = 2
for (i in 1:nsnp){
    SNP = rownames(data_split$zstat)[i]
    zstat_i = data_split$zstat[i,]
    sampsize_i = data_split$sampsize[i,]
    
    ind = (!is.na(zstat_i)) & (!is.na(sampsize_i))
    # Pre-screening is only implemented when data for >=2 traits are available
    if (sum(ind)>=2){
        zstat_i = zstat_i[ind]
        sampsize_i = sampsize_i[ind]
        
        traits_i = names(zstat_i)
        betahat = zstat_i/sqrt(sampsize_i)
        Neff = sampsize_i
        SE = 1/sqrt(Neff)
        
        # blockwise orthogonalization
        dt_scr = scr.transform.block(betahat, SE, SNP, traits_i, cor = ldscintmat, block, scr_pthr = 0.05)
        if (sum(dt_scr$betahat_orth>0)>16 | sum(dt_scr$betahat_orth<=0)>16){
            error_warning[i] = 6 # Many traits pass screening
        }
        # 16 is the genome-wide significance threshold for binomial probabilities
        # > pbinom(16,116,0.025,lower.tail=FALSE)
        # [1] 5.640032e-09
        
        betahat_orth = dt_scr$betahat_orth
        sigma_orth = dt_scr$sigma_orth
        traits_i=names(betahat_orth)
        beta_scr[[i]] = betahat_orth
        sigma_scr[[i]] = sigma_orth
        rm(dt_scr)
        
        if (length(traits_i)>=2 & is.na(error_warning[i])){
            res.asset = try(h.traits(snp.vars=SNP, traits.lab=traits_i,
                                     beta.hat=betahat_orth, sigma.hat=sigma_orth, ncase=2/sigma_orth^2,                                    
                                     ncntl=2/sigma_orth^2, cor=ldscintmat[traits_i,traits_i], side=side, search=side, cor.numr = FALSE))
            if (class(res.asset)!="try-error"){
                pval_ast[i] = res.asset$Subset.2sided$pval
                pheno_pos[[i]] = colnames(res.asset$Subset.2sided$pheno.1)[res.asset$Subset.2sided$pheno.1]
                pheno_neg[[i]] = colnames(res.asset$Subset.2sided$pheno.2)[res.asset$Subset.2sided$pheno.2]
                error_warning[i] = 0 # ASSET implemented without error
                
                ast1s_info[i,"pval.1"] = res.asset$Subset.2sided$pval.1
                ast1s_info[i,"beta.1"] = res.asset$Subset.2sided$beta.1
                ast1s_info[i,"sd.1"] = res.asset$Subset.2sided$sd.1
                ast1s_info[i,"sd.1.meta"] = res.asset$Subset.2sided$sd.1.meta
                ast1s_info[i,"pval.2"] = res.asset$Subset.2sided$pval.2
                ast1s_info[i,"beta.2"] = res.asset$Subset.2sided$beta.2
                ast1s_info[i,"sd.2"] = res.asset$Subset.2sided$sd.2
                ast1s_info[i,"sd.2.meta"] = res.asset$Subset.2sided$sd.2.meta
                
            } else{
                error_warning[i] = 1 # Error when implementing ASSET
                print(paste(i,"SNP",SNP,"error;",length(traits_i),"traits pass threshold"))
            }
        } else if (length(traits_i)==1){
            error_warning[i] = 2 # Only one trait is left after screening
            pval_ast[i] = as.numeric(2*pnorm(-abs(betahat_orth/sigma_orth)))
            if (betahat_orth>0){
                pheno_pos[[i]] = traits_i
                pheno_neg[[i]] = character(0)
            } else{
                pheno_pos[[i]] = character(0)
                pheno_neg[[i]] = traits_i
            }
        } else if (length(traits_i)==0){
            error_warning[i] = 3 # No trait left after screening
            pheno_pos[[i]] = pheno_neg[[i]] = character(0)
        }
    } else if (sum(ind)==1){
        error_warning[i] = 4 # Only one trait is not NA before screening
        pval_ast[i] =  as.numeric(2*pnorm(-abs(zstat_i[ind])))
        if (zstat_i[ind]>0){
            pheno_pos[[i]] = colnames(ldscintmat)[ind]
            pheno_neg[[i]] = character(0)
        } else{
            pheno_pos[[i]] = character(0)
            pheno_neg[[i]] = colnames(ldscintmat)[ind]
        }
    } else{
        error_warning[i] = 5 # Data missing for all traits
        pheno_pos[[i]] = pheno_neg[[i]] = character(0)
    }
    
    if (i%%1000==0){
        print(i)
        res = list(pval_ast = pval_ast, ast1s_info = ast1s_info,
                   beta_scr = beta_scr, sigma_scr = sigma_scr,
                   pheno_pos = pheno_pos, pheno_neg = pheno_neg, error_warning = error_warning)
        save(res, file = paste0("res_split",split,"_",grp,".rda"))
    }
}    
res = list(pval_ast = pval_ast, ast1s_info = ast1s_info,
           beta_scr = beta_scr, sigma_scr = sigma_scr,
           pheno_pos = pheno_pos, pheno_neg = pheno_neg, error_warning = error_warning)
save(res, file = paste0("res_split",split,"_",grp,".rda"))
