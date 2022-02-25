### Below is orthogonalization function for a single SNP
# Only implement two-side orthogonalization - it is unlikely for the SNP to have effects on many traits in the same direction
## beta.hat, sigma.hat, n.eff, ncase, ncntl are matrices
scr.transform.singlesnp <- function(snp.vars, traits.lab, beta.hat, sigma.hat, cor=NULL, n.eff=NULL, p.bound=1)
{
    k <- length(traits.lab)
    if((mode(beta.hat)!="numeric") || (length(beta.hat) != k)) 
        stop("Error beta.hat should be numeric vector of length #trait")
    if((mode(sigma.hat)!="numeric") || (length(sigma.hat) != k)) 
        stop("Error sigma.hat should be numeric vector of length #trait")
    names(beta.hat) <- traits.lab
    names(sigma.hat) <- traits.lab
    if(!is.null(n.eff))
    {
        ord <- order(n.eff)
        beta.mat <- beta.hat[ord]
        sigma.mat <- sigma.hat[ord]
        t.lab <- traits.lab[ord]
    }	
    if(p.bound < 1){
        c0 <- qnorm(p.bound/2, lower.tail=FALSE)
    } else {
        return(list(beta.mat=beta.mat, sigma.mat=sigma.mat))
    }
    
    if(!is.null(cor))
    {
        cor.ord = cor[ord,ord]
        U <- chol(cor.ord)
        Ui <- solve(U)
    } else {
        Ui <- diag(1, k)
    }
    zz <- beta.mat/sigma.mat
    zz.o <- zz %*% Ui

    # Orthogonalization for two-sided ASSET
    # pp1: P(|Z|>z | |Z|>c0)
    pp1 <- pmax(pmin(2 * pnorm(abs(zz.o), lower.tail=FALSE)/p.bound, 0.9999999), .Machine$double.xmin)
    scr <- (abs(zz.o) > c0)
    zz0 = sign(zz.o)*qnorm(pp1/2, lower.tail=FALSE)
    
    colnames(scr) <- t.lab
    
    zzn <- zz0 %*% U
    
    beta1 <- as.vector(zzn * sigma.mat)
    names(beta1) <- t.lab
    
    beta.scr <- beta1[scr]
    sigma.scr <- sigma.mat[scr]
    
    return(list(beta.scr=beta.scr, sigma.scr=sigma.scr))
}

### Blocked diagonal pattern
### Check if the non missing traits has block structures (estimated below)
### If so, do blockwise orthogonalization
scr.transform.block = function(betahat, SE, SNP, traits_i, cor, block, scr_pthr = 0.05){
    traits_i_block = sapply(block, function(x) traits_i[traits_i%in%x])
    traits_i_block = traits_i_block[sapply(traits_i_block,length)>=2]
    
    betahat_orth = vector(length = 0)
    sigma_orth = vector(length = 0)
    
    if (length(traits_i_block)>0){
        for (j in 1:length(traits_i_block)){
            scr = scr.transform.singlesnp(SNP, traits_i_block[[j]], beta.hat = betahat[traits_i_block[[j]]],
                                          sigma.hat = SE[traits_i_block[[j]]], n.eff = Neff[traits_i_block[[j]]],
                                          cor = cor[traits_i_block[[j]],traits_i_block[[j]]], p.bound = scr_pthr)
            # print(j)
            # print(scr)
            betahat_orth = c(betahat_orth, scr$beta.scr)
            sigma_orth = c(sigma_orth, scr$sigma.scr)
        }
    }
    
    betahat = betahat[!(names(betahat)%in%unlist(traits_i_block))]
    SE = SE[!(names(SE)%in%unlist(traits_i_block))]

    # two-sided correction
    scr = (2*pnorm(abs(betahat/SE),lower.tail=FALSE))<scr_pthr
    betahat = betahat[scr]
    SE = SE[scr]
    pp1 = pmin(2*pnorm(abs(betahat/SE),lower.tail=FALSE)/scr_pthr, 0.9999999)
    betahat_orth = c(betahat_orth, SE*sign(betahat)*abs(qnorm(pp1/2,lower.tail=FALSE)))
    sigma_orth = c(sigma_orth, SE)
    
    traits_i=names(betahat_orth)
    
    return(list(betahat_orth = betahat_orth, sigma_orth = sigma_orth))
}

