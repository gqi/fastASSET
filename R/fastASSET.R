scr_transform_in_block <- function(snp, traits.lab, beta.hat, sigma.hat, cor=NULL, n.eff=NULL, p.bound=1)
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

scr_transform = function(betahat, SE, SNP, traits_i, cor, block, Neff, scr_pthr = 0.05){
    traits_i_block = lapply(block, function(x) traits_i[traits_i%in%x])
    traits_i_block = traits_i_block[sapply(traits_i_block,length)>=2]

    betahat_orth = vector(length = 0)
    sigma_orth = vector(length = 0)

    if (length(traits_i_block)>0){
        for (j in 1:length(traits_i_block)){
            scr = scr_transform_in_block(SNP, traits_i_block[[j]], beta.hat = betahat[traits_i_block[[j]]],
                                          sigma.hat = SE[traits_i_block[[j]]], n.eff = Neff[traits_i_block[[j]]],
                                          cor = cor[traits_i_block[[j]],traits_i_block[[j]]], p.bound = scr_pthr)
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

#' Fast ASSET using pre-screening
#' @description This function implements a fast version of ASSET for analysis of a large number of traits. It accelerates the computation by restricting subset search among the traits with suggestive level of associations. The traits are selected by a liberal p-value threshold. The input is GWAS summary statistics for one SNP. See details for more information.
#' @param snp A character string giving the SNP name to be analyzed. No default.
#' @param traits.lab A character vector giving the names/identifiers of the k studies/traits being analyzed. The order of this vector must match the columns of beta.hat and sigma.hat No default.
#' @param beta.hat A numeric vector of length k giving the coefficients obtained from the analysis of that SNP across the k studies/traits. No default.
#' @param sigma.hat A vector of same dimension as beta.hat, giving the corresponding standard errors. No default.
#' @param Neff A numeric vector of same dimension as beta.hat, giving the effective sample size of each study. For continuous traits, the effective sample size is the total sample size; for binary traits, the effective sample size is Ncase*Ncontrol/(Ncase+Ncontrol).
#' @param cor A matrix of dimension k by k storing the correlation of z-statistics under the null hypothesis. It can be estimated using LD score regression (https://github.com/bulik/ldsc). For example, element \code{(j,k)} is the LD score regression intercept between traits \code{j} and {k}. Diagonal elements are equal to 1. Set to identity matrix if different studies do not have overlapping samples.
#' @param block A list of character vectors. Each element is a vector of study/trait names that fall into the same correlation block based on \code{cor}. Studies from different blocks are treated as independent. Studies that are not in \code{block} are treated as independent from other studies.
#' @param scr_pthr P-value threshold for pre-screening. Only traits with p-value < scr_pthr are retained. Default to 0.05.
#'
#' @details The standard ASSET (\code{h.traits}) which searches through all subsets can be computationally intractable for analyzing a large number of traits. \code{fast_asset} reduces the computational burden by the following procedure: (1) De-correlate the z-statistics associated with different traits using the correlation matrix (\code{cor}); (2) select the set of traits using a liberal threshold (p<0.05 by default) of the de-correlated z-statistics; (3) adjust the de-correlated z-statistics for the pre-selection independently across different traits; (4) re-introduce the correlation using the matrix \code{cor} and supply the adjusted statistics to standard ASSET.
#'
#' @return A list of the same format as the output of \code{h.traits}.
#'
#' @references Guanghao Qi, Surya Chhetri, Elizaveta Naydanova, Debashree Ray, Diptavo Dutta, Alexis Battle, Samsiddhi Bhattacharjee and Nilanjan Chatterjee. "Genome-wide Multi-trait Analysis Across 116 Studies Identifies Complex Patterns of Pleiotropy and Unique Trait-Specific Loci". In preparation (2021).
#'
#' @examples
#' data("example_rs6678982", package="fastASSET")
#' block <- create_blocks(ldscintmat)
#' test <- fast_asset(snp=SNP, traits.lab=traits, beta.hat=betahat, sigma.hat=SE,
#' Neff=Neff, cor=ldscintmat, block=block, scr_pthr=0.05)
#' @import ASSET
#' @export
fast_asset <- function(snp, traits.lab, beta.hat, sigma.hat, Neff, cor, block, scr_pthr=0.05){
    ind <- (!is.na(beta.hat)) & (!is.na(sigma.hat))
    beta.hat <- beta.hat[ind]
    sigma.hat <- sigma.hat[ind]
    Neff <- Neff[ind]
    traits_i <- traits.lab[ind]
    names(beta.hat) <- names(sigma.hat) <- names(Neff) <- traits_i

    # Pre-screening is only implemented when data for >=2 traits are available
    if (sum(ind)>=2){
        cormat <- cor[ind,ind]
        rownames(cormat) <- colnames(cormat) <- traits_i

        # Pre-screen and adjust the summary statistics
        dt_scr <- scr_transform(betahat=beta.hat, SE=sigma.hat, SNP=snp, traits_i=traits_i,
                                cor = cormat, block, Neff=Neff, scr_pthr = 0.05)
        if (sum(dt_scr$betahat_orth>0)>16 | sum(dt_scr$betahat_orth<=0)>16){
            print("Too many traits passed pre-screening, subset search can be slow. Consider using a lower scr_pthr.")
        }

        betahat_orth <- dt_scr$betahat_orth
        sigma_orth <- dt_scr$sigma_orth
        traits_i <- names(betahat_orth)
        rm(dt_scr)

        if (length(traits_i)>=2){
            res.asset <- h.traits(snp.vars=snp, traits.lab=traits_i,
                                 beta.hat=betahat_orth, sigma.hat=sigma_orth, ncase=2/sigma_orth^2,
                                 ncntl=2/sigma_orth^2, cor=cormat[traits_i,traits_i],
                                 side=2, search=2, cor.numr = FALSE)
        } else if (length(traits_i)==1){
            res.asset <- h.traits(snp.vars=snp, traits.lab=traits_i,
                                 beta.hat=betahat_orth, sigma.hat=sigma_orth, ncase=2/sigma_orth^2,
                                 ncntl=2/sigma_orth^2, side=2, search=2)
        } else if (length(traits_i)==0){
            stop("No trait passed pre-screening. Use a higher scr_pthr or standard ASSET.")
        }
    } else if (sum(ind)==1){
        res.asset <- h.traits(snp.vars=SNP, traits.lab=traits_i,
                             beta.hat=beta.hat, sigma.hat=sigma.hat, ncase=2/sigma.hat^2,
                             ncntl=2/sigma.hat^2, side=2, search=2)
    } else{
        stop("Data are missing for every trait.")
    }

    return(res.asset)
}

#' Create blocks from correlation matrix
#' @description Partition the set of traits into blocks of correlated traits using hierarchical clustering.
#' @param cormat A named k by k correlation matrix of z-statistics with row and column names being the names of studies/traits.
#' @param cor_thr Threshold of correlation to divide the set of traits into blocks. Z-statistics with correlation less than \code{cor_thr} are treated as uncorrelated.
#' @return A list of character vectors. Each element is a vector of study/trait names that fall into the same correlation block based on \code{cormat}. Studies from different blocks are treated as independent. Studies that are not in the list are treated as independent from other studies.
#' @examples
#' data("example_rs6678982", package="fastASSET")
#' block <- create_blocks(ldscintmat)
#' @export
create_blocks <- function(cormat, cor_thr=0.2){
    # Hierarchical clustering
    corrdist <- as.dist(1-abs(cormat))
    hc <- hclust(corrdist)
    htree <- cutree(hc, h=1-cor_thr) # Criterion: corr>0.2
    block <- as.integer(names(table(htree))[table(htree)>=2])
    block <- lapply(block, function(x) names(htree)[htree==x])

    return(block)
}


