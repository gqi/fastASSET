#' dag: data frame of at least 4 columns: SNP, gene, tissue, slope
eqtl_connection_plot <- function(dag1t, gtex.colcode,
                                 hflen.snp=5.5, hflen.gene=7, hfht=1, slope_ratio_thr=0.7,
                                 filename = "eqtl_combined.png", height = 4000, width = 4100, res = 300,
                                 trtcol=c(brewer.pal(n=8, name = "Set2"),"steelblue","cyan","pink")){
    gene.numtissue <- dag1t %>% group_by(SNP,gene) %>% summarize(numtissue = length(tissue))
    dag1t <- dag1t %>%
        group_by(gene) %>% 
        filter(abs(slope)>slope_ratio_thr*max(abs(slope))) %>%
        mutate(edge_col = ifelse(abs(slope)==max(abs(slope)),"grey30","grey60"),
               edge_lwd = ifelse(abs(slope)==max(abs(slope)),1,0.5),
               edge_lty = ifelse(abs(slope)==max(abs(slope)),2,1)) %>%
        ungroup() 
    
    dag1t.snpgene <- dag1t %>% select(SNP,gene,pheno,colgrp) %>% distinct() %>%
        left_join(gene.numtissue, by=c("SNP","gene"))
    
    # trtcol <- brewer.pal(n=8, name = "Set2")
    names(trtcol) <- unique(dag1t.snpgene$colgrp)
    dag1t.snpgene$col <- trtcol[dag1t.snpgene$colgrp]
    dag1t.snp <- dag1t.snpgene %>% select(SNP,pheno,col) %>% distinct()
    
    y.gene <- seq(5,95,length.out = length(dag1t.snpgene$col))
    y.gene <- y.gene[length(y.gene):1]
    names(y.gene) <- dag1t.snpgene$gene
    
    dag1t.snpgene$y <- y.gene
    temp <- dag1t.snpgene %>% group_by(SNP) %>% summarize(y.snp = mean(y))
    dag1t.snp$y <- temp$y.snp[match(dag1t.snp$SNP,temp$SNP)]
    y.snp <- dag1t.snp$y
    names(y.snp) <- dag1t.snp$SNP
    
    dag1t.pheno <- dag1t.snp %>% group_by(pheno) %>% summarize(y=mean(y), col=unique(col))
    
    alltissue <- unique(dag1t$tissue)
    y.tissue <- seq(5,95,length.out = length(alltissue))
    y.tissue <- y.tissue[length(y.tissue):1]
    names(y.tissue) <- alltissue
    
    png(filename = filename, height = height, width = width, res = res)
    
    plot(1, type="n", xlab="", ylab="", xlim=c(-2,107), ylim=c(-2,107),
         bty = 'n', xaxt = 'n', yaxt = 'n')
    rect(xleft=20-hflen.snp,xright=20+hflen.snp, ytop=y.snp+hfht, ybottom = y.snp-hfht, col = dag1t.snp$col, border=NA)
    rect(xleft=50-hflen.gene,xright=50+hflen.gene, ytop=y.gene+hfht, ybottom = y.gene-hfht, col = dag1t.snpgene$col, border=NA)
    points(x=rep(80,length(y.tissue)), y=y.tissue, 
           col=gtex.colcode$col[match(names(y.tissue),gtex.colcode$tsname_indata)], pch=16, cex = 1.5)
    
    text(x=20, y=y.snp, labels=names(y.snp), cex=0.8, font=2)
    text(x=50, y=y.gene, labels=paste0(dag1t.snpgene$gene," (",dag1t.snpgene$numtissue,")"), cex=0.8, font=4)
    text(x=83, y=y.tissue, labels=gtex.colcode$tsname[match(names(y.tissue),gtex.colcode$tsname_indata)], cex=0.9, adj=0)
    segments(x0=20+hflen.snp, y0=y.snp[dag1t.snpgene$SNP], x1=50-hflen.gene, y1=y.gene[dag1t.snpgene$gene],
             lwd = 0.8, col="grey50")
    segments(x0=50+hflen.gene, y0 = y.gene[dag1t$gene], x1=79, y1=y.tissue[dag1t$tissue],
             lwd = dag1t$edge_lwd, col=dag1t$edge_col, lty=dag1t$edge_lty)
    text(x=c(20,50,85), y=99, labels=c("SNP","Gene","Tissue"), font=2, cex=1.5)
    text(x=2, y=dag1t.pheno$y, labels=dag1t.pheno$pheno, col=dag1t.pheno$col, font=2, cex=1.1)
    
    dev.off()
}
