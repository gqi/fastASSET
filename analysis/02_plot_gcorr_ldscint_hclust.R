rm(list=ls())
library(dplyr)
library(data.table)
library(readxl)
library(corrplot)

traitdf = read_excel("../support_files/supptab1_final_data_list_w_DISPLAYNAME.xlsx")

# Load LDSC intercept matrix: correlation of z stats
load("../support_files/ldsc_log/ldscintmat.rda")
ldscintmat = ldscintmat[traitdf$trait,traitdf$trait]
load("../support_files/gcorrmat_lowh2higcorr_rm.rda")
gcorrmat_lowh2higcorr_rm = gcorrmat_lowh2higcorr_rm[traitdf$trait,traitdf$trait]
all.equal(rownames(ldscintmat),rownames(gcorrmat_lowh2higcorr_rm))
all.equal(rownames(ldscintmat),traitdf$trait)

gcorr <- gcorrmat_lowh2higcorr_rm
gcorr[upper.tri(gcorr)] <- ldscintmat[upper.tri(gcorr)]
gcorr[is.na(gcorr)] <- 0
gcorr[gcorr>1] <- 1
colnames(gcorr) <- rownames(gcorr) <- traitdf$Display_name
col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", 
                          "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", 
                          "#4393C3", "#2166AC", "#053061"))(200)
col <- col[length(col):1]
png(filename = "gcorr_ldscint.png", width = 3000, height = 3000, res = 300, type = "cairo")
corrplot(gcorr, method="color", tl.col="black", tl.cex = 0.4, col = col)
dev.off()

# Hierarchical clustering plot
ldscintmat = solve(diag(sqrt(diag(ldscintmat)))) %*% ldscintmat %*% solve(diag(sqrt(diag(ldscintmat))))
colnames(ldscintmat) = rownames(ldscintmat) = traitdf$Display_name
# Hierarchical clustering
corrdist = as.dist(1-abs(ldscintmat))
hc = hclust(corrdist)

png(filename = "hclust.png", width = 3500, height = 1500, res = 300, type = "cairo")
plot(hc, cex = 0.6, ylab = "1-|correlation of z statistics|", 
     xlab = "Trait", main = "")
abline(h=0.8, lty = 2, col=2, cex = 1)
dev.off()
