############################################################
# RNA-seq data analysis                                    #
############################################################

# Load libraries
library(edgeR)
library(ggplot2)
library(org.Hs.eg.db)
library(pheatmap)
library(dendsort)
library(RColorBrewer)

############################################################
# Differential expression                                  #
############################################################

# Creating groups
group <- c("HT29_WT","HT29_WT","HT29_WT","HT29_WT","HT29_KO","HT29_KO","HT29_KO","HT29_KO")
group <- factor(group)

table(group)

# Loading read counts in edgeR
setwd("~/Dropbox/GU")

GenewiseCounts <- read.delim("readcount.xls", row.names="geneID")

dim(GenewiseCounts)

y <- DGEList(GenewiseCounts, group=group,genes=rownames(GenewiseCounts))

# Adding gene annotation
y$genes$Symbol <- mapIds(org.Hs.eg.db, rownames(y), keytype="ENSEMBL", column="SYMBOL")
y$genes$EntrezID <- mapIds(org.Hs.eg.db, rownames(y), keytype="ENSEMBL", column="ENTREZID")

head(y$gene)

y <- y[!is.na(y$genes$Symbol), ]

dim(y)

# Filtering to remove low counts
keep <- rowSums(cpm(y) > 0.5) >= 4
table(keep)

y <- y[keep, , keep.lib.sizes=FALSE]

# Normalization for composition bias
y <- calcNormFactors(y)
options(digits=3)
y$samples

plotMD(y, column=1)
abline(h=0, col="red", lty=2, lwd=2)

# Exploring differences between libraries
pch <- c(17,16)
colors <- c("orange","purple")
plotMDS(y, col=colors[group], pch=pch[group], cex = 2)
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)
## 800 x 800

# Design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# Dispersion estimation
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

plotQLDisp(fit)

summary(fit$df.prior)

# Testing for differential expression
H.LvsP <- makeContrasts(HT29_KO-HT29_WT, levels=design) ####
res <- glmQLFTest(fit, contrast=H.LvsP)

topTags(res)

is.de <- decideTestsDGE(res)
summary(is.de)

plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")

# Differential expression above a fold-change threshold
tr <- glmTreat(fit, contrast=H.LvsP, lfc=log2(1.5))

is.de <- decideTestsDGE(tr, p.value = 0.05)
summary(is.de)

plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")
## 499+262 genes
## 1200 x 800

# Export differentially expressed genes
DE_F1_A <- as.data.frame(is.de)
DE_F1_B <- subset(DE_F1_A, `1*HT29_KO -1*HT29_WT` != 0)
DE_F1_B$include <- rownames(DE_F1_B)

DE_F1_C <- as.data.frame(topTags(tr,n=Inf))

DE_results <- subset(DE_F1_C, genes %in% DE_F1_B$include)

#write.csv(file="edgeR_DE_results_Treat_761.csv",DE_results,row.names = FALSE)

############################################################
# Gene ontologies goana()                                  #
############################################################

y1 <- y[!duplicated(y$genes$EntrezID),]

rownames(y1) <- y1$genes$EntrezID

fit1 <- glmQLFit(y1, design, robust=TRUE)

H.LvsP2 <- makeContrasts(HT29_KO-HT29_WT, levels=design)

tr1 <- glmTreat(fit1, contrast=H.LvsP2, lfc=log2(1.5))

go <- goana(tr1,species="Hs",FDR=0.05)
go_all <- go[which(go$P.Up<0.01 | go$P.Down<0.01), ]
topGO <- topGO(go_all, n=Inf)

# write.csv(file="topGO_all.csv",topGO,row.names = TRUE)

# Selected GO categories

## 1. GO:0071357 cellular response to type I interferon Up
## 2. GO:0048870 cell motility Down
## 3. GO:0007155 cell adhesion  Down
## 4. GO:0045087 innate immune response  Up
## 5. GO:0009636 response to toxic substance Down
## 6. GO:0048566 embryonic digestive tract development Up
## 7. GO:0001525 angiogenesis Down
## 8. GO:0060429 epithelium development Down
## 9. GO:0045597 positive regulation of cell differentiation Down
## 10. GO:0004872	receptor activity Down
## 11. GO:0071944	cell periphery
## 12. GO:0005576	extracellular region

extract_by_GO <- function(GO_term,DE_table,name){
  hsGO <- org.Hs.egGO2ALLEGS
  Rkeys(hsGO) <- GO_term
  EG <- mappedLkeys(hsGO)
  genes_ENSEMBL = mapIds(org.Hs.eg.db, EG, keytype="ENTREZID", column="ENSEMBL")
  result <- DE_table[DE_table$genes %in% unname(genes_ENSEMBL),]
  result$GO <- name
  return(result)
}

GO_1 <- extract_by_GO("GO:0071357",DE_results,"cellular response to type I interferon") # cellular response to type I interferon Up
GO_2 <- extract_by_GO("GO:0048870",DE_results,"cell motility") # cell motility Down
GO_3 <- extract_by_GO("GO:0007155",DE_results,"cell adhesion") # cell adhesion  Down
GO_4 <- extract_by_GO("GO:0045087",DE_results,"innate immune response") # innate immune response  Up
GO_5 <- extract_by_GO("GO:0009636",DE_results,"response to toxic substance") # response to toxic substance Down
GO_6 <- extract_by_GO("GO:0048566",DE_results,"embryonic digestive tract development") # embryonic digestive tract development Up
GO_7 <- extract_by_GO("GO:0001525",DE_results,"angiogenesis") # angiogenesis Down
GO_8 <- extract_by_GO("GO:0060429",DE_results,"epithelium development") # epithelium development Down
GO_9 <- extract_by_GO("GO:0045597",DE_results,"positive regulation of cell differentiation") # positive regulation of cell differentiation Down
GO_10 <- extract_by_GO("GO:0007186",DE_results,"G-protein coupled receptor signaling pathway") # G-protein coupled receptor signaling pathway

GO_top_genes <- rbind(GO_1,GO_2,GO_3,GO_4,GO_5,GO_6,GO_7,GO_8,GO_9,GO_10)
GO_top_genes
#write.csv(file="GO_top_genes.csv",GO_top_genes,row.names = FALSE)

############################################################
# Plot gene ontologies                                     #
############################################################

# load data
setwd("~/Dropbox/MIEN1/Data")

raw_data_up <- read.csv("sel_GO_up.csv",sep = ",",stringsAsFactors = FALSE)

GO_to_plot <- data.frame("GO" = raw_data_up$GO,"Up" = -log10(raw_data_up$P.Up))
positions <- rev(GO_to_plot$GO)

ggplot(data=GO_to_plot, aes(x=GO, y=Up, fill=GO)) +
  geom_bar(stat="identity", width=0.7,alpha=0.5,fill="red") +
  xlab("GO Term") + ylab("-log10 (Fisher Classic)") +
  theme(plot.title = element_text(color="black", size=24,hjust = 0.5),
        axis.title = element_text(size=24),
        axis.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="none",
        axis.line = element_line(size=0.5))+
  coord_flip() +
  scale_y_continuous(position = "bottom",limits = c(0,10)) +
  scale_x_discrete(limits = positions)

raw_data_down <- read.csv("sel_GO_down.csv",sep = ",",stringsAsFactors = FALSE)

GO_to_plot <- data.frame("GO" = raw_data_down$GO, "Down" = -log10(raw_data_down$P.Down))
positions <- rev(GO_to_plot$GO)

ggplot(data=GO_to_plot, aes(x=GO, y=Down, fill=GO)) +
  geom_bar(stat="identity", width=0.7,alpha=0.5,fill="blue") +
  xlab("GO Term") + ylab("-log10 (Fisher Classic)") +
  theme(plot.title = element_text(color="black", size=24,hjust = 0.5),
        axis.title = element_text(size=24),
        axis.text = element_text(size=22),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="none",
        axis.line = element_line(size=0.5))+
  coord_flip() +
  scale_y_continuous(position = "bottom",limits = c(0,10)) +
  scale_x_discrete(limits = positions)

## 1000 x 600

############################################################
# Heatmap selected genes                                   #
############################################################

# Load data
raw_data <- read.csv("sel_genes.csv",sep = ",",stringsAsFactors = FALSE,row.names = 1)

logCPM <- as.data.frame(cpm(y, prior.count=2, log=TRUE))
logCPM_selected_genes = merge(raw_data,logCPM,by="row.names",sort=FALSE)
logCPM_selected_genes_matrix = logCPM_selected_genes[]
matrix_def = logCPM_selected_genes[,11:18]
row.names(matrix_def) <- raw_data$Symbol

# Dendogram
mat_cluster_cols <- hclust(dist(matrix_def))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")

# Annotations_for_columns
annotation_row = raw_data[,3:7]
rownames(annotation_row) = rownames(matrix_def)
annotation_row

# Colours
Var1 <- c("#00aedb", "white")
names(Var1) <- c("Yes", "No")
Var2 <- c("#ffc425", "white")
names(Var2) <- c("Yes", "No")
Var3 <- c("#d11141", "white")
names(Var3) <- c("Yes", "No")
Var4 <- c("#8874a3", "white")
names(Var4) <- c("Yes", "No")
Var5 <- c("#00b159", "white")
names(Var5) <- c("Yes", "No")

anno_colors <- list(Interferon = Var1,Embryonic=Var2,Angiogenesis=Var3,
                    Cell_adhesion=Var4,Cell_motility=Var5)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}

breaksList = seq(0, 15, by = 1)

pheatmap(mat = matrix_def,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         border_color = "gray", show_colnames = TRUE, show_rownames = TRUE,  drop_levels = TRUE,
         annotation_row = annotation_row,fontsize = 11, cluster_cols = FALSE,cluster_rows = FALSE,
         annotation_colors = anno_colors)
## 800 x 520

############################################################
# Molecular signature MSigDB                               #
############################################################

## Load databases
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c1_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c6_v5p2.rdata"))

# Create indexes
# H hallmark, 1 Positional, 2 curated, 3 motif, 4 computational, 5 GO, 6 oncogenic, 7 immunogenic
idx_H <- ids2indices(Hs.H,id=rownames(y1))
idx_c1 <- ids2indices(Hs.c1,id=rownames(y1))
idx_c2 <- ids2indices(Hs.c2,id=rownames(y1))
idx_c6 <- ids2indices(Hs.c6,id=rownames(y1))

BvsL.v <- makeContrasts(HT29_KO-HT29_WT, levels=design)
cam_H <- camera(y1, idx_H, design, contrast=BvsL.v, inter.gene.cor=0.01)
cam_c1 <- camera(y1, idx_c1, design, contrast=BvsL.v, inter.gene.cor=0.01)
cam_c2 <- camera(y1, idx_c2, design, contrast=BvsL.v, inter.gene.cor=0.01)
cam_c6 <- camera(y1, idx_c6, design, contrast=BvsL.v, inter.gene.cor=0.01)
options(digits=2)

camera_H <- head(cam_H,10)
camera_c1 <- head(cam_c1,10)
camera_c2 <- head(cam_c2,10)
camera_c6 <- head(cam_c6,10)

camera_H
camera_c1
camera_c2
camera_c6

# Plot results
res <- glmQLFTest(fit, contrast=BvsL.v)
barcodeplot(res$table$logFC,index=idx_H[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
            labels=c("HT29 WT","HT29 KO"),main="HALLMARK_INTERFERON_ALPHA_RESPONSE",alpha=1)

barcodeplot(res$table$logFC,index=idx_c2[["RICKMAN_METASTASIS_UP"]],index2 = idx_c2[["RICKMAN_METASTASIS_DN"]],
            labels=c("HT29 WT","HT29 KO"),main="RICKMAN_METASTASIS",alpha=1)

## 800 x 500