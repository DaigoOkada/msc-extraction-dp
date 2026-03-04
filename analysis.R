#Biological analysis using msc extraction
#source("/Users/dokada/Dropbox/analysis/2025.4/msc_extract_anal0302.R") #Run2026.3.4
out_path <- "/Users/dokada/Desktop/work/msc_extract_anal0302/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#Set data_path
data_path <- "/Users/dokada/Desktop/work/msc_extract0225/"

#Pulpitis vs healthy
library(effsize)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(DESeq2)
dp_control <- read.csv(paste0(data_path, "DP_atlas_msc_pseudobulk_counts.csv"), row.names=1)
dp_inflamed <- read.csv(paste0(data_path, "DP_inflamed_msc_pseudobulk_counts.csv"), row.names=1)
common_genes <- intersect(rownames(dp_control), rownames(dp_inflamed))
tmp_count1 <- as.matrix(dp_control[common_genes, ])
tmp_count2 <- as.matrix(dp_inflamed[common_genes, ])
int_counts <- cbind(tmp_count1, tmp_count2)
condition <- c(rep("control", ncol(tmp_count1)), rep("inflamed", ncol(tmp_count2)))

#DEG analysis for inflamed vs control
set.seed(1234)
clst_name <- "inflamed_vs_control"
colData <- data.frame(row.names = colnames(int_counts), condition = factor(condition, levels = c("control","inflamed")))
keep <- which(rowSums(int_counts) > 10)
counts_f <- int_counts[keep,,drop=FALSE]
dds <- DESeqDataSetFromMatrix(countData=counts_f,colData=colData,design=~condition)
dds <- DESeq(dds)
res <- results(dds,contrast=c("condition","inflamed","control"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
pulpitis_deg <- res_df 
write.csv(pulpitis_deg, file = paste0(out_path, clst_name, "_deg_results.csv"))

#GSEA for inflamed vs control
set.seed(1234)
gene_ranks <- res_df$log2FoldChange
names(gene_ranks) <- res_df$gene
gene_ranks <- gene_ranks[is.finite(gene_ranks)]
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
gsea_h <- GSEA(geneList = gene_ranks, TERM2GENE = msig_h[, c("gs_name", "gene_symbol")], pvalueCutoff = Inf, nPermSimple = 10000)
write.csv(gsea_h@result, file = paste0(out_path, clst_name, "_gsea_results.csv"))
#Viz of GSEA
df_inf <- gsea_h@result
df_sig <- df_inf[df_inf$p.adjust < 0.05, ] #26 gene sets
tops <- min(10, nrow(df_sig))
df_top10 <- df_sig[order(abs(df_sig$NES), decreasing=T)[tops:1], ]
df_top10$Description <- factor(df_top10$Description, levels = df_top10$Description)
df_top10$Description <- gsub("^HALLMARK_", "", df_top10$Description)
png(paste0(out_path, clst_name, "_GSEA.png"), width=1400, height=900)
par(mar=c(10,55,10,2), mgp=c(5,2.2,0), tck=-0.02)
barplot(df_top10$NES, horiz=TRUE, names.arg=df_top10$Description, las=1, xlab="", cex.axis=3.5, cex.names=3.0)
mtext("NES", side=1, line=6.5, cex=4)
dev.off()

#Volcano plot
pulpitis_gsea <- df_inf
rownames(pulpitis_gsea) <- pulpitis_gsea$ID
df <- pulpitis_deg[!is.na(pulpitis_deg$padj) & !is.na(pulpitis_deg$log2FoldChange), ]
df$significance <- "NS"
df$significance[df$padj < 0.05 & df$log2FoldChange > 0]  <- "Up"
df$significance[df$padj < 0.05 & df$log2FoldChange < 0] <- "Down"
ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red",
                                "Down" = "blue",
                                "NS" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  theme(
    axis.text  = element_text(size = 20),   # 目盛り
    axis.title = element_text(size = 24),   # 軸ラベル
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    title = ""
  )

ggsave(paste0(out_path, "volcano_plot.png"),width = 8, height = 6, dpi = 300)

#Aging analysi
cpm <- read.csv("/Users/dokada/Desktop/work/msc_extract0225/integrated_msc_cpm.csv", row.names=1)
tmp <- read.csv("/Users/dokada/Desktop/work/msc_extract0225/integrated_sample_ages.csv", row.names=1)
donor_ages <- tmp[,1]
msc_log_expr <- log1p(as.matrix(cpm))
sds <- apply(msc_log_expr, 1, sd) 
means <- apply(msc_log_expr, 1, mean)
msc_expr_sub <- msc_log_expr[means > 1, ]
aging_cor_genes <- rep(NA, nrow(msc_expr_sub))
p_genes <- rep(NA, nrow(msc_expr_sub))
names(aging_cor_genes) <- rownames(msc_expr_sub)
names(p_genes) <- rownames(msc_expr_sub)
for(i in 1:nrow(msc_expr_sub)){
    test_res <- cor.test(msc_expr_sub[i,], donor_ages, method="spearman")
    aging_cor_genes[i] <- test_res$estimate
    p_genes[i] <- test_res$p.value
    #cat(i, "\n")
}
adj_pvals <- p.adjust(p_genes, method="BH") #No genes adjp < 0.05
res <- data.frame("spearman_cor" = aging_cor_genes, "p_value" = p_genes, "adj_p_value" = adj_pvals)
write.csv(res, file=paste0(out_path, "msc_gene_age_correlation.csv"))
#GSEA
set.seed(1234)
clst_name <- "aging"
gene_ranks <- aging_cor_genes
names(gene_ranks) <- names(aging_cor_genes)
gene_ranks <- gene_ranks[is.finite(gene_ranks)]
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
gsea_h <- GSEA(geneList = gene_ranks, TERM2GENE = msig_h[, c("gs_name", "gene_symbol")], pvalueCutoff = Inf, nPermSimple = 10000)
aging_gsea <- gsea_h@result
rownames(aging_gsea) <- aging_gsea$ID
write.csv(aging_gsea, file = paste0(out_path, clst_name, "_gsea_results.csv"))

#Viz of GSEA
df_sig <-aging_gsea[aging_gsea$p.adjust < 0.05, ]
tops <- min(10, nrow(df_sig))
df_top10 <- df_sig[order(abs(df_sig$NES), decreasing=T)[tops:1], ]
df_top10$Description <- factor(df_top10$Description, levels = df_top10$Description)
df_top10$Description <- gsub("^HALLMARK_", "", df_top10$Description)
png(paste0(out_path, clst_name, "_GSEA.png"), width=1400, height=900)
par(mar=c(10,55,10,2), mgp=c(5,2.2,0), tck=-0.02)
barplot(df_top10$NES, horiz=TRUE, names.arg=df_top10$Description, las=1, xlab="", cex.axis=3.5, cex.names=3.0)
mtext("NES", side=1, line=6.5, cex=4)
dev.off()

#Comparison
stopifnot(all(sort(rownames(pulpitis_gsea)) == sort(rownames(aging_gsea))))
genesets <- sort(rownames(pulpitis_gsea))
df1 <- aging_gsea[genesets, ]
df2 <- pulpitis_gsea[genesets, ]
x <- df1$NES
names(x) <- rownames(df1)
y <- df2$NES
names(y) <- rownames(df2)
stopifnot(all(names(x) == names(y)))
x_padj <- df1$p.adjust
y_padj <- df2$p.adjust
compare1 <- "Pulpitis"
compare2 <- "Donor aging"
png(paste0(out_path, compare1, "_vs_", compare2, "_GSEA.png"), width = 960, height = 960)
par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
cols <- rep("black", length(x))
names(cols) <- names(x)
cols[x_padj < 0.05 & y_padj < 0.05] <- "purple"
cols[x_padj < 0.05 & y_padj >= 0.05] <- "red"
cols[x_padj >= 0.05 & y_padj < 0.05] <- "blue"
plot(x, y, xlab="", ylab="", cex.axis=4, cex.lab=4, cex=2, pch=16, col=cols)
xlab <- compare1
ylab <- compare2
mtext(xlab, side = 1, line = 6, cex = 4)
mtext(ylab, side = 2, line = 6, cex = 4)
abline(h=0, lwd=4, col="red")
abline(v=0, lwd=4, col="red")
dev.off()
xy <- data.frame("NES_aging" = x, "NES_pulpitis" = y, "color" = cols)
write.csv(xy, file=paste0(out_path, "GSEA_comparison.csv"))


#save 
sink(paste0(out_path, "session_info.txt")) 
print(sessionInfo()) 
sink()
