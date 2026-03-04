#MSC extract function
#k_can start 5, end 20, step 1
msc_extract <- function(obj, out_path, sample_name,  msc_markers = c("ENG","THY1","NT5E","FRZB","NOTCH3","MCAM"), seed=1, k_can=seq(5,20,by=1)){
    #Create new directory
    set.seed(seed)
    out_dir <- paste0(out_path, sample_name, "/")
    if(file.exists(out_dir)){
        unlink(out_dir, recursive=TRUE)
        dir.create(out_dir)
    }else{
        dir.create(out_dir)
    }
    #Preprocrssing
    DefaultAssay(obj) <- "RNA"
    obj[["percent.mt"]] <- PercentageFeatureSet(obj,pattern="^MT-")
    obj <- subset(obj,subset=nFeature_RNA>200 & percent.mt<10)
    obj <- NormalizeData(obj,normalization.method="LogNormalize",scale.factor=10000)
    obj <- FindVariableFeatures(obj,selection.method="vst",nfeatures=2000)
    obj <- ScaleData(obj,features=VariableFeatures(obj))
    obj <- RunPCA(obj,features=VariableFeatures(obj),npcs=50,verbose=FALSE)
    obj <- RunUMAP(obj,reduction="pca",dims=1:50)
    #Add MSC module score
    gs_use <- intersect(msc_markers,rownames(obj))
    obj <- AddModuleScore(obj,features=list(gs_use),name="MSC")

    #Search pptimal clustering
    pc_mat <- Embeddings(obj,reduction="pca")[,1:50,drop=FALSE]
    opt_score_vec <- c()
    k_vec <- c()
    for(tmp_k in k_can){
        set.seed(tmp_k)
        km <- kmeans(pc_mat, centers = tmp_k, nstart = 100, iter.max = 200)  # iter.max増やす
        cl <- km$cluster
        msc_scores_by_cluster <- tapply(obj$MSC1, cl, mean, na.rm=TRUE)
        top2 <- names(sort(msc_scores_by_cluster, decreasing=TRUE))[1:2]
        x1 <- obj$MSC1[cl == top2[1]]
        x2 <- obj$MSC1[cl == top2[2]]
        opt_score <- cohen.d(x1, x2, paired=FALSE)$estimate
        opt_score_vec <- c(opt_score_vec, opt_score)
    }
    opt_k <- k_can[which.max(opt_score_vec)]

    # get same situation
    set.seed(opt_k)
    opt_km <- kmeans(pc_mat, centers = opt_k, nstart = 100, iter.max = 200)
    opt_cl <- opt_km$cluster
    msc_scores_by_cluster <- tapply(obj$MSC1, opt_cl, mean, na.rm=TRUE)
    top2 <- names(sort(msc_scores_by_cluster, decreasing=TRUE))[1:2]
    x1 <- obj$MSC1[opt_cl == top2[1]]
    x2 <- obj$MSC1[opt_cl == top2[2]]
    opt_score <- cohen.d(x1, x2, paired=FALSE)$estimate
    stopifnot(abs(opt_score - max(opt_score_vec)) < 1e-6) #同じスコアが得られることを確認
    obj$optimal_clusters <- factor(opt_cl)

    # Score trace plot    
    png(paste0(out_dir, sample_name, "_k_selection.png"), width=1200, height=900)
    par(mar=c(9,9,9,2), mgp=c(2.5,2.5,0))
    plot(k_can, opt_score_vec, xlab="", ylab="", cex.axis=4, cex.lab=4, cex=2, pch=16, lwd=4, type="b", col="blue")
    mtext("Number of clusters", side=1, line=6, cex=4)
    mtext("Optimized score", side=2, line=6, cex=4)
    abline(v=opt_k, col="red", lty=4, lwd=4)
    dev.off()
 
    #Annot MSC marker expression
    png(paste0(out_dir, sample_name, ".MSC.boxplot.png"), width = 960, height = 960)
    par(mar = c(9, 9, 9, 2), mgp = c(2.5, 2.5, 0))
    boxplot(obj$MSC1 ~ obj$optimal_clusters, xlab = "", ylab = "", cex.axis = 4, cex.lab = 4, cex = 2)
    xlab <- "Cluster labels"
    ylab <- "MSC Module score"
    main <- paste0("Optimized score: ", round(opt_score, 3))
    mtext(xlab, side = 1, line = 6, cex = 4)
    mtext(ylab, side = 2, line = 6, cex = 4)
    mtext(main, side = 3, line = 3, cex = 4)
    dev.off()
    #Feature plot for obj$MSC1
    p <- FeaturePlot(obj, features = "MSC1", reduction = "umap") +
    ggtitle("MSC module score") +
    theme(
        plot.title = element_text(size=40, hjust=0.5, face="bold"),
        axis.text = element_text(size=34),
        axis.title = element_text(size=38),
        legend.text = element_text(size=20),
        legend.title = element_text(size=34)
    ) +
    xlab("UMAP1") +
    ylab("UMAP2")
    ggsave(paste0(out_dir, sample_name, ".MSC.feature.png"), p, width=12, height=10, dpi=300)
    #Identify MSC cluster
    means <- tapply(obj$MSC1, obj$optimal_clusters, mean, na.rm = TRUE)
    msc_cluster <- names(means)[which.max(means)]
    cells_msc <- colnames(obj)[as.character(obj$optimal_clusters) == msc_cluster]
    obj$msc_or_not <- ifelse(as.character(obj$optimal_clusters) == msc_cluster, "MSC cluster", "Other clusters")
    obj_msc <- subset(obj, cells = cells_msc)
    expr <- as.matrix(GetAssayData(obj_msc, assay="RNA", layer="counts"))
    pseudo_bulk_counts <- rowSums(expr)
    n_msc <- ncol(obj_msc)
    #MSC cluster or not
    highlight_cells <- WhichCells(obj,expression=optimal_clusters==msc_cluster)
    stopifnot(length(highlight_cells) == n_msc)
    p <- DimPlot(obj,reduction="umap",cells.highlight=highlight_cells,cols.highlight="red",cols="grey85",pt.size=1) +
    theme(
        plot.title=element_text(size=40,hjust=0.5,face="bold"),
        axis.text=element_text(size=34),
        axis.title=element_text(size=38),
        legend.text=element_text(size=20),
        legend.title=element_text(size=34),
        legend.position="none"
    ) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    ggtitle(paste0("MSC cluster (n=", n_msc, ")"))
    ggsave(paste0(out_dir, sample_name, ".MSC_highlight_umap.png"),p,width=10,height=8,dpi=300)
    #Save data
    saveRDS(obj, file = paste0(out_dir, sample_name, "_obj.rds"))
    write.csv(pseudo_bulk_counts, file = paste0(out_dir, sample_name, "_pseudo_bulk_counts.csv"))
    result <- list(pseudo_bulk_counts = pseudo_bulk_counts, obj = obj, n_msc = n_msc)
    return(result)
}

