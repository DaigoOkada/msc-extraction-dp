#source("/Users/dokada/Dropbox/analysis/2025.4/msc_extract0225.R") #2026.3.1
out_path <- "/Users/dokada/Desktop/work/msc_extract0225/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#library
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(effsize)

#set functions
source("/Users/dokada/Dropbox/analysis/2025.4/msc_extract_function3.R")

#Set marker
markers <- c("ENG","THY1","NT5E","FRZB","NOTCH3","MCAM")
marker3 <- c("ENG","THY1","NT5E")

#Import sample data: "GSE202476" #13-year-old
tmp_dir <- "/Users/dokada/Desktop/work/shizui_pub_data/GSE202476_RAW/"
data <- Read10X(data.dir=tmp_dir)
obj1 <- CreateSeuratObject(counts=data, min.cells=3)
result <- msc_extract(obj1, msc_markers = markers, out_path=out_path, sample_name="example1_markers")
pseudo_bulk_counts_ex1 <- result$pseudo_bulk_counts
obj1_6m <- result$obj
result <- msc_extract(obj1, msc_markers = marker3, out_path=out_path, sample_name="example1_marker3")
obj1_3m <- result$obj

#Other dataset: "GSM7029373_Tooth-M" #27 yaers old
counts <- read.table("/Users/dokada/Desktop/work/shizui_pub_data/GSM7029373_Tooth-M.counts.tsv", header=TRUE, row.names=1) #列名重複なし
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"
obj2 <- CreateSeuratObject(counts=counts, min.cells=3)
result <- msc_extract(obj2, msc_markers = markers, out_path=out_path, sample_name="example2_markers")
pseudo_bulk_counts_ex2 <- result$pseudo_bulk_counts
obj2_6m <- result$obj
result <- msc_extract(obj2, msc_markers = marker3, out_path=out_path, sample_name="example2_marker3")
obj2_3m <- result$obj



#Human dental pulp Atlas (Healthy controls)
data_path <- "/Users/dokada/Desktop/work/shizui_data/GSE164157_RAW/"
sample_names <- c("Pulp1", "Pulp2", "Pulp3", "Pulp4", "Pulp5")
gsms <- c("GSM4998457", "GSM4998458", "GSM4998459", "GSM4998460", "GSM4998461")
expr_list <- list()
common_genes <- c()
n_cells_vec <- c()
for(i in 1:length(sample_names)){
    tmp_dir <- paste0(data_path, sample_names[i], "/")
    setwd(tmp_dir)
    old_name1 <- paste0(gsms[i],"_", sample_names[i], "_barcodes.tsv.gz")
    if (file.exists(old_name1)) file.rename(old_name1, "barcodes.tsv.gz")
    old_name2 <- paste0(gsms[i],"_", sample_names[i], "_genes.tsv.gz")
    if (file.exists(old_name2)) file.rename(old_name2, "features.tsv.gz")
    old_name3 <- paste0(gsms[i],"_", sample_names[i], "_matrix.mtx.gz")
    if (file.exists(old_name3)) file.rename(old_name3, "matrix.mtx.gz")
    #Analysis
    dat <- Read10X(data.dir = tmp_dir)
    obj <- CreateSeuratObject(counts = dat)
    tmp <- msc_extract(obj, msc_markers = markers, out_path = out_path, sample_name = sample_names[i])
    pseudo_bulk_counts <- tmp$pseudo_bulk_counts
    n_cells <- tmp$n_msc
    n_cells_vec <- c(n_cells_vec, n_cells)
    expr_list[[sample_names[i]]] <- pseudo_bulk_counts
    #Common genes
    if(i == 1){
        common_genes <- names(pseudo_bulk_counts)
    }else{
        common_genes <- intersect(common_genes, names(pseudo_bulk_counts))
    }
    cat(sample_names[i], "done\n")
}
mat <- cbind(expr_list[["Pulp1"]][common_genes], expr_list[["Pulp2"]][common_genes], expr_list[["Pulp3"]][common_genes], expr_list[["Pulp4"]][common_genes], expr_list[["Pulp5"]][common_genes])
colnames(mat) <- sample_names
write.table(data.frame("sample"=sample_names, "n_msc_cells"=n_cells_vec), file=paste0(out_path, "DP_atlas_msc_n_cells.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.csv(mat, file=paste0(out_path, "DP_atlas_msc_pseudobulk_counts.csv"))


#Pulpitis
data_path <- "/Users/dokada/Desktop/work/shizui_data/GSE280528_RAW/"
sample_names <- c("Pulpitis1", "Pulpitis2", "Pulpitis3")
gsms <- c("GSM8599816", "GSM8599817", "GSM8599818")
expr_list <- list()
common_genes <- c()
n_cells_vec <- c()
for(i in 1:length(sample_names)){
    tmp_dir <- paste0(data_path, sample_names[i], "/")
    setwd(tmp_dir)
    old_name1 <- paste0(gsms[i],"_", sample_names[i], "_barcodes.tsv.gz")
    if (file.exists(old_name1)) file.rename(old_name1, "barcodes.tsv.gz")
    old_name2 <- paste0(gsms[i],"_", sample_names[i], "_features.tsv.gz")
    if (file.exists(old_name2)) file.rename(old_name2, "features.tsv.gz")
    old_name3 <- paste0(gsms[i],"_", sample_names[i], "_matrix.mtx.gz")
    if (file.exists(old_name3)) file.rename(old_name3, "matrix.mtx.gz")
    #Analysis
    dat <- Read10X(data.dir = tmp_dir)
    obj <- CreateSeuratObject(counts = dat)
    tmp <- msc_extract(obj, msc_markers = markers, out_path = out_path, sample_name = sample_names[i])
    pseudo_bulk_counts <- tmp$pseudo_bulk_counts
    n_cells <- tmp$n_msc
    n_cells_vec <- c(n_cells_vec, n_cells)
    expr_list[[sample_names[i]]] <- pseudo_bulk_counts
    #Common genes
    if(i == 1){
        common_genes <- names(pseudo_bulk_counts)
    }else{
        common_genes <- intersect(common_genes, names(pseudo_bulk_counts))
    }
    cat(sample_names[i], "done\n")
}
mat <- cbind(expr_list[["Pulpitis1"]][common_genes], expr_list[["Pulpitis2"]][common_genes], expr_list[["Pulpitis3"]][common_genes])
colnames(mat) <- sample_names
write.csv(mat, file=paste0(out_path, "DP_inflamed_msc_pseudobulk_counts.csv"))
write.table(data.frame("sample"=sample_names, "n_msc_cells"=n_cells_vec), file=paste0(out_path, "DP_inflamed_msc_n_cells.txt"), sep="\t", row.names=FALSE, quote=FALSE)

#Donor aging dataset
ages <- c("18y", "22y", "26y", "35y", "45y", ">48y")
age_values <- c(18, 22, 26, 35, 45, 48)
samples <- c("Control2", "Control1", "Control3", "Control6", "Control5", "Control7")
geo <- c("GSM8452929", "GSM8452930", "GSM8452931", "GSM8452932", "GSM8452933", "GSM8452934")
data_dirs <- paste0("/Users/dokada/Desktop/work/detarmin_data/GSE274562_RAW/", samples, "/")
n_samples <- length(data_dirs)
n_cells_vec <- c()
msc_expr <- NULL
for(i in 1:n_samples){
    tmp_dir <- data_dirs[i]
    data <- Read10X(data.dir = tmp_dir)
    obj <- CreateSeuratObject(counts = data)
    #analysis
    tmp <- msc_extract(obj, out_path, samples[i])
    pseudo_bulk_counts <- tmp$pseudo_bulk_counts
    n_cells <- tmp$n_msc
    n_cells_vec <- c(n_cells_vec, n_cells)
    msc_expr <- cbind(msc_expr, pseudo_bulk_counts)
    cat(i, "done\n")
}
annot <- data.frame("sample" = samples, "age" = age_values, "geo" = geo)
colnames(msc_expr) <- samples
write.csv(annot, file=paste0(out_path, "sample_annotation.csv"))
write.csv(msc_expr, file=paste0(out_path, "aging_msc_pseudobulk_counts.csv"))


#Integrated count matrix
counts_aging <- as.matrix(msc_expr)
common_genes <- intersect(intersect(names(pseudo_bulk_counts_ex1), names(pseudo_bulk_counts_ex2)), rownames(counts_aging))
msc_aging_counts <- cbind(pseudo_bulk_counts_ex1[common_genes], pseudo_bulk_counts_ex2[common_genes], counts_aging[common_genes, ])
ann <- read.csv(paste0(out_path,"sample_annotation.csv"),row.names=1,check.names=FALSE)
age1 <- c(13, 27, ann$age)
msc_expr_cpm <- msc_aging_counts
for(i in 1:ncol(msc_expr_cpm)){
    total_counts <- sum(msc_expr_cpm[,i ])
    msc_expr_cpm[, i] <- (msc_expr_cpm[, i] / total_counts) * 1e6
}
write.csv(msc_expr_cpm, file=paste0(out_path, "integrated_msc_cpm.csv"))
write.csv(age1, file=paste0(out_path, "integrated_sample_ages.csv"))

#save 
sink(paste0(out_path, "session_info.txt")) 
print(sessionInfo()) 
sink()
