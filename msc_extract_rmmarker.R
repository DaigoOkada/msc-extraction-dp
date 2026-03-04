
#source("/Users/dokada/Dropbox/analysis/2025.4/msc_extract0225_rmmarker.R")
out_path <- "/Users/dokada/Desktop/work/msc_extract0225_rmmarker/"
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
#markers <- c("ENG","THY1","NT5E","FRZB","NOTCH3","MCAM")
markers_FRZB_rm <- c("ENG","THY1","NT5E","NOTCH3","MCAM")
markers_NOTCH3_rm <- c("ENG","THY1","NT5E","FRZB","MCAM")
markers_MCAM_rm <- c("ENG","THY1","NT5E","FRZB","NOTCH3")

#Import sample data: "GSE202476" #13-year-old
tmp_dir <- "/Users/dokada/Desktop/work/shizui_pub_data/GSE202476_RAW/"
data <- Read10X(data.dir=tmp_dir)
obj1 <- CreateSeuratObject(counts=data, min.cells=3)
result <- msc_extract(obj1, msc_markers = markers_FRZB_rm, out_path=out_path, sample_name="GSE202476_FRZB_rm")
result <- msc_extract(obj1, msc_markers = markers_NOTCH3_rm, out_path=out_path, sample_name="GSE202476_NOTCH3_rm")
result <- msc_extract(obj1, msc_markers = markers_MCAM_rm, out_path=out_path, sample_name="GSE202476_MCAM_rm")

#Other dataset: "GSM7029373_Tooth-M" #27 yaers old
counts <- read.table("/Users/dokada/Desktop/work/shizui_pub_data/GSM7029373_Tooth-M.counts.tsv", header=TRUE, row.names=1) #列名重複なし
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"
obj2 <- CreateSeuratObject(counts=counts, min.cells=3)
result <- msc_extract(obj2, msc_markers = markers_FRZB_rm, out_path=out_path, sample_name="GSM7029373_FRZB_rm")
result <- msc_extract(obj2, msc_markers = markers_NOTCH3_rm, out_path=out_path, sample_name="GSM7029373_NOTCH3_rm")
result <- msc_extract(obj2, msc_markers = markers_MCAM_rm, out_path=out_path, sample_name="GSM7029373_MCAM_rm")

#Montage
files <- c(
    "GSE202476_FRZB_rm",
    "GSE202476_NOTCH3_rm",
    "GSE202476_MCAM_rm",
    "GSM7029373_FRZB_rm",
    "GSM7029373_NOTCH3_rm",
    "GSM7029373_MCAM_rm"
)
for(f in files){
    tmp_dir <- paste0(out_path, f, "/")
    setwd(tmp_dir)
    img1 <- paste0(f, ".MSC.feature.png")
    img2 <- paste0(f, "_k_selection.png")
    img3 <- paste0(f, ".MSC.boxplot.png")
    img4 <- paste0(f, ".MSC_highlight_umap.png")
    command <- paste0("montage ", img1, " ", img2, " ", img3, " ", img4, " -resize 1200x1200 -tile 4x1 -geometry +0+0 ", out_path, f, "_combined.png")
    system(command)
}

#save
sink(paste0(out_path, "session_info.txt")) 
print(sessionInfo()) 
sink()
