library(Matrix)
library(Seurat)
library(ggplot2)
library(biomaRt)
library(sctransform)
library(svglite)
library(RColorBrewer)

# # check UMI and percent.mito plots, write barcodes to file for citeseq count input
# labels <- c('B1','B2','B3','B4','B5A','B5B','B6A','B6B','B7A','B7B','B8A','B8B')
# 
# for (i in 1:length(labels)) {
# 
#   T <- Read10X(data.dir = paste("~/Documents/Bladder/",labels[i],"_GEX/",sep=""))
#   T <- CreateSeuratObject(counts = T, project = "BladderNuclei", assay = "RNA", min.cells = 0, min.features = 0)
# 
#   mito.genes <- grep(pattern = "^MT-", x = rownames(x = T@assays$RNA@counts), value = TRUE)
#   percent.mito <- Matrix::colSums(T@assays$RNA@counts[mito.genes, ])/Matrix::colSums(T@assays$RNA@counts)
#   T <- AddMetaData(object = T, metadata = percent.mito, col.name = "percent.mito")
# 
#   p <- VlnPlot(T, features = "nCount_RNA", y.max = 50000)
#   ggsave(filename = paste("~/Documents/Bladder/",labels[i],"_GEX/UMI_VlnPlot.eps",sep=""), plot = p, device = 'eps', dpi = 300)
#   p <- VlnPlot(T, features = "percent.mito", y.max = 0.2)
#   ggsave(filename = paste("~/Documents/Bladder/",labels[i],"_GEX/mito_VlnPlot.eps",sep=""), plot = p, device = 'eps', dpi = 300)
# 
#   write.table(colnames(T@assays$RNA@counts), sep = "\n", file = paste("~/Documents/Bladder/",labels[i],"_GEX/filtered_cellbarcodes.txt",sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE)
#   rm(T)
# 
# }
# 
# rm(list=ls())

# add HTO counts to samples, remove doublets and negatives, and combine all samples together
labels <- c('B1','B2','B3','B4','B5A','B5B','B6A','B6B','B7A','B7B','B8A','B8B')

for (i in 1:length(labels)) {

  T <- Read10X(data.dir = paste("~/Documents/Bladder/seuratOutput/",labels[i],"_GEX/",sep=""))

  T <- CreateSeuratObject(counts = T, project = "BladderNuclei", assay = "RNA", min.cells = 0, min.features = 0)

  mito.genes <- grep(pattern = "^MT-", x = rownames(x = T@assays$RNA@counts), value = TRUE)
  percent.mito <- Matrix::colSums(T@assays$RNA@counts[mito.genes, ])/Matrix::colSums(T@assays$RNA@counts)
  T <- AddMetaData(object = T, metadata = percent.mito, col.name = "percent.mito")

  matrix_dir = paste("~/Documents/Bladder/citeseqcountOutput/",labels[i],"_HTO/umi_count/",sep="")
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat_hto <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat_hto) = barcode.names$V1
  rownames(mat_hto) = feature.names$V1

  GOI <- grep(pattern = 'Hash', x = rownames(mat_hto))
  mat_hto <- mat_hto[GOI,]
  BCs <- mat_hto@Dimnames[2]
  T <- T[,unlist(BCs)]
  T[["HASH"]] <- CreateAssayObject(mat_hto)
  T <- AddMetaData(T, metadata = labels[i], col.name = "sampleID")

  T <- NormalizeData(T, assay = "HASH", normalization.method = "CLR")

  if(labels[i] == "B1"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.97, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B2"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.96, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B3"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.9999, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B4"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.95, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B5A"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.9999, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B5B"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.9999, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B6A"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.98, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B6B"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.98, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B7A"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.99, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B7B"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.999, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B8A"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.999, kfunc = "kmeans", verbose = TRUE)
  }else if (labels[i] == "B8B"){
    T <- HTODemux(T, assay = "HASH", positive.quantile = 0.999, kfunc = "kmeans", verbose = TRUE)
  }

  print(labels[i])
  print(table(T$HASH_classification.global))
  write.csv(T@meta.data, file = paste("~/Documents/Bladder/seuratOutput/",labels[i],"_preDemux_metadata.csv",sep = ""))
  write.csv(as.matrix(T@assays$HASH@data), file = paste("~/Documents/Bladder/seuratOutput/",labels[i],"_preDemux_HASH_data.csv",sep = ""))
  T <- subset(T, idents = c("Negative","Doublet"), invert = TRUE)
  
  if (labels[i]=="B1" || labels[i]=="B2" || labels[i]=="B3" || labels[i]=="B4"){
    samples <- c('54','674','702','752','831','844','896','912')
    
    Idents(T) <- T$HASH_classification
    oldA <- c("Hash1-TTCCTGCCATTACTA","Hash2-CCGTACCTCATTGTT","Hash3-GGTAGATGTCCTCAG","Hash4-TGGTGTCATTCTTGA","Hash5-ATGATGAACAGCCAG","Hash6-CTCGAACGCTTATCG")
    newA <- c(paste("Bladder",samples[2*i-1],"_1",sep=""),paste("Bladder",samples[2*i-1],"_2",sep=""),paste("Bladder",samples[2*i-1],"_3",sep=""),paste("Bladder",samples[2*i],"_1",sep=""),paste("Bladder",samples[2*i],"_2",sep=""),paste("Bladder",samples[2*i],"_3",sep=""))
    Idents(T) <- plyr::mapvalues(x = Idents(T), from = oldA, to = newA)
    T[["Patient_Rep"]] <- Idents(T)
    
    Idents(T) <- T$Patient_Rep
    oldA <- c(paste("Bladder",samples[2*i-1],"_1",sep=""),paste("Bladder",samples[2*i-1],"_2",sep=""),paste("Bladder",samples[2*i-1],"_3",sep=""),paste("Bladder",samples[2*i],"_1",sep=""),paste("Bladder",samples[2*i],"_2",sep=""),paste("Bladder",samples[2*i],"_3",sep=""))
    newA <- c(paste("Bladder",samples[2*i-1],sep=""),paste("Bladder",samples[2*i-1],sep=""),paste("Bladder",samples[2*i-1],sep=""),paste("Bladder",samples[2*i],sep=""),paste("Bladder",samples[2*i],sep=""),paste("Bladder",samples[2*i],sep=""))
    Idents(T) <- plyr::mapvalues(x = Idents(T), from = oldA, to = newA)
    T[["Patient"]] <- Idents(T)
    
    T <- AddMetaData(T, metadata = "OLD", col.name = "cohort")
  }
  else if (labels[i]=="B5A" || labels[i]=="B5B"){
    Idents(T) <- T$HASH_classification
    oldA <- c("Hash1-TTCCTGCCATTACTA","Hash2-CCGTACCTCATTGTT","Hash3-GGTAGATGTCCTCAG","Hash4-TGGTGTCATTCTTGA","Hash5-ATGATGAACAGCCAG","Hash6-CTCGAACGCTTATCG")
    newA <- c("Bladder914","Bladder994","Bladder1245","Bladder1250","Bladder1192","Bladder1217")
    Idents(T) <- plyr::mapvalues(x = Idents(T), from = oldA, to = newA)
    T[["Patient"]] <- Idents(T)
    
    T <- AddMetaData(T, metadata = "NEW", col.name = "cohort")
  }
  else if (labels[i]=="B6A" || labels[i]=="B6B"){
    Idents(T) <- T$HASH_classification
    oldA <- c("Hash1-TTCCTGCCATTACTA","Hash2-CCGTACCTCATTGTT","Hash3-GGTAGATGTCCTCAG","Hash4-TGGTGTCATTCTTGA","Hash5-ATGATGAACAGCCAG","Hash6-CTCGAACGCTTATCG")
    newA <- c("Bladder489","Bladder59","Bladder590","Bladder824","Bladder36","Bladder72")
    Idents(T) <- plyr::mapvalues(x = Idents(T), from = oldA, to = newA)
    T[["Patient"]] <- Idents(T)
    
    T <- AddMetaData(T, metadata = "NEW", col.name = "cohort")
  }
  else if (labels[i]=="B7A" || labels[i]=="B7B"){
    Idents(T) <- T$HASH_classification
    oldA <- c("Hash1-TTCCTGCCATTACTA","Hash2-CCGTACCTCATTGTT","Hash3-GGTAGATGTCCTCAG","Hash4-TGGTGTCATTCTTGA","Hash5-ATGATGAACAGCCAG","Hash6-CTCGAACGCTTATCG")
    newA <- c("Bladder739","Bladder763","Bladder913","Bladder1126","Bladder1204","Bladder371")
    Idents(T) <- plyr::mapvalues(x = Idents(T), from = oldA, to = newA)
    T[["Patient"]] <- Idents(T)
    
    T <- AddMetaData(T, metadata = "NEW", col.name = "cohort")
  }
  else if (labels[i]=="B8A" || labels[i]=="B8B"){
    Idents(T) <- T$HASH_classification
    oldA <- c("Hash1-TTCCTGCCATTACTA","Hash2-CCGTACCTCATTGTT","Hash3-GGTAGATGTCCTCAG","Hash4-TGGTGTCATTCTTGA","Hash5-ATGATGAACAGCCAG","Hash6-CTCGAACGCTTATCG")
    newA <- c("Bladder593","Bladder419","Bladder446","Bladder485","Bladder518","Bladder8")
    Idents(T) <- plyr::mapvalues(x = Idents(T), from = oldA, to = newA)
    T[["Patient"]] <- Idents(T)
    
    T <- AddMetaData(T, metadata = "NEW", col.name = "cohort")
  }
  
  if (i == 1) {
    combinedSobj <- RenameCells(object = T, add.cell.id = paste(labels[i],"_",sep=""))
  }
  else{
    T <- RenameCells(object = T, add.cell.id = paste(labels[i],"_",sep=""))
    combinedSobj <- merge(combinedSobj, T)
  }

  rm(T,GOI)

}

# deal with first sample
T <- Read10X(data.dir = paste("~/Documents/Bladder/seuratOutput/B0_GEX/",sep=""))
T <- CreateSeuratObject(counts = T, project = "BladderNuclei", assay = "RNA", min.cells = 0, min.features = 0)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = T@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(T@assays$RNA@counts[mito.genes, ])/Matrix::colSums(T@assays$RNA@counts)
T <- AddMetaData(object = T, metadata = percent.mito, col.name = "percent.mito")

matrix_dir = paste("~/Documents/Bladder/citeseqcountOutput/B0_HTO/umi_count/",sep="")
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat_hto <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat_hto) = barcode.names$V1
rownames(mat_hto) = feature.names$V1

GOI <- grep(pattern = 'Hash', x = rownames(mat_hto))
mat_hto <- mat_hto[GOI,]
BCs <- mat_hto@Dimnames[2]
T <- T[,unlist(BCs)]
T[["HASH"]] <- CreateAssayObject(mat_hto)
T <- AddMetaData(T, metadata = "B0", col.name = "sampleID")

T <- NormalizeData(T, assay = "HASH", normalization.method = "CLR")
T <- HTODemux(T, assay = "HASH", positive.quantile = 0.9, kfunc = "kmeans", verbose = TRUE)
write.csv(T@meta.data, file = paste("~/Documents/Bladder/seuratOutput/B0_preDemux_metadata.csv",sep = ""))
write.csv(as.matrix(T@assays$HASH@data), file = paste("~/Documents/Bladder/seuratOutput/B0_preDemux_HASH_data.csv",sep = ""))
print(table(T$HASH_classification.global))
T <- subset(T, idents = c("Negative","Doublet"), invert = TRUE)

T <- AddMetaData(T, metadata = "Bladder1246", col.name = "Patient")
T <- AddMetaData(T, metadata = "B0", col.name = "sampleID")

Idents(T) <- T$HASH_classification
oldA <- c("Hash1-TTCCTGCCATTACTA","Hash2-CCGTACCTCATTGTT","Hash3-GGTAGATGTCCTCAG")
newA <- c("Bladder1246_1","Bladder1246_2","Bladder1246_3")
Idents(T) <- plyr::mapvalues(x = Idents(T), from = oldA, to = newA)
T[["Patient_Rep"]] <- Idents(T)

T <- AddMetaData(T, metadata = "OLDEST", col.name = "cohort")

T <- RenameCells(object = T, add.cell.id = "B0_")

combinedSobj <- merge(combinedSobj, T)

rm(T,GOI)

saveRDS(combinedSobj, file = "~/Documents/Bladder/seuratOutput/combinedSobj_filtered.rds")
write.csv(combinedSobj@meta.data, file = "~/Documents/Bladder/seuratOutput/combinedSobj_filtered_metadata.csv")