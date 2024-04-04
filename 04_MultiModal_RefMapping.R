# Run Serial Multimodal Reference Mapping

# Install packages if needed
# BiocManager::install('multtest')
# install.packages('Seurat')
# remotes::install_github("mojaveazure/seurat-disk")
# install.packages("Matrix", repos="http://R-Forge.R-project.org")

# Load required libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sceasy)
library(data.table)

# Prep the reference
reference <- LoadH5Seurat('pbmc_multimodal.h5seurat', verbose = FALSE)

# Check annotations in the reference dataset
# DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l3", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# Prep the query dataset
# Convert h5ad to seurat
# sceasy::convertFormat('AsymDENV_QC.h5ad', from="anndata", to="seurat", outFile='AsymDENV_QC_srat.rds')

DENV <- readRDS('AsymDENV_QC_srat.rds')

DENV.list <- SplitObject(DENV, split.by = 'batch')

DENV.list <- lapply(X = DENV.list, FUN = SCTransform)

# Mapping and Annotation
for (i in 1:25) {
  anchors <- FindTransferAnchors(
    reference = reference,
    query = DENV.list[[i]],
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )

  DENV.list[[i]] <- MapQuery(
    anchorset = anchors,
    query = DENV.list[[i]],
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      celltype.l3 = "celltype.l3",
      predicted_ADT = "ADT" #predicted protein levels based on CITE-seq
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
}

# Save the annotated list
saveRDS(DENV.list, 'AsymDENV_SCT_list.rds')

# Uncomment if you need to load the list again
# DENV.list = readRDS(DENV.list, "AsymDENV_SCT_list.rds")

# Aggregate metadata for cell type predictions
Celltype <- NULL              
for (i in 1:25) {
  Celltype <- rbind(Celltype, DENV.list[[i]]@meta.data %>% 
                    select(predicted.celltype.l1.score, predicted.celltype.l1, predicted.celltype.l2.score, predicted.celltype.l2, predicted.celltype.l3.score, predicted.celltype.l3))
}

# Combine embeddings for all batches
refspca <- NULL
refumap <- NULL
for (i in 1:25) {
  refspca <- rbind(refspca, as.data.frame(Embeddings(DENV.list[[i]], reduction = 'ref.spca')))
  refumap <- rbind(refumap, as.data.frame(Embeddings(DENV.list[[i]], reduction = 'ref.umap')))
}

# Save the results to CSV files
fwrite(Celltype, 'MMRP_celltype.csv', row.names = TRUE)
fwrite(as.data.frame(refspca), 'MMRP_spca.csv', row.names = TRUE)
fwrite(as.data.frame(refumap), 'MMRP_umap.csv', row.names = TRUE)
