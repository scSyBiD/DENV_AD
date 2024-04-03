library(Seurat)
library(SoupX)
library(DropletUtils)

#Assemble file path names
path_fltr <- paste('GEX/filtered/', list.files(path = 'GEX/filtered/'), sep = '')
path_raw <- paste('GEX/raw/', list.files(path = 'GEX/raw/'), sep = '')

#Load the matrices (both filtered and raw to run SoupX)
mtx_fltr <- NULL
mtx_raw <- NULL
for (i in 1:length(path_fltr)) {
  mtx_fltr[[i]] <- Read10X(path_fltr[[i]])
  mtx_raw[[i]] <- Read10X(path_raw[[i]])
}

sample_id <- gsub('_filtered_feature_bc_matrix', '', list.files(path = 'GEX/filtered/'))

for (i in 1:length(path_fltr)) {
  mtx_fltr[[i]]@Dimnames[[2]] <- paste(sample_id[i], mtx_fltr[[i]]@Dimnames[[2]], sep = '_')
  mtx_fltr[[i]]@Dimnames[[2]] <- gsub('-1$', '', mtx_fltr[[i]]@Dimnames[[2]])
  
  mtx_raw[[i]]@Dimnames[[2]] <- paste(sample_id[i], mtx_raw[[i]]@Dimnames[[2]], sep = '_')
  mtx_raw[[i]]@Dimnames[[2]] <- gsub('-1$', '', mtx_raw[[i]]@Dimnames[[2]])
}

srat <- NULL
soup.channel <- NULL
meta <- NULL
umap <- NULL
adj.matrix <- NULL

for (i in 1:length(path_fltr)) {
  srat[[i]] <- CreateSeuratObject(mtx_fltr[[i]])
  srat[[i]] <- SCTransform(srat[[i]], verbose = F)
  srat[[i]] <- RunPCA(srat[[i]], verbose = F)
  srat[[i]] <- RunUMAP(srat[[i]], dims = 1:30, verbose = F)
  srat[[i]] <- FindNeighbors(srat[[i]], dims = 1:30, verbose = F)
  srat[[i]] <- FindClusters(srat[[i]], verbose = F)

#Set up soup channels
  soup.channel[[i]] <- SoupChannel(mtx_raw[[i]], mtx_fltr[[i]])
  meta[[i]] <- srat[[i]]@meta.data
  umap[[i]] <- srat[[i]]@reductions$umap@cell.embeddings
  soup.channel[[i]] <- setClusters(soup.channel[[i]], setNames(meta[[i]]$seurat_clusters, rownames(meta[[i]])))
  soup.channel[[i]] <- setDR(soup.channel[[i]], umap[[i]])
  soup.channel[[i]] <- autoEstCont(soup.channel[[i]])
  
  head(soup.channel[[i]]$soupProfile[order(soup.channel[[i]]$soupProfile$est, decreasing = T), ], n = 5)
  
  adj.matrix[[i]] <- adjustCounts(soup.channel[[i]], roundToInt = T)
  DropletUtils::write10xCounts(path = paste('SoupCorrectedGEX/', sample_id[i], sep = ''), adj.matrix[[i]])
}

for (i in 3:9) {
  srat[[i]] <- CreateSeuratObject(mtx_fltr[[i]])
  srat[[i]] <- SCTransform(srat[[i]], verbose = F)
  srat[[i]] <- RunPCA(srat[[i]], verbose = F)
  srat[[i]] <- RunUMAP(srat[[i]], dims = 1:30, verbose = F)
  srat[[i]] <- FindNeighbors(srat[[i]], dims = 1:30, verbose = F)
  srat[[i]] <- FindClusters(srat[[i]], verbose = F)
  
  #Set up soup channels
  soup.channel[[i]] <- SoupChannel(mtx_raw[[i]], mtx_fltr[[i]])
  meta[[i]] <- srat[[i]]@meta.data
  umap[[i]] <- srat[[i]]@reductions$umap@cell.embeddings
  soup.channel[[i]] <- setClusters(soup.channel[[i]], setNames(meta[[i]]$seurat_clusters, rownames(meta[[i]])))
  soup.channel[[i]] <- setDR(soup.channel[[i]], umap[[i]])
  soup.channel[[i]] <- autoEstCont(soup.channel[[i]])
  
  head(soup.channel[[i]]$soupProfile[order(soup.channel[[i]]$soupProfile$est, decreasing = T), ], n = 5)
  
  adj.matrix[[i]] <- adjustCounts(soup.channel[[i]], roundToInt = T)
  DropletUtils::write10xCounts(path = paste('SoupCorrectedGEX/', sample_id[i], sep = ''), adj.matrix[[i]])
}

