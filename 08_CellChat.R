# CellChat Comparison

## Load the libraries

library(CellChat)
library(patchwork)

## Create a directory to save figures

data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)

getwd()

## Load CellChat object of each dataset and then merge together

cellchatAD <- readRDS('/home/jovyan/AsymDENV/cellchat_AD.rds')
cellchatDHF <- readRDS('/home/jovyan/AsymDENV/cellchat_DHF.rds')

object.list <- list(AD = cellchatAD, DHF = cellchatDHF)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

## Part I: Predict general principles of cell-cell communication

### Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

### Differential number of interactions or interaction strength among different cell populations

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

### Part II: Identify the upgulated and down-regulated signaling ligand-receptor pairs

netVisual_bubble(cellchat, sources.use = 20, targets.use = c(1:13),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = 20, targets.use = c(13:13),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in DHF", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 20, targets.use = c(13:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in DHF", angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(13:13),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in DHF", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(13:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in DHF", angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(13:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in DHF: B memory", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(15), targets.use = c(13:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in DHF: CD16 Mono", angle.x = 45, remove.isolate = T)
gg1 + gg2

ggsave('BmemCD16_CD8TEM.pdf')

gg1 <- netVisual_bubble(cellchat, sources.use = 20, targets.use = c(13:13),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in DHF", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 20, targets.use = c(13:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in DHF", angle.x = 45, remove.isolate = T)
gg1 + gg2

### Identify dysfunctional signaling by using differential expression analysis

pos.dataset = "DHF"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "DHF",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "AD",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 3, targets.use = c(13:13), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 3, targets.use = c(13:13), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

## Part III: Visually compare cell-cell communication for IL-10

pathways.show <- c("IL10") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
