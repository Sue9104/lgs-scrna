library(CellChat)
library(patchwork)
library(glue)

indir <- "~/projects/seq3/20230530/4_scrna/counts/8_communication"
pro <- readRDS(glue("{indir}/RIF.isoforms.pro.interaction.cellchat.system.rds"))
sec <- readRDS(glue("{indir}/RIF.isoforms.sec.interaction.cellchat.system.rds"))
object.list <- list(pro=pro, sec=sec)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat




# Part I: Predict general principles of cell-cell communication
## total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

## Compare the number of interactions among different cell populations
#### red: increased; blue: decreased
#### chord
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#### heatmap
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
## interactions or interaction strength between any two cell types in each dataset
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
## differential interaction among cell types
group.cellType <- c(rep("LYM", 2), rep("MYE", 2), c("ENDO"), 
                    rep("EPI", 4), rep("STR", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "EPI", "ENDO", "LYM", "MYE"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T)


## differential number of interactions or interaction strength
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T)
## Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
## specific signaling changes of Inflam.DC and cDC1 between NL and LS
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "eStr1", signaling.exclude = "MK")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "eStr2", signaling.exclude = "MK")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "eStr3", signaling.exclude = "MK")
gg1
gg2
gg3
data1 <-  gg1$data %>% filter(specificity != "Shared") %>% mutate(celltype = "eStr1")
data2 <-  gg2$data %>% filter(specificity != "Shared") %>% mutate(celltype = "eStr2")
data3 <-  gg3$data %>% filter(specificity != "Shared") %>% mutate(celltype = "eStr3")
data_outgoing <- bind_rows(data1, data2, data3) %>% 
  filter(specificity.out.in == "Outgoing specific") %>% 
  pivot_wider(id_cols = c(labels, specificity), 
              names_from = celltype, values_from = outgoing) %>% 
  arrange(specificity)
data_outgoing
data_incoming <- bind_rows(data1, data2, data3) %>% 
  filter(specificity.out.in == "Incoming specific") %>% 
  pivot_wider(id_cols = c(labels, specificity), 
              names_from = celltype, values_from = incoming) %>% 
  arrange(specificity) 
data_incoming

#patchwork::wrap_plots(plots = list(gg1,gg2, gg3))


# Part II: Identify the conserved and context-specific signaling pathways
## Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5, 
                            title = "functional similarity")
#netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

## Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5,
                            title = "structure similarity")
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

## Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")

## Identify and visualize the conserved and context-specific signaling pathways
### Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
### Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
#### outgoing
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(
  object.list[[i]], pattern = "outgoing", signaling = pathway.union, 
  title = names(object.list)[i], width = 5, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(
  object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, 
  title = names(object.list)[i+1], width = 5, height = 16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#### incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#### overall
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 16, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 16, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



# Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#### cEpi 6, gEpi 7, lEpi 9, eStr1 13, eStr2 12, eStr3 11 
#### T 1, NK 2, pDC 3, Mac 4, Endo 5, pro 8, Endo 10
netVisual_bubble(cellchat, sources.use = c(7,9), targets.use = c(13,12,11),  comparison = c(1, 2), angle.x = 45)
## upgulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs
gg1 <- netVisual_bubble(cellchat, sources.use = c(7,9), targets.use = c(13,12,11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(7,9), targets.use = c(13,12,11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
gg1 + gg2
## Identify dysfunctional signaling by using differential expression analysis
pos.dataset = "pro"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "pro",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "sec",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2
#### chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], targets.use = c(13,12,11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 13, targets.use = c(12,11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
computeEnrichmentScore(net.up, species = 'human')
computeEnrichmentScore(net.down, species = 'human')

# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
#### heatmap
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
#### Chord diagram at celltype level
pathways.show <- c("PRL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
#### Chord diagram at gene level
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(5:8), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}

## compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2, 3, 4), targets.use = c(8,10),  
                       title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
}
## show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}

# Part V: Compare the signaling gene expression distribution between different datasets
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("pro", "sec")) # set factor level
plotGeneExpression(cellchat, signaling = "PARs", split.by = "datasets", colors.ggplot = T)

# Interested
obj.pro <- object.list[[1]]
obj.sec <- object.list[[2]]
a <- levels(obj.pro@idents)
names(a) <- sequence(length(a))
a
obj.pro@netP$pathways

pdf(glue::glue("{indir}/RIF.stage_comm.pdf"))
vertex.receiver = c(11,12,13) # a numeric vector. 
for (pathways.show in intersect(obj.pro@netP$pathways, obj.sec@netP$pathways)){
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_aggregate(obj.pro, signaling = pathways.show, 
                      vertex.receiver = vertex.receiver, layout = "hierarchy")
  netVisual_aggregate(obj.sec, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "hierarchy")  
}

ggplot() + ggtitle("Only pro")
for (pathways.show in setdiff(obj.pro@netP$pathways, obj.sec@netP$pathways)){
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_aggregate(obj.pro, signaling = pathways.show, 
                      vertex.receiver = vertex.receiver, layout = "hierarchy")
}

ggplot() + ggtitle("Only sec")
for (pathways.show in setdiff(obj.sec@netP$pathways, obj.pro@netP$pathways)){
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_aggregate(obj.sec, signaling = pathways.show, vertex.receiver = vertex.receiver, layout = "hierarchy")  
}

for (targets in seq(10)){
  p1 <- netVisual_bubble(obj.pro, sources.use = c(11,12,13), targets.use = targets, remove.isolate = FALSE)
  p2 <- netVisual_bubble(obj.sec, sources.use = c(11,12,13), targets.use = targets, remove.isolate = FALSE)
  p1 / p2
}
dev.off()
