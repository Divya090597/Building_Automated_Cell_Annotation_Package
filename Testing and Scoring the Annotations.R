

# Load Required Libraries -------------------------------------------------

library(SingleR)
library(data.table)
library(celldex)
library(ggparty)
library(igraph)
library(tidyverse)
library(data.tree)

# SingleR -----------------------------------------------------------------

# Load in seurat object

seu_obj = readRDS("~/data/Kidney/ccRCC_Updated_Seurat_Object_Cellcycle_score_SI_18854_Tumor.rds")

#  Load Reference Data (HumanPrimaryCellAtlas)

ref <- celldex::HumanPrimaryCellAtlasData()

View(as.data.frame(colData(ref)))

# Extract the count data from Seurat object to use for SingleR predictions
 
pred = GetAssayData(seu_obj, slot = "counts")

# Run SingleR for cell type prediction using the reference data

pred_data = SingleR(test = pred, ref = ref, labels = ref$label.main)
pred_data
View(as.data.frame(pred_data))

# Save the SingleR prediction results as a CSV file

fwrite(as.data.frame(pred_data), row.names = T, "~/data/Kidney/SingleR/Output.csv" )


# Visualize the Results 

plotScoreHeatmap(pred_data)
plotDeltaDistribution(pred_data)



# ScGate Annotation -------------------------------------------------------



# Visualizing gating models

# Define scGate Models for Cell Types

scGate::plot_tree(scGate_models_DB$human$generic$Myeloid)

# Define scGate Models for Cell Types

##...............  1. B Cell Model

my_scGate_model <- gating_model(name = "Bcell", signature = c("MS4A1"))

# Test the model performance on B cells

panBcell.performance <- test_my_model(my_scGate_model, target = "Bcell")

# Apply the scGate model to the Seurat object

seu_obj <- scGate(data = seu_obj, model = my_scGate_model, verbose = T)

# Define color palette and plot the results

palette <- c(list(Impure = "gray", Pure = "green"))
DimPlot(seu_obj_1, cols = palette) + theme(aspect.ratio = 1) + ggtitle("Bcell")


##............... 2.CD8 T Cell Model

my_scGate_model <- gating_model(name = "CD8T", signature = c("CD8A"))

obj <- scGate(data = seu_obj, model = my_scGate_model)

DimPlot(obj, cols = palette ) + theme(aspect.ratio = 1) + ggtitle("CD8 Tcell")

my_scGate_model <- gating_model(name = "immune", signature = c("PTPRC"), level = 1)  # initialize model with one positive signature

my_scGate_model <- gating_model(model = my_scGate_model, name = "CD8T_Cell", signature = c("CD8A"), level = 2) # add positive signature at second step
                                                                                                   

# Apply the scGate model to the Seurat object

seu_obj <- scGate(data = seu_obj, model = my_scGate_model, save.levels = T)

# Plot the hierarchical gating tree for Macrophages

wrap_plots(plot_levels(seu_obj))+plot_annotation(title = "Hierarchichal Gating of TAM")


# Visualize the hierarchical gating results

DimPlot(seu_obj, cols = palette) + theme(aspect.ratio = 1)+ ggtitle("CD8T")


# Performance Metrics

scGate::performance.metrics(obj$scGate_multi %in% c("CD8T"), obj$is.pure == "Pure")
FeaturePlot(obj, features = c("CD8T")) + scale_color_viridis(option = "D") +
  theme(aspect.ratio = 1)



##.................3.Macrophage Model with Hierarchical Gating


# Define the hierarchical gating model for Macrophages

my_scGate_model <- gating_model(name = "immune", signature = c("PTPRC"), level = 1)  # initialize model with one positive signature

my_scGate_model <- gating_model(model = my_scGate_model, name = "Macrophage", signature = c("CD68",
                                                                                            "FCGR1A"), level = 2)  # add positive signature at second step
my_scGate_model <- gating_model(model = my_scGate_model, name = "LA-TAMS", signature = c("APOC1",
                                                                                      "APOE"), level = 3)


# Apply the scGate model to the Seurat object

seu_obj <- scGate(data = seu_obj, model = my_scGate_model, save.levels = T)

# Visualize the hierarchical gating results

DimPlot(seu_obj, cols = palette) + theme(aspect.ratio = 1)+ ggtitle("TAM")

# Plot the hierarchical gating tree for Macrophages

plots <- scGate::plot_levels(seu_obj)

wrap_plots(plot_levels(seu_obj))+plot_annotation(title = "Hierarchichal Gating of TAM")



# ScType Annotation -------------------------------------------------------


# prepare edges

cL_resutls <- cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes

nodes_lvl1 <- sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss <- c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes <- rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db <- openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph

gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

DimPlot(seu_obj, reduction = "umap",group.by = "RNA_snn_res.0.5", repel = TRUE, cols = ccolss)+ gggr


# Check for duplicate vertex names in the nodes dataframe

duplicates <- nodes[duplicated(nodes$cluster), ]
print(duplicates)

# If you want to append a unique number to the duplicate names

nodes$cluster <- make.unique(as.character(nodes$cluster))

# Create the graph again after handling duplicates

mygraph <- graph_from_data_frame(edges, vertices = nodes)

# load auto-detection function
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")


# guess a tissue type
#tissue_guess <- auto_detect_tissue_type(path_to_db_file = db_, seuratObject = seu_obj, scaled = TRUE, assay = "RNA")  # if scaled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used         






