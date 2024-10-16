

# Load in Required Libraries ----------------------------------------------

library(Seurat)          
library(ggplot2)         
library(dplyr)           
library(data.table)      
library(future)          
library(scGate)          
library(patchwork)


# Import the seurat object

seurat_obj = readRDS("~/data/Kidney/ccRCC_seurat_object_GSE159115_SI_18854_Tumor.rds")


fix(my_scGate_model)  # edit locally in R
write.table(my_scGate_model, "Scgate_T_cell_example_genelist_2024_10_04.tsv", sep = "\t")  # export and then edit on Excel

# Load my scGate Model

my_scGate_model_gene_list1 <- scGate::load_scGate_model("Scgate_T_cell_example_genelist_2024_10_04.tsv")
my_scGate_model_gene_list2 <- scGate::load_scGate_model("Scgate_Myeloidcell_example_genelist_2024_10_04.tsv.txt")
# Apply scGate Model to Seurat Object

# Only for single model Testimng
obj <- scGate(seurat_obj, model = my_scGate_model_gene_list2, save.levels = TRUE)

# Visualize Pure and Impure Cells with UMAP Plot

palette <- c(list(Impure = "gray", Pure = "green"))
DimPlot(obj, cols = palette) + theme(aspect.ratio = 1)

# Visualize scGate Levels
plots <- scGate::plot_levels(obj)
wrap_plots(plots, ncol = 2)

# Remove non-ASCII characters (including <a0> or others)
my_scGate_model_gene_list1$signature <- iconv(my_scGate_model_gene_list1$signature, "UTF-8", "ASCII", sub = "") 
my_scGate_model_gene_list1$signature <- gsub(" ", "", my_scGate_model_gene_list1$signature)  # Remove all spaces
my_scGate_model_gene_list1$signature <- gsub(",", ";", my_scGate_model_gene_list1$signature)  # Replace commas with semicolons
my_scGate_model_gene_list1$signature <- gsub("hi", "", my_scGate_model_gene_list1$signature)  # Remove "hi"


##### scGate Annoatation #####
# List of model names and corresponding scGate models

# List of models
models <- list(

  Divya_model_Myeloid = my_scGate_model_gene_list2
)

# Create a list to store new columns
new_columns_list <- list()

# Loop through each model
for (model_name in names(models)) {
  # Get the current model
  my_scGate_model <- models[[model_name]]
  
  # Copy the original Seurat object to avoid modification
  seu_obj_copy <- seurat_obj
  
  # Apply scGate to the copied Seurat object
  seu_obj_processed <- scGate(seu_obj_copy, model = my_scGate_model, pos.thr = 0.2, neg.thr = 0.2, ncores = 4)
  
  # Assign the processed Seurat object to a variable named according to the model name
  assign(paste0("seu_obj_", model_name), seu_obj_processed)
  
  # Extract metadata
  meta_data <- seu_obj_processed@meta.data
  
  ##### Process "is.pure" Column #####
  # Initialize the new column to store the cell types based on "is.pure"
  new_col_name <- paste0("scGate_", model_name, "_Celltype")
  meta_data[[new_col_name]] <- NA
  
  # Get barcodes of cells marked as "Pure"
  pure_barcodes <- rownames(meta_data)[meta_data$is.pure == "Pure"]
  
  # Annotate these cells with the respective feature
  for (barcode in pure_barcodes) {
    # Get the relevant UCell scores (using columns that end with "_UCell")
    cell_type_columns <- grep("_UCell$", colnames(meta_data), value = TRUE)
    
    # Identify the cell type with the highest UCell score
    highest_score <- -Inf
    highest_cell_type <- NA
    
    for (col in cell_type_columns) {
      score <- meta_data[barcode, col]
      if (!is.na(score) && score > highest_score) {
        highest_score <- score
        highest_cell_type <- sub("_UCell$", "", col)  # Remove "_UCell" suffix to get the cell type
      }
    }
    
    # Assign the highest cell type to the new column
    if (!is.na(highest_cell_type)) {
      meta_data[barcode, new_col_name] <- highest_cell_type
    }
  }
  
  # Store the new "Celltype" column in the list
  new_columns_list[[new_col_name]] <- meta_data[, new_col_name, drop = F]
  
  ##### Process "scGate_multi" Column #####
  if ("scGate_multi" %in% colnames(meta_data)) {
    # Rename and update the column only if it exists
    multi_col_name <- paste0("scGate_multi_", model_name)
    
    # Initialize the scGate_multi column if it is NA or "Target"
    meta_data$scGate_multi[is.na(meta_data$scGate_multi) | meta_data$scGate_multi == "Target"] <- NA
    
    # Update scGate_multi column with the new cell type (from the highest scoring UCell)
    for (barcode in pure_barcodes) {
      cell_type <- meta_data[barcode, new_col_name]
      if (!is.na(cell_type)) {
        if (is.na(meta_data[barcode, "scGate_multi"]) || meta_data[barcode, "scGate_multi"] == "Target") {
          meta_data[barcode, "scGate_multi"] <- cell_type
        } else {
          # Concatenate multiple features if already assigned
          meta_data[barcode, "scGate_multi"] <- paste(meta_data[barcode, "scGate_multi"], cell_type, sep = " ")
        }
      }
    }
    
    # Only store the column if it exists and has been modified
    if (multi_col_name %in% colnames(meta_data)) {
      new_columns_list[[multi_col_name]] <- meta_data[, multi_col_name, drop = F]
    }
  } else {
    cat("scGate_multi column not found in metadata for model:", model_name, "\n")
  }
  
  # Save the updated metadata to a CSV file
  write.csv(meta_data, file = paste0(dir, "/meta_data_", model_name, ".csv"))
  
  # View the metadata in the RStudio viewer or console
  View(meta_data)
}



new_columns_df <- do.call(cbind, new_columns_list)

# Convert the combined list into a data frame
new_columns_df <- as.data.frame(new_columns_df)

# View the final data frame
View(new_columns_df)

new_columns_df$cellnames <- rownames(new_columns_df)
seu_obj_Divya_model_Myeloid@meta.data$cellnames = rownames(seu_obj_Divya_model_Myeloid@meta.data)

New_seu_obj = merge(seu_obj_Divya_model_Myeloid@meta.data, new_columns_df, by = "cellnames", rownames = T)


# Set the 'cellnames' as rownames
rownames(New_seu_obj) <- New_seu_obj$cellnames

# Assign the updated metadata back to the Seurat object
seu_obj_processed@meta.data = New_seu_obj

# View the updated Seurat object metadata
View(seu_obj_processed@meta.data)


DimPlot(seu_obj_processed, reduction = "umap", group.by = "scGate_Divya_model_Myeloid_Celltype")+ 
  ggtitle("Divya_Model_Cell_Type_ccRCC")

# Save the entire environment to a file called "my_environment.RData"
save.image(file = "my_environment_Divya_Sc_gate.RData")

#Save updated seurat object

saveRDS(seu_obj_processed, "~/data/Kidney/ScGate_project/Updated_Seurat_object_Divya_Cell_type_column_10_06_2024.rds")


# List of marker genes
marker_genes= c("CCL2","CCL3","CCL4","CCL7","CCL8", "CCL3L3", "CCL4L2","CD274", "CD40","CXCL1", "CXCL2","CXCL3","CXCL5","CXCL8","CXCL9","CXCL10","CXCL11", "IDO1", "IFI6", "IFIT1","IFIT2","IFIT3", "IFITM1","IFITM3", "IRF1", "IRF7", "ISG15", "LAMP3", "PDCD1LG2", "TNFSF10", "C1QA","C1QC", "CD38", "IL4I1", "IFI44L","CCL5","CCL20","G0S2", "IL1B","IL1RN","INHBA", "KLF2", "KLF6", "NEDD9", "PMAIP1", "S100A8","S100A9","SPP1","ACP5", "APOC1", "ATF1", "C1QB", "CCL18", "CD163", "CD36", "CD63", "CHI3L1", "CTSB","CTSD","CTSL", "F13A1", "FABP5", "FOLR2", "GPNMB", "IRF3", "LGALS3","NRF1", "NUPR1", "TFDP2", "TREM2")
marker_genes_2 = c("APOC1","APOE","ACP5","FABP5","VEGFA","SPP1","IL1B","CCL3","CXCL1","CXCL2","CXCL3","CXCL5","IDO1","ISG15","CXCL8","CXCL9","CXCL10","MRC1","CD274","CX3CR1","MKI67","CDK1","LYVE1","HES1","FOLR2","CD80","CD40","CD86","HLA-DQA1","HLA-DQB1","CHIT1","ITGA4","TGFB2","HLA-DRA","HLA-DRB1","HLA-DRB5","ICOSLG","LGALS9")

# Check which of the marker genes are present in the dataset
present_genes <- marker_genes %in% rownames(seu_obj_processed)
missing_genes <- marker_genes[!present_genes]  # Identify missing genes
data.frame(marker_genes, present_genes)  # Show which genes are present

# Replace NA values in the scGate_Divya_model_Celltype column with "Unknown"
seu_obj_processed@meta.data$scGate_Divya_model_Celltype[is.na(seu_obj_processed@meta.data$scGate_Divya_model_Celltype)] <- "Unknown"

# Generate the dot plot with the updated column
dotplot <- DotPlot(seu_obj_processed, features = markers, group.by = "scGate_Divya_model_Celltype") + 
  RotatedAxis() + 
  scale_color_gradientn(colors = c("lightblue", "blue", "red")) +
  theme(axis.text.x = element_text(size = 10))

# Display the dotplot
print(dotplot)


markers = c("CD14","S100A4","S100A6","CD36","CEBPD","NR4A1","NLRP3","NR4A2","AREG","EGR1","TCF7L2","IFITM2","SERPINA1","SIGLEC10","RHOC","MS4A7","ZFP36L2","CXCL10","LAMP3", "PDCD1LG2","SLC2A1","ISG15","IDO1","IRF1","CD86","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DQA1","HLA-DQB1","MRC1","CX3CR1","TREM2","CXCL1", "CXCL2","CXCL3","CXCL8","CCL3","NFKB1","NFE2L2","REL","VEGFA","SPP1","VCAN","FCN1","THBS1","APOE","APOC1","ACP5","FABP5","LYVE1","HES1","FOLR2","MARCO","CCNA2", "CDC45", "CDK1","MKI67","C1QA","C1QC","C1QB","CD274")

# Create your base dot plot
dotplot <- DotPlot(seu_obj_processed, features = markers, group.by = "scGate_Divya_model_Celltype") + 
  RotatedAxis() + 
  scale_color_gradientn(colors = c("lightblue", "blue", "red")) +
  theme(axis.text.x = element_text(size = 10))

# Define the x-axis ranges for different sets of markers (column positions in your plot)
# For example, if markers from columns 1-10 are for "Cell Type A"
rect_data <- data.frame(
  xmin = c(0.5, 10.5, 20.5),  # Define starting x positions for rectangles
  xmax = c(10.5, 20.5, 30.5),  # Define ending x positions for rectangles
  ymin = -Inf,                 # ymin to cover the full y-range
  ymax = Inf                   # ymax to cover the full y-range
)

# Add the rectangles to the plot using geom_rect
dotplot + 
  geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", size = 1)  # Customize the border color and size

# Create the base dot plot first
dotplot <- DotPlot(seu_obj_processed, features = markers, group.by = "scGate_Divya_model_Celltype") + 
  RotatedAxis() + 
  scale_color_gradientn(colors = c("lightblue", "blue", "red")) +
  theme(axis.text.x = element_text(size = 10))

# Add rectangles manually using annotation_custom
dotplot + 
  annotation_custom(grob = rectGrob(gp = gpar(col = "black", fill = NA)), xmin = 0.5, xmax = 10.5, ymin = -Inf, ymax = Inf) +
  annotation_custom(grob = rectGrob(gp = gpar(col = "black", fill = NA)), xmin = 10.5, xmax = 20.5, ymin = -Inf, ymax = Inf) +
  annotation_custom(grob = rectGrob(gp = gpar(col = "black", fill = NA)), xmin = 20.5, xmax = 30.5, ymin = -Inf, ymax = Inf)
