library(data.table)
library(pheatmap)

dt <- fread(snakemake@input[[1]], sep = "\t")

# 2. Filter for colocalising hits
h4_threshold <- 0.8
dt_coloc <- dt[PP.H4.abf >= h4_threshold]

# 3. Create a wide format matrix for the heatmap
# Formula: Row ~ Column, value.var: the cell contents
# We use fill = 0 for cases where a SNP-Trait pair wasn't tested or didn't meet criteria
dt_wide <- dcast(
  dt_coloc, 
  `Isotype's lead SNP` ~ `Non-Ig trait`, 
  value.var = "Pearson correlation", 
  fill = 0
)

# Convert to matrix and set row names
mat <- as.matrix(dt_wide[, -1])
rownames(mat) <- dt_wide[[1]]

# 4. Prepare Row Annotation (to group by Isotype)
# Get unique mappings of SNP to Isotype
row_info <- unique(dt_coloc[, .(`Isotype's lead SNP`, Isotype)])

# Set factor levels to ensure specific ordering: IgA, IgG, IgM
row_info[, Isotype := factor(Isotype, levels = c("IgA", "IgG", "IgM"))]

# Sort the SNPs based on the Isotype factor
setorder(row_info, Isotype)
snp_order <- row_info[[1]]

# Reorder the matrix and create the annotation dataframe for pheatmap
mat <- mat[snp_order, , drop = FALSE]
row_ann <- data.frame(Isotype = row_info$Isotype)
rownames(row_ann) <- row_info[[1]]

# 5. Define Heatmap Colors
# Using a blue-white-red palette centered at 0 to show correlation sign
max_val <- max(abs(mat), na.rm = TRUE)
palette_colors <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-1, 1, length.out = 101)
legend_breaks <- seq(-1, 1, by = 0.5)
color_map_breaks <- seq(-1, 1, length.out = 101)

# 6. Plot Heatmap
png(snakemake@output[[1]], width = 7, height = 10, units = "in", res = 300)
pheatmap(
  mat,
  annotation_row = row_ann,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = palette_colors,
  breaks = color_map_breaks,
  legend_breaks = legend_breaks,
  legend_labels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
  
  # --- FIXING THE CLIPPING ---
  # Increase the left margin by adding padding to the row names area
  # or using the 'silent' return to manually adjust the viewport.
  # However, the most direct way in pheatmap is:
  
  treeheight_row = 0,      # Remove row dendrogram space since we aren't clustering rows
  annotation_names_row = TRUE, # Ensure "Isotype" label appears above the color bar
  
  # --- SPACE ALLOCATION ---
  # If row names are long, pheatmap might clip the annotation.
  # We can decrease the font size further or use margins:
  fontsize_row = 9, 
  
  # This parameter handles the width of the annotation color bar area
  annotation_legend = TRUE, 
  
  # legend_breaks can help if the Pearson scale is too wide
  legend = TRUE,
 
  # Manually set the width of the row annotation area if needed
  # (pheatmap doesn't have a direct 'margin' argument, so we optimize the layout)
  cellwidth = 25,    # Smaller cells give more room to the margins
  cellheight = 9
)
dev.off()
