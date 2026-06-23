library(dplyr)
library(ConsensusClustering)
library(pheatmap)
library(factoextra)
library(ggpubr)
library(uwot)
library(tidyr)
library(cluster)
library(ggsci)


setwd("C:/Users/Anna/Documents/pancancer/korean_3CA_tisch_2/")
load("consensus_matrix_CRC_ImmuneFilter_compact_clr_temp2.RData")

# ==========
LogitScore = Result[["LogitScore"]]
PAC = Result[["PAC"]]
deltaA = Result[["deltaA"]]
CMavg = Result[["CMavg"]]
Kopt = Result[["Kopt_LogitScore"]]
message(paste0("The optimum number of clusters = ", Kopt))

pheatmap1<-pheatmap::pheatmap(CM[[Kopt]], 
                              show_colnames = FALSE,
                              show_rownames = FALSE)
ggsave("./CLR_nonzero/Heatmap_Kopt.png", pheatmap1, width = 6, height = 5)

clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = Kopt)
clusters_df <- as.data.frame(clusters)
rownames(clusters_df) <- colnames(prop)
clusters_df$Sample <- row.names(clusters_df)
write.csv(clusters_df, "./CLR_nonzero/cluster_df_464.csv")

cluster_labels <- clusters_df$clusters
names(cluster_labels) <- clusters_df$Sample

d <- dist(clr_mat)
sil <- silhouette(
  cluster_labels,
  d
)
mean(sil[,3])

cluster_palette <- pal_nejm("default")(Kopt)
names(cluster_palette) <- 1:Kopt

clr_mat <- t(clr_mat)
# heatmap
annotation_col <- data.frame(Cluster = as.factor(clusters_df$clusters))
rownames(annotation_col) <- clusters_df$Sample
annotation_colors1 <- list(Cluster = cluster_palette)

pheatmap1 <- pheatmap(
  clr_mat,
  show_colnames = FALSE,
  clustering_method = "ward.D2",
  treeheight_col = 100,
  treeheight_row = 50,
  cluster_rows = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors1)
ggsave("./CLR_nonzero/Heatmap_CLR.png", pheatmap1, width = 7, height = 5)

annotation_col_order <- data.frame(
  Cluster = clusters_df[order(clusters_df$clusters),]$clusters, 
  row.names = clusters_df[order(clusters_df$clusters),]$Sample
)
combined_order <- clr_mat[, rownames(annotation_col_order)]
pheatmap_order <- pheatmap(
  combined_order,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors1
)
ggsave("./CLR_nonzero/Heatmap_CLR_order.png", pheatmap_order, width = 7, height = 4)

# Heatmap prop

prop_order <- prop[, rownames(annotation_col_order)]
pheatmap_prop_order <- pheatmap(
  prop_order,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors1
)
ggsave("./CLR_nonzero/Heatmap_Proportions_order.png", pheatmap_prop_order, width = 7, height = 4)

# PCA
res.pca_clr <- prcomp(t(clr_mat))
res.pca_prop <- prcomp(t(prop))

create_pca_plots <- function(pca_obj, data_type, palette_colors) {
  pca_ind_df <- as.data.frame(pca_obj$x)
  pca_ind_df$Sample <- rownames(pca_ind_df)
  pca_ind_df$clusters <- clusters_df[match(rownames(pca_ind_df), clusters_df$Sample), "clusters"]
  cluster_levels <- as.numeric(names(palette_colors))  # 1,2,3,4,5
  pca_ind_df$clusters <- factor(pca_ind_df$clusters, levels = cluster_levels)
  
  pca_ind_plot <- fviz_pca_ind(
    pca_obj,
    geom = "point",
    col.ind = pca_ind_df$clusters,
    palette = palette_colors,
    title = paste0("PCA - Individuals (", data_type, ")"),
    repel = TRUE) +
    theme_minimal()
  pca_var_plot <- fviz_pca_var(
    pca_obj,
    labelsize = 3,
    repel = TRUE,
    title = paste0("PCA - Variables (", data_type, ")")
  ) +
    theme_minimal()
  return(list(ind = pca_ind_plot, var = pca_var_plot, df = pca_ind_df))
}

pca_clr <- create_pca_plots(res.pca_clr, "CLR", cluster_palette)
ggsave("./CLR_nonzero/pca_ind_CLR.png", pca_clr$ind, width = 6, height = 5)
ggsave("./CLR_nonzero/pca_var_CLR.png", pca_clr$var, width = 8, height = 6)

pca_prop <- create_pca_plots(res.pca_prop, "Proportions", cluster_palette)
ggsave("./CLR_nonzero/pca_ind_Proportions.png", pca_prop$ind, width = 6, height = 5)
ggsave("./CLR_nonzero/pca_var_Proportions.png", pca_prop$var, width = 8, height = 6)

# Umap

create_umap <- function(pca_df, data_type, palette_colors, n_neighbors = 100) {
  
  matrix_umap <- umap(
    pca_df[, 1:3], 
    n_neighbors = n_neighbors, 
    min_dist = 0.01, 
    verbose = TRUE, 
    seed = 123
  )
  
  umap_df <- data.frame(
    UMAP1 = as.data.frame(matrix_umap)$V1,
    UMAP2 = as.data.frame(matrix_umap)$V2, 
    pca_df, 
    stringsAsFactors = FALSE
  )
  
  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = clusters)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = palette_colors) +
    theme_minimal() +
    labs(
      color = "Subtypes"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(list(plot = umap_plot, df = umap_df))
}

umap_clr <- create_umap(pca_clr$df, "CLR", cluster_palette)
ggsave("./CLR_nonzero/Umap_CLR.png", umap_clr$plot, width = 8.5, height = 5)

umap_prop <- create_umap(pca_prop$df, "Proportions", cluster_palette)
ggsave("./CLR_nonzero/Umap_Proportions.png", umap_prop$plot, width = 8.5, height = 5)

# ===================================================
df_long_combined <- data.frame()

for(i in unique(clusters_df$clusters)) {
  sample_list <- clusters_df %>%
    filter(clusters == i) %>%
    pull(Sample)
  cell_df_subset <- prop[, colnames(prop) %in% sample_list, drop = FALSE]
  cell_df_subset$CellType <- rownames(cell_df_subset)
  
  df_long <- cell_df_subset %>%
    pivot_longer(cols = -CellType,
                 names_to = "Sample",
                 values_to = "Count") %>%
    mutate(Cluster = i)
  
  df_long_combined <- bind_rows(df_long_combined, df_long)
}

c21 <- c(
  "dodgerblue2", "#E31A1C", "green4", "#6A3D9A",
  "#FF7F00", "gold1", "skyblue2", "#FB9A99",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "purple1", "red1", "green3"
)

box_plot <- ggplot(df_long_combined, aes(x = CellType, y = Count*100, fill = CellType)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.8) +
  scale_fill_manual(values = c21) +
  facet_wrap(~ Cluster, scales = "fixed", ncol = 2) + 
  theme_minimal() +
  labs(x = "Cell Type",
       y = "Percentage (%) of cell types",
       fill = "Cell Type") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(
          size = 14, 
          face = "bold",
          margin = margin(b = 10)
        ),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        panel.spacing = unit(1, "cm"),
        plot.margin = margin(20, 20, 20, 20, unit = "pt"))

print(box_plot)
ggsave("./CLR_nonzero/Boxplot_cell_types_by_cluster.png", box_plot, width = 15, height = 9)


df_long_combined <- data.frame()

for(i in unique(clusters_df$clusters)) {
  sample_list <- clusters_df %>%
    filter(clusters == i) %>%
    pull(Sample)
  cell_df_subset <- as.data.frame(clr_mat[, colnames(clr_mat) %in% sample_list, drop = FALSE])
  cell_df_subset$CellType <- as.vector(rownames(cell_df_subset))
  
  df_long <- cell_df_subset %>%
    pivot_longer(cols = -CellType,
                 names_to = "Sample",
                 values_to = "Count") %>%
    mutate(Cluster = i)
  
  df_long_combined <- bind_rows(df_long_combined, df_long)
}
box_plot <- ggplot(df_long_combined, aes(x = CellType, y = Count*100, fill = CellType)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.8) +
  scale_fill_manual(values = c21) +
  facet_wrap(~ Cluster, scales = "fixed", ncol = 2) + 
  theme_minimal() +
  labs(x = "Cell Type",
       y = "Percentage (%) of cell types",
       fill = "Cell Type") +
  guides(fill = guide_legend(ncol = 1)) + 
  theme(text = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(
          size = 14, 
          face = "bold",
          margin = margin(b = 10)
        ),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        panel.spacing = unit(1, "cm"),
        plot.margin = margin(20, 20, 20, 20, unit = "pt"))

print(box_plot)
ggsave("./CLR_nonzero/Boxplot_cell_types_by_cluster_clr_score.png", box_plot, width = 15, height = 9)
