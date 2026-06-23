library(limma)
library(ggplot2)
library(glmnet)
library(dplyr)
library(tidyverse)
library(enrichR)
library(patchwork)
library(zCompositions)
library(ggsci)

setwd("C:/Users/Anna/Documents/pancancer/korean_3CA_tisch_2")
data <- read.csv("cell_composition_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
prop <- sweep(data, 2, colSums(data), "/")
prop <- prop + min(prop[prop>0])/2
clr_mat <- t(as.matrix(clr(t(prop))))

expr_counts <- read.csv("pseudobulk_malignant.csv", header = TRUE, row.names = 1, check.names = FALSE)
genes_464 <- read.table("genes_464.txt")
expr_counts <- expr_counts[rownames(expr_counts) %in% genes_464$x,]
metadata <- read.csv("metadata_all.csv", row.names=2)
metadata <- metadata[colnames(expr_counts), ]

cluster_df <- read.csv("./CLR_nonzero/cluster_df_464.csv", header = TRUE, row.names = 1, 
                       check.names = FALSE)
cluster_df$Sample <- trimws(cluster_df$Sample)

expr_counts_log <- log2(expr_counts + 1)
expr_norm <- normalizeBetweenArrays(expr_counts_log, method="quantile")

clr_mat <- clr_mat[,colnames(expr_counts)]
tcell_prop <- as.numeric(clr_mat["Tem/Trm cytotoxic T cells", ])
metadata$tcell_prop <- tcell_prop

# Pathology
phathology <- read.table(unz("C:/Users/Anna/Documents/pancancer/korean_3CA_thisch/cancer_data.tsv.zip", "cancer_data.tsv"), sep = "\t", header = T)
phathology <- phathology[phathology$Cancer!="lymphoma",]
unknow_genes <- row.names(expr_norm)[!(row.names(expr_norm) %in% phathology$Gene.name)]

phathology_res <- phathology %>%
  group_by(Gene.name) %>%
  summarise(
    sum_high = sum(High, na.rm = TRUE),
    sum_medium = sum(Medium, na.rm = TRUE),
    sum_low = sum(Low, na.rm = TRUE),
    sum_not_detected = sum(Not.detected, na.rm = TRUE)
  ) %>%
  mutate(
    statistic = ifelse(
      is.na(sum_high) | is.na(sum_medium) | is.na(sum_low) | is.na(sum_not_detected),
      NA,
      ((sum_high + sum_medium + sum_low) / (sum_high + sum_medium + sum_low + sum_not_detected))*100
    )
  ) %>%
  dplyr::select(Gene.name, statistic)

immune_genes<-phathology_res[is.na(phathology_res$statistic) == FALSE & phathology_res$statistic < 5,]$Gene.name

gene_stats <- read.csv("malignant_gene_expression_stats.csv", check.names = FALSE, header = TRUE)

hist(gene_stats$pct_cells_expressed, breaks = 100)
hist(gene_stats[gene_stats$pct_cells_expressed < 1, ]$pct_cells_expressed, breaks = 100)
hist(log(gene_stats$mean_expression), breaks = 100)

banned_genes <- gene_stats[(gene_stats$pct_cells_expressed <= 0.3 |
                              log(gene_stats$mean_expression) < -5) |
                             gene_stats$gene %in% immune_genes,]$gene
nrow(expr_norm[!(rownames(expr_norm) %in% banned_genes), ])
expr_norm <- expr_norm[!(rownames(expr_norm) %in% banned_genes), ]

# 4. ELASTIC NET
cancer_types <- as.factor(metadata$unified_cancer_type)
design_cancer <- model.matrix(~ cancer_types - 1)
colnames(design_cancer) <- gsub("cancer_types", "", colnames(design_cancer))
x_combined <- cbind(t(expr_norm), design_cancer)
y <- metadata$tcell_prop
penalty_factor <- c(
  rep(1, ncol(t(expr_norm))),
  rep(0, ncol(design_cancer))
)
set.seed(42)
fit <- cv.glmnet(x_combined, y, alpha = 0.5,
                 penalty.factor = penalty_factor)

coef_mat <- as.matrix(coef(fit, s = "lambda.min"))
coef_df <- data.frame(
  feature = rownames(coef_mat),
  beta = coef_mat[, 1],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
coef_df <- coef_df[coef_df$feature != "(Intercept)", ]
coef_df <- coef_df[!(coef_df$feature %in% cancer_types), ]
coef_df <- coef_df[coef_df$beta != 0, ]
coef_df <- coef_df[order(abs(coef_df$beta), decreasing = TRUE), ]

signature_genes <- coef_df$feature
signature_weights <- coef_df$beta
names(signature_weights) <- signature_genes

signature_score <- colSums(expr_norm[signature_genes, ] * signature_weights)
metadata$signature_score <- signature_score

# Visualizing
# CD8 T
df_cor <- data.frame(
  signature_score = signature_score,
  tcell_prop = metadata$tcell_prop
)

p1 <- ggplot(df_cor, aes(signature_score, tcell_prop)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", color = "purple3", se = TRUE) +
  annotate("text", 
           x = min(signature_score) + 0.1 * diff(range(signature_score)),
           y = max(tcell_prop) - 0.05 * diff(range(tcell_prop)),
           label = paste("r =", round(cor(signature_score, tcell_prop), 3)),
           size = 5, hjust = 0) +
  theme_bw() +
  labs(title = "Signature Score vs CD8 T cell Proportion",
       x = "Signature Score", 
       y = "CD8 T cell Proportion (CLR)")

print(p1)
ggsave("./CLR_nonzero/signature_score.png", p1, width = 5, height = 5)

# =========== 5. Clinical variables ===========
clinical_vars <- c('unified_cancer_type', 'source')

metadata$signature_score <- signature_score

plot_list <- list()

for (var in clinical_vars) {
  if (var %in% colnames(metadata)) {
    temp_df <- data.frame(
      signature_score = metadata$signature_score,
      var_value = metadata[[var]],
      row_names = rownames(metadata),
      stringsAsFactors = FALSE
    )
    temp_df <- temp_df[!is.na(temp_df$var_value) & temp_df$var_value != "", ]
    if (nrow(temp_df) > 0 && length(unique(temp_df$var_value)) > 1) {
      p <- ggplot(temp_df, aes(x = var_value, y = signature_score, fill = var_value)) +
        geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.title = element_text(size = 11, face = "bold"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "none"
        ) +
        labs(
          x = gsub("_", " ", var),
          y = "Signature Score"
        ) +
        scale_fill_viridis_d()
      
      if (length(unique(temp_df$var_value)) >= 2) {
        anova_test <- aov(signature_score ~ var_value, data = temp_df)
        p_value <- summary(anova_test)[[1]][["Pr(>F)"]][1]
        
        if (p_value < 0.0001) {
          p_label <- format(p_value, scientific = TRUE, digits = 3)
        } else {
          p_label <- round(p_value, 4)
        }
        
        p <- p + annotate("text", 
                          x = 0.5, 
                          y = max(temp_df$signature_score, na.rm = TRUE) * 1.05,
                          label = paste("ANOVA p =", p_label),
                          size = 3.5,
                          hjust = 0)
      }
      
      plot_list[[var]] <- p
      
      stats <- temp_df %>%
        group_by(var_value) %>%
        summarise(
          n = n(),
          mean_score = mean(signature_score, na.rm = TRUE),
          sd_score = sd(signature_score, na.rm = TRUE),
          median_score = median(signature_score, na.rm = TRUE)
        )
    }
  }
}

n_plots <- length(plot_list)
n_cols <- ifelse(n_plots >= 3, 3, n_plots)

p_combined <- wrap_plots(plot_list, ncol = n_cols, guides = "collect")
print(p_combined)

ggsave("./CLR_nonzero/signature_by_vars.png", p_combined, width = 15, height = 10)
  
# 5.3. Clusters
cluster_df$Sample_clean <- trimws(cluster_df$Sample)
metadata$Sample_clean <- trimws(rownames(metadata))

merged_data <- merge(metadata, cluster_df, by.x = "Sample_clean", by.y = "Sample_clean")
merged_data$clusters <- as.factor(merged_data$clusters)

p_cluster <- ggplot(merged_data, aes(x = clusters, y = signature_score, fill = clusters)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  scale_fill_nejm() +
  theme_bw() +
  labs(x = "Cluster",
       y = "Signature Score") +
  guides(fill = FALSE)

print(p_cluster)
ggsave("./CLR_nonzero/signature_by_cluster.png", p_cluster, width = 7, height = 4)

cluster_stats <- merged_data %>%
  group_by(clusters) %>%
  summarise(
    n = n(),
    mean_score = mean(signature_score, na.rm = TRUE),
    sd_score = sd(signature_score, na.rm = TRUE),
    median_score = median(signature_score, na.rm = TRUE),
    min_score = min(signature_score, na.rm = TRUE),
    max_score = max(signature_score, na.rm = TRUE)
  )
print(cluster_stats)

if (length(unique(merged_data$clusters)) > 1) {
  anova_res <- aov(signature_score ~ clusters, data = merged_data)
  summary(anova_res)
}

# 6.  ========= ENRICHMENT ANALYSIS ==============
pos_genes <- coef_df$feature[coef_df$beta > 0]
neg_genes <- coef_df$feature[coef_df$beta < 0]

background_genes <- rownames(expr_norm)
background_genes <- unique(na.omit(background_genes))
write.table(background_genes, "./CLR_nonzero/background_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

run_enrichment <- function(genes, background_genes, name_prefix = "") {

  dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
           "GO_Biological_Process_2021",
           "KEGG_2026",
           "Reactome_2022")
  
  enrich_res <- enrichr(genes, dbs, background = background_genes)
  
  res <- do.call(rbind, lapply(names(enrich_res), function(db_name) {
    df <- enrich_res[[db_name]]
    if (!is.null(df) && nrow(df) > 0) {
      df$Database <- db_name
      return(df)
    }
    return(NULL)
  }))
  
  overlap_counts <- sapply(res$Genes, function(x) {
    genes_list <- strsplit(as.character(x), ";")[[1]]
    return(length(genes_list))
  })
  res <- res[overlap_counts >= 2, ]
  res <- res[res$P.value < 0.01, ]
  res <- res[order(res$P.value), ]
  
  message(paste("Found", nrow(res), "terms", name_prefix))
  
  return(res)
}

pos_enrich <- run_enrichment(pos_genes, background_genes, "positive")

neg_enrich <- run_enrichment(neg_genes, background_genes, "negative")

# ===========  7. Genes networks =========== 
create_network <- function(enrich_res, genes, betas, name = "") {
  if (is.null(enrich_res) || nrow(enrich_res) == 0) {
    message(paste("Нет данных для сети:", name))
    return(NULL)
  }
  
  term2genes <- strsplit(as.character(enrich_res$Genes), ";")
  names(term2genes) <- enrich_res$Term
  
  edges <- data.frame(
    source = character(),
    target = character(),
    term_type = character(),
    p_value = numeric(),
    adjusted_p = numeric(),
    genes_count = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(enrich_res)) {
    term <- enrich_res$Term[i]
    genes_in_term <- term2genes[[term]]
    
    genes_in_term <- genes_in_term[genes_in_term %in% genes]
    
    if (length(genes_in_term) == 0) next
    
    term_info <- enrich_res[i, ]
    
    for (gene in genes_in_term) {
      edges <- rbind(edges, data.frame(
        source = term,
        target = gene,
        term_type = term_info$Database,
        p_value = term_info$P.value,
        adjusted_p = ifelse("Adjusted.P.value" %in% colnames(term_info), 
                            term_info$Adjusted.P.value, NA),
        genes_count = length(genes_in_term),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(edges) == 0) {
    message(paste("Нет связей для сети:", name))
    return(NULL)
  }
  
  term_nodes <- data.frame(
    id = enrich_res$Term,
    type = "term",
    beta = NA,
    sign = NA,
    database = enrich_res$Database,
    p_value = enrich_res$P.value,
    adjusted_p = ifelse("Adjusted.P.value" %in% colnames(enrich_res), 
                        enrich_res$Adjusted.P.value, NA),
    neg_log_p = -log10(enrich_res$P.value),
    genes_count = sapply(strsplit(as.character(enrich_res$Genes), ";"), length),
    stringsAsFactors = FALSE
  )
  
  genes_in_edges <- unique(edges$target)
  gene_nodes <- data.frame(
    id = genes_in_edges,
    type = "gene",
    beta = betas[genes_in_edges],
    sign = ifelse(betas[genes_in_edges] > 0, "positive", "negative"),
    database = NA,
    p_value = NA,
    adjusted_p = NA,
    neg_log_p = NA,
    genes_count = NA,
    stringsAsFactors = FALSE
  )
  
  nodes <- rbind(gene_nodes, term_nodes)
  
  nodes$size <- ifelse(nodes$type == "gene", 
                       abs(nodes$beta) * 10 + 20,
                       nodes$neg_log_p * 2 + 30)
  
  nodes$color <- ifelse(nodes$type == "term",
                        ifelse(nodes$database == "GO_Biological_Process_2021", "#95D7AE",
                               ifelse(nodes$database == "KEGG_2026", "#F4D03F",
                                      ifelse(nodes$database == "Reactome_2022", "#E8A87C",
                                             ifelse(nodes$database == "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "#C38D9E", "#999999")))),
                        "#CCCCCC")
  
  nodes[is.na(nodes)] <- ""
  
  return(list(nodes = nodes, edges = edges))
}

create_combined_network <- function(pos_enrich, neg_enrich, pos_genes, neg_genes, 
                                    betas, color_genes = TRUE) {
  
  combined_enrich <- rbind(pos_enrich, neg_enrich)
  
  if (is.null(combined_enrich) || nrow(combined_enrich) == 0) {
    message("Нет данных для комбинированной сети")
    return(NULL)
  }
  
  combined_enrich <- combined_enrich[order(combined_enrich$P.value), ]
  combined_enrich <- combined_enrich[!duplicated(combined_enrich$Term), ]
  
  term2genes <- strsplit(as.character(combined_enrich$Genes), ";")
  names(term2genes) <- combined_enrich$Term
  
  edges <- data.frame(
    source = character(),
    target = character(),
    term_type = character(),
    p_value = numeric(),
    adjusted_p = numeric(),
    genes_count = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(combined_enrich)) {
    term <- combined_enrich$Term[i]
    genes_in_term <- term2genes[[term]]
    
    all_genes <- unique(c(pos_genes, neg_genes))
    genes_in_term <- genes_in_term[genes_in_term %in% all_genes]
    
    if (length(genes_in_term) == 0) next
    
    term_info <- combined_enrich[i, ]
    
    for (gene in genes_in_term) {
      edges <- rbind(edges, data.frame(
        source = term,
        target = gene,
        term_type = term_info$Database,
        p_value = term_info$P.value,
        adjusted_p = ifelse("Adjusted.P.value" %in% colnames(term_info), 
                            term_info$Adjusted.P.value, NA),
        genes_count = length(genes_in_term),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(edges) == 0) {
    message("Нет связей для комбинированной сети")
    return(NULL)
  }
  
  term_nodes <- data.frame(
    id = combined_enrich$Term,
    type = "term",
    beta = NA,
    sign = NA,
    database = combined_enrich$Database,
    p_value = combined_enrich$P.value,
    adjusted_p = ifelse("Adjusted.P.value" %in% colnames(combined_enrich), 
                        combined_enrich$Adjusted.P.value, NA),
    neg_log_p = -log10(combined_enrich$P.value),
    genes_count = sapply(strsplit(as.character(combined_enrich$Genes), ";"), length),
    stringsAsFactors = FALSE
  )
  
  genes_in_edges <- unique(edges$target)
  
  gene_nodes <- data.frame(
    id = genes_in_edges,
    type = "gene",
    beta = betas[genes_in_edges],
    sign = ifelse(betas[genes_in_edges] > 0, "positive", "negative"),
    database = NA,
    p_value = NA,
    adjusted_p = NA,
    neg_log_p = NA,
    genes_count = NA,
    stringsAsFactors = FALSE
  )
  
  nodes <- rbind(gene_nodes, term_nodes)
  
  nodes$size <- ifelse(nodes$type == "gene", 
                       abs(nodes$beta) * 10 + 20,
                       nodes$neg_log_p * 2 + 30)
  
  if (color_genes) {
    nodes$color <- ifelse(nodes$type == "gene",
                          ifelse(nodes$sign == "positive", "#FF6B6B", "#4D9DE0"),
                          ifelse(nodes$database == "GO_Biological_Process_2021", "#95D7AE",
                                 ifelse(nodes$database == "KEGG_2026", "#F4D03F",
                                        ifelse(nodes$database == "Reactome_2022", "#E8A87C",
                                               ifelse(nodes$database == "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "#C38D9E", "#999999")))))
  } else {
    nodes$color <- ifelse(nodes$type == "gene",
                          "#CCCCCC",
                          ifelse(nodes$database == "GO_Biological_Process_2021", "#95D7AE",
                                 ifelse(nodes$database == "KEGG_2026", "#F4D03F",
                                        ifelse(nodes$database == "Reactome_2022", "#E8A87C",
                                               ifelse(nodes$database == "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "#C38D9E", "#999999")))))
  }
  
  nodes[is.na(nodes)] <- ""
  
  return(list(nodes = nodes, edges = edges))
}

# =============================================================================
pos_network <- create_network(pos_enrich, pos_genes, signature_weights, "positive")
write.table(pos_network$nodes, "./CLR_nonzero/nodes_positive_network.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pos_network$edges, "./CLR_nonzero/edges_positive_network.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


neg_network <- create_network(neg_enrich, neg_genes, signature_weights, "negative")
write.table(neg_network$nodes, "./CLR_nonzero/nodes_negative_network.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(neg_network$edges, "./CLR_nonzero/edges_negative_network.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

combined_network_colored <- create_combined_network(
  pos_enrich, neg_enrich, 
  pos_genes, neg_genes, 
  signature_weights,
  color_genes = TRUE
)

write.table(combined_network_colored$nodes, "./CLR_nonzero/nodes_combined_colored_network.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(combined_network_colored$edges, "./CLR_nonzero/edges_combined_colored_network.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
enrichment_summary <- data.frame()

pos_summary <- data.frame(
  Gene_Set = "Positive",
  Term = pos_enrich$Term,
  Database = pos_enrich$Database,
  P_value = pos_enrich$P.value,
  Adjusted_P = ifelse("Adjusted.P.value" %in% colnames(pos_enrich), 
                      pos_enrich$Adjusted.P.value, NA),
  Genes = pos_enrich$Genes,
  Gene_Count = sapply(strsplit(as.character(pos_enrich$Genes), ";"), length),
  stringsAsFactors = FALSE
)
enrichment_summary <- rbind(enrichment_summary, pos_summary)

neg_summary <- data.frame(
  Gene_Set = "Negative",
  Term = neg_enrich$Term,
  Database = neg_enrich$Database,
  P_value = neg_enrich$P.value,
  Adjusted_P = ifelse("Adjusted.P.value" %in% colnames(neg_enrich), 
                      neg_enrich$Adjusted.P.value, NA),
  Genes = neg_enrich$Genes,
  Gene_Count = sapply(strsplit(as.character(neg_enrich$Genes), ";"), length),
  stringsAsFactors = FALSE
)
enrichment_summary <- rbind(enrichment_summary, neg_summary)

write.table(enrichment_summary, "./CLR_nonzero/enrichment_summary_all.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
# 8. Save
# =============================================================================

gene_coef_final <- coef_df
gene_stats_final <- gene_stats[gene_stats$gene %in% gene_coef_final$feature, ]
colnames(gene_coef_final)[colnames(gene_coef_final) == "feature"] <- "gene"
combined_final <- merge(gene_coef_final, gene_stats_final, by.x = "gene", by.y = "gene")

write.table(combined_final, "./CLR_nonzero/Elastic_penalized_betas_stats_CLR.txt", 
            quote = FALSE, row.names = FALSE)
write.table(pos_genes, "./CLR_nonzero/Elastic_penalized_pos_CLR.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(neg_genes, "./CLR_nonzero/Elastic_penalized_neg_CLR.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

gene_cor <- sapply(signature_genes, function(g) {
  cor(expr_norm[g, ], metadata$tcell_prop)
})

df_cor <- data.frame(
  gene = signature_genes,
  cor = gene_cor,
  beta = signature_weights
)

cor_value <- cor(signature_score, metadata$tcell_prop)

p_cor <- ggplot(df_cor, aes(x = reorder(gene, cor), y = cor, fill = beta > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = cor_value, color = "red", linewidth = 1.2) +
  annotate("text", 
           x = 0.8, 
           y = cor_value,
           label = paste("r =", round(cor_value, 3)),
           hjust = 1.1, 
           vjust = -0.5,
           color = "black", 
           size = 5) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c("TRUE" = "#FF6B6B", "FALSE" = "#4D9DE0")) +
  ylab("Correlation with CD8 proportion") +
  xlab("Gene") +
  theme(axis.text.y = element_blank(),
        legend.position = "none")

ggsave("./CLR_nonzero/gene_correlations.png", p_cor, width = 10, height = 8)

signature_results <- data.frame(
  sample = names(signature_score),
  signature_score = signature_score,
  tcell_prop = metadata$tcell_prop,
  cluster = merged_data$clusters[match(names(signature_score), merged_data$Sample_clean)]
)
write.csv(signature_results, "./CLR_nonzero/signature_scores.csv", row.names = FALSE)
