setwd("C:/Users/Anna/Documents/pancancer/korean_3CA_tisch_2")

lr_data <- read.csv("liana_all.csv", header = TRUE, check.names = FALSE, row.names = 1)
lr_data <- lr_data[, colnames(lr_data) != ""]
lr_malignant <- lr_data[lr_data$source == "Malignant", ]

pancancer_sig <- read.table("C:/Users/Anna/Documents/pancancer/korean_3CA_tisch_2/CLR_nonzero/Elastic_penalized_betas_stats_CLR.txt", 
                            header = TRUE)

signature_genes <- pancancer_sig$gene

lr_malignant_sig <- lr_malignant[lr_malignant$ligand_complex %in% signature_genes, ]
lr_malignant_sig$beta <- pancancer_sig$beta[match(lr_malignant_sig$ligand_complex, pancancer_sig$gene)]
lr_malignant_sig <- lr_malignant_sig[!is.na(lr_malignant_sig$beta), ]

# ============

lr_malignant_unique <- lr_malignant_sig %>%
  group_by(ligand_complex, receptor_complex, target) %>%
  summarise(
    mean_lrscore = mean(lrscore, na.rm = TRUE),
    mean_lr_means = mean(lr_means, na.rm = TRUE),
    mean_cellphone_pvals = mean(cellphone_pvals, na.rm = TRUE),
    mean_expr_prod = mean(expr_prod, na.rm = TRUE),
    mean_scaled_weight = mean(scaled_weight, na.rm = TRUE),
    mean_lr_logfc = mean(lr_logfc, na.rm = TRUE),
    mean_spec_weight = mean(spec_weight, na.rm = TRUE),
    mean_specificity_rank = mean(specificity_rank, na.rm = TRUE),
    mean_magnitude_rank = mean(magnitude_rank, na.rm = TRUE),
    beta = first(beta),
    n_samples = n(),
    samples = paste(unique(sample_name), collapse = "; "),
    datasets = paste(unique(dataset), collapse = "; "),
    clusters = paste(unique(cluster_label), collapse = "; ")
  ) %>%
  ungroup() %>%
  arrange(desc(mean_lrscore))

print(head(lr_malignant_unique[, c("ligand_complex", "receptor_complex", "target", 
                                   "mean_lrscore", "beta", "n_samples", "mean_cellphone_pvals")], 20))

write.csv(lr_malignant_unique, "./CLR_nonzero/malignant_LR_interactions_unique.csv", row.names = FALSE)

ligand_summary <- lr_malignant_unique %>%
  group_by(ligand_complex) %>%
  summarise(
    n_interactions = n(),
    mean_lrscore = mean(mean_lrscore, na.rm = TRUE),
    max_lrscore = max(mean_lrscore, na.rm = TRUE),
    mean_beta = mean(beta, na.rm = TRUE),
    targets = paste(unique(target), collapse = ", ")
  ) %>%
  arrange(desc(mean_lrscore))

print(head(ligand_summary, 10))

library(ggplot2)

lr_malignant_sorted <- lr_malignant_unique %>%
  arrange(desc(beta > 0), desc(abs(beta)))

p <- ggplot(lr_malignant_sorted, 
            aes(x = reorder(paste(ligand_complex, "->", target), -abs(beta)), 
                y = mean_lrscore, 
                fill = beta > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#FF6B6B", "FALSE" = "#4D9DE0"),
                    labels = c("TRUE" = "Positive beta", "FALSE" = "Negative beta")) +
  labs(x = "Ligand -> Target", y = "Mean LR Score") +
  theme_bw() +
  theme(legend.position = "bottom")

print(p)
ggsave("./CLR_nonzero/top20_malignant_LR_interactions_sorted.png", p, width = 5, height = 7)
