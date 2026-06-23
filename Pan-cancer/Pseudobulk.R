library(dplyr)
setwd("C:/Users/Anna/Documents/pancancer/korean_3CA_thisch")

pseudobulk_korean <- read.csv("all_sum_matrix.txt", row.names = 1, check.names = FALSE)
pseudobulk_tisch <- read.csv("pseudobulk_malignant_3.csv", row.names = 1, check.names = FALSE)
pseudobulk_tisch <- round(pseudobulk_tisch, digits = 0)
pseudobulk_tisch <- as.data.frame(t(pseudobulk_tisch))
pseudobulk_3CA <- read.csv("Merged_Gene_Expression_Matrix_3.csv", row.names = 1, check.names = FALSE)

# Delete genes without detected expression in malignant
phathology <- read.table(unz("pathology.tsv.zip", "pathology.tsv"), sep = "\t", header = T)
phathology <- phathology[phathology$Cancer!="lymphoma",]
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
  select(Gene.name, statistic)

hist(phathology_res$statistic, breaks = 100)
immune_genes<-phathology_res[is.na(phathology_res$statistic) == FALSE & phathology_res$statistic < 5,]$Gene.name

pseudobulk_korean<-pseudobulk_korean[(rownames(pseudobulk_korean)%in%immune_genes)==FALSE,]
pseudobulk_tisch<-pseudobulk_tisch[(rownames(pseudobulk_tisch)%in%immune_genes)==FALSE,]
pseudobulk_3CA<-pseudobulk_3CA[(rownames(pseudobulk_3CA)%in%immune_genes)==FALSE,]

genes_korean <- rownames(pseudobulk_korean)
genes_tisch <- rownames(pseudobulk_tisch)
genes_3CA <- rownames(pseudobulk_3CA)
common_elements <- intersect(intersect(genes_korean, genes_tisch), genes_3CA)
common_elements[grep("^CXCL", common_elements)]
common_elements[grep("^CD", common_elements)]

pseudobulk_korean_f <- pseudobulk_korean[rownames(pseudobulk_korean) %in% common_elements,]
pseudobulk_tisch_f <- pseudobulk_tisch[rownames(pseudobulk_tisch) %in% common_elements,]
pseudobulk_3CA_f <- pseudobulk_3CA[rownames(pseudobulk_3CA) %in% common_elements,]

immune_genes <- c("CXCL13","CCL4","XCL1","CCL5","XCL2","CCL15","CCL22","CXCL16","CCL23",
         "CXCL9","CCL17","CXCL10","CXCL11","CCL16","CXCL8","CXCL6","CCL19","CCL18",
         "CXCL3","CCL2","CCL24","CXCL12","CCL14","CCL8","CCL21","CCL13","CXCL5",
         "CCL7","CXCL1","CCL20","CXCL2","CCL28","CXCL17","CXCL14","CX3CL1","CCL26",
         "IL1A","IL1B","IL23A","IL10","IL12A","TNF","IFNG", "IL6R", "TNFRSF1A", "CXCR4",
         "TLR3","TLR5","TLR10","TLR7","TLR2","TLR4","TLR1","TLR6", "C3","C4BPA","C4BPB",
         "C5","C7","C8G", "CTLA4","PDCD1","PDCD1LG2","CD274","LAG3","TIGIT","HAVCR2",
         "SIGLEC15","CD86","CD40","ICOS","HLA-DQA1","HLA-DOA","HLA-DPA1","HLA-DPB1",
         "HLA-DQB1","HLA-DRA","HLA-DQA2","HLA-DRB1","HLA-F","HLA-DRB5","HLA-DMB","HLA-E","HLA-B",
         "HLA-DMA","HLA-A","HLA-DOB","HLA-C")
setdiff(immune_genes, common_elements)


dfs <- list(pseudobulk_korean_f, pseudobulk_tisch_f, pseudobulk_3CA_f) %>% 
  map(~as.data.frame(.) %>% rownames_to_column("gene"))
big_mat <- reduce(dfs, full_join, by = "gene")
rownames(big_mat) <- big_mat$gene
big_mat <- big_mat %>% select(-gene)
# write.csv(big_mat, "pseudobulk_tisch_korean_3ca.csv")
# big_mat <- read.csv("pseudobulk_tisch_korean_3ca.csv", check.names = FALSE, row.names = 1)
binary_mat <- big_mat %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, 1)))
binary_mat <- binary_mat[rownames(binary_mat) %in% immune_genes,]
immun_genes_samples <- as.data.frame(colSums(binary_mat))
samples_464 <- rownames(immun_genes_samples)[immun_genes_samples$`colSums(binary_mat)`
                                             %in% c(80)]

big_mat_f <- big_mat[,colnames(big_mat) %in% samples_464]
big_mat_f <- na.omit(big_mat_f)
write.table(big_mat_f, "big_mat_464.txt", quote = FALSE)
# big_mat_f <- read.table("big_mat_464.txt", check.names = FALSE, row.names = 1)

setdiff(immune_genes, rownames(big_mat_f))
samples_464 <- colnames(big_mat_f)
write.table(samples_464, "samples_464.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)
samples_464_tisch <- colnames(pseudobulk_tisch)[colnames(pseudobulk_tisch) %in% samples_464]
samples_464_3ca <- colnames(pseudobulk_3CA)[colnames(pseudobulk_3CA) %in% samples_464]

# samples_464 <- setdiff(rownames(clusters_df), colnames(big_mat_f))
# metadata_korean <- read.csv("combined_h5ad.csv")
# metadata_3CA <- read.csv("merged_metadata_3CA.csv")
# metadata_tisch <- read.csv("samples_metadata_tisch.csv")
# 
# metadata_3CA_586_464_19 <- metadata_3CA[metadata_3CA$sample_id_final %in% samples_586_464,]
# metadata_tisch_586_464_103 <- metadata_tisch[metadata_tisch$sample_id_final %in% samples_586_464,]
# write.csv(metadata_3CA_586_464_19, "metadata_3CA_586_464_19.csv", row.names = FALSE)
# write.csv(metadata_tisch_586_464_103, "metadata_tisch_586_464_103.csv", row.names = FALSE)

