gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "PTEN", "KRAS")
control    <- c(5.2, 7.8, 3.1, 9.4, 6.0, 4.5)
treated    <- c(8.9, 7.6, 6.7, 12.1, 2.3, 9.8)
chromosome <- c("17", "17", "7", "8", "10", "12")

fold_change   <- treated / control
status_char   <- ifelse(fold_change > 1.2, "upregulated",
                        ifelse(fold_change < 0.8, "downregulated", "stable"))
classification <- factor(status_char, levels = c("upregulated", "downregulated", "stable"))

gene_info <- list(
  names = gene_names,
  stats = list(control = control, treated = treated),
  source = "RNA-Seq Gene Expression Study - Control vs Treated"
)

selected_source <- gene_info$source
print(selected_source)

gene_info$date_analyzed <- "2026-06-06"

gene_df <- data.frame(
  gene_names     = gene_names,
  control        = control,
  treated        = treated,
  chromosome     = chromosome,
  fold_change    = round(fold_change, 2),
  classification = classification,
  stringsAsFactors = FALSE
)

chr17_subset <- gene_df[gene_df$chromosome == "17", ]
print(chr17_subset)

gene_df$log2_fc <- round(log2(gene_df$fold_change), 2)

sorted_indices <- order(gene_df$fold_change, decreasing = TRUE)
gene_df_sorted <- gene_df[sorted_indices, ]
print(gene_df_sorted)

write.csv(gene_df_sorted, file = "gene_df_sorted.csv", row.names = FALSE)

point_colors <- c("upregulated" = "red", "downregulated" = "blue", "stable" = "gray")

png(filename = "expression_plot.png", width = 800, height = 600)
plot(
  x = gene_df$control, 
  y = gene_df$treated,
  col = point_colors[gene_df$classification],
  pch = 19,
  cex = 1.5,
  xlab = "Control Expression",
  ylab = "Treated Expression",
  main = "Gene Expression: Control vs Treated"
)
legend("topleft", legend = levels(gene_df$classification), col = point_colors, pch = 19)
dev.off()