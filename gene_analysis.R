gene_names <- c("BRCA1", "TP53", "EGFR", "MYC", "PTEN", "KRAS")
control    <- c(5.2, 7.8, 3.1, 9.4, 6.0, 4.5)
treated    <- c(8.9, 7.6, 6.7, 12.1, 2.3, 9.8)

names(control) <- gene_names
names(treated) <- gene_names

fold_change <- treated / control

high_fc_genes <- fold_change[fold_change > 1.5]
print("Genes with Fold Change > 1.5:")
print(high_fc_genes)

expr_matrix <- cbind(Control = control, Treated = treated)
rownames(expr_matrix) <- gene_names
avg_expression <- rowMeans(expr_matrix)

status_char <- ifelse(fold_change > 1.2, "upregulated",
                      ifelse(fold_change < 0.8, "downregulated", "stable"))

gene_status <- factor(status_char, levels = c("upregulated", "downregulated", "stable"))

results_df <- data.frame(
  Control = control,
  Treated = treated,
  FoldChange = round(fold_change, 2),
  MeanExpr = round(avg_expression, 2),
  Status = gene_status
)

write.table(results_df, file = "results.txt", sep = "\t", quote = FALSE, col.names = NA)