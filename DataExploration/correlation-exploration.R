load("../4-IDs-metadata/metadata.RData")

biblioteca_rnaseq <- data.frame(
  sample_id = metadata$sample_id,
  barcode = metadata$barcode,
  reads = metadata$reads
)

df_data_clones <- read.csv("df_data_clones.csv")
df_data_reads <- read.csv("df_data_reads.csv")

biblioteca_rnaseq_amostrasFaltantes <- data.frame(
  sample_id = df_data_clones$sample_id[c(1,2,78)],
  barcode = df_data_clones$TCGA_barcode[c(1,2,78)],
  reads = c(60343619, 34138897, 49076175)
)

biblioteca_rnaseq_full <- rbind(biblioteca_rnaseq, biblioteca_rnaseq_amostrasFaltantes)

idx <- match(df_data_clones$sample_id, biblioteca_rnaseq_full$sample_id)
df_data_clones$reads <- biblioteca_rnaseq_full$reads[idx]

idx <- match(df_data_reads$sample_id, biblioteca_rnaseq_full$sample_id)
df_data_reads$reads <- biblioteca_rnaseq_full$reads[idx]

df_data_clones$tcr_counts <- rowSums(df_data_clones[,c("TRA","TRB","TRD","TRG")])
df_data_clones$bcr_counts <- rowSums(df_data_clones[,c("IGH","IGK","IGL")])
df_data_clones$sum_bcr_tcr <- rowSums(df_data_clones[,c("tcr_counts","bcr_counts")])

df_data_reads$tcr_counts <- rowSums(df_data_reads[,c("TRA","TRB","TRD","TRG")])
df_data_reads$bcr_counts <- rowSums(df_data_reads[,c("IGH","IGK","IGL")])
df_data_reads$sum_bcr_tcr <- rowSums(df_data_reads[,c("tcr_counts","bcr_counts")])

data_cor <- df_data_clones[,c("IGH","IGK","IGL","TRA","TRB","TRD","TRG",
                              "bcr_counts","tcr_counts","sum_bcr_tcr","reads")]

library(corrplot)

# Calcular correlação Spearman
cor_matrix <- cor(data_cor, method = "spearman", use = "pairwise.complete.obs")

# Exibir matriz de correlação
print(cor_matrix)

corrplot(cor_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)

library(Hmisc)

# Calcular matriz de correlação Spearman com p-valores
cor_test <- rcorr(as.matrix(data_cor), type = "spearman")

# Extrair matriz de correlação e p-valores
cor_matrix <- cor_test$r
p_matrix <- cor_test$P

library(ggcorrplot)

# Gerar gráfico de correlação com p-valores
pdf("correlation_plot.pdf", width = 8, height = 6)  # Define tamanho do PDF
ggcorrplot(cor_matrix, method = "square", type = "upper",
           p.mat = p_matrix, sig.level = 0.05, lab = TRUE)
dev.off()

write.csv(df_data_reads, file = "df_data_reads.csv")
write.csv(df_data_clones, file = "df_data_clones.csv")
