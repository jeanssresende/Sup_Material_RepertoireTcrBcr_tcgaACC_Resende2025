library(tidyr)
library(ggplot2)
library(dplyr)

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")
dir_metrics <- "metricsTrust4/results_pipeline_report/"

arquivos <- list.files(dir_metrics, pattern = "\\.tsv$", full.names = TRUE)

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados$Clonality))
  rownames(dados_abundance) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_abundance)
}

colnames(data_metrics) <- dados$chain

data_metrics[is.na(data_metrics)] <- 0

# adicionando os casos que o TRUST4 nao encontrou TCR e BCR
mat_add <- data.frame(IGH = rep(0,10),
                      IGK = rep(0,10), 
                      IGL = rep(0,10), 
                      TRA = rep(0,10), 
                      TRB = rep(0,10), 
                      TRG = rep(0,10), 
                      TRD = rep(0,10))

rownames(data_clones) <- data_clones$sample_id

rownames(mat_add) <- rownames(data_clones)[rownames(data_clones) %in% 
                                             rownames(data_metrics) == FALSE]

data_metrics <- rbind(data_metrics, mat_add)

idx <- match(data_clones$sample_id, rownames(data_metrics))
data_metrics <- data_metrics[idx,]

data_metrics$steroid <- data_clones$steroid

# Supondo que seu dataframe se chama 'data_metrics'
data_long_log <- data_metrics %>%
  pivot_longer(cols = c(IGH, IGK, IGL, TRA, TRB, TRG, TRD),
               names_to = "gene",
               values_to = "count") %>%
  mutate(log2_count = log2(count + 1))

# Realizar o teste de Mann-Whitney U e corrigir os p-valores
p_values_log <- data_long_log %>%
  group_by(gene) %>%
  summarize(p_value_log = wilcox.test(log2_count ~ steroid, data = cur_data())$p.value) %>%
  mutate(p_adjust_log = p.adjust(p_value_log, method = "BH"))

# Juntar os p-valores ajustados aos dados transformados
data_long_log <- data_long_log %>%
  left_join(p_values_log, by = "gene")

# Abrir o dispositivo PDF
pdf(file = "boxplot_comparacao_Clonality.pdf", width = 8, height = 6)

# Criar e imprimir os boxplots
for (current_gene in unique(data_long_log$gene)) {
  # Filtrar os dados para o gene atual
  data_gene_log <- data_long_log %>% filter(gene == current_gene)
  
  # Obter o p-valor ajustado para o gene atual
  current_p_value_adj_log <- unique(data_gene_log$p_adjust_log)
  
  # Criar o boxplot com os dados transformados e adicionar o p-valor ajustado
  p <- ggplot(data_gene_log, aes(x = steroid, y = log2_count, fill = steroid)) +
    geom_boxplot() +
    labs(title = paste("Boxplot de log2(", current_gene, " + 1)"),
         x = "Steroid",
         y = "log2(Count + 1)") +
    theme_minimal() +
    annotate("text", x = 1.5, y = max(data_gene_log$log2_count) * 1.1,
             label = paste("p.adj =", format.pval(current_p_value_adj_log, digits = 3)))
  
  print(p) # Imprimir o gráfico para o PDF
}

# Fechar o dispositivo PDF
dev.off()

print("Gráficos salvos em boxplot_comparacao_Clonality.pdf")

for (current_gene in unique(data_long_log$gene)) {
  # Filtrar os dados para o gene atual
  data_gene_log <- data_long_log %>% filter(gene == current_gene)
  
  # Obter o p-valor ajustado para o gene atual
  current_p_value_adj_log <- unique(data_gene_log$p_adjust_log)
  
  # Nome do arquivo para este gene
  filename <- paste0("boxplot_", current_gene, "_clonality.pdf")
  
  # Abrir o dispositivo PDF para este arquivo
  pdf(file = filename, width = 8, height = 6)
  
  # Criar o boxplot
  p <- ggplot(data_gene_log, aes(x = steroid, y = log2_count, fill = steroid)) +
    geom_boxplot() +
    labs(title = paste("Boxplot de log2(", current_gene, " + 1)"),
         x = "Steroid",
         y = "log2(Count + 1)") +
    theme_minimal() +
    annotate("text", x = 1.5, y = max(data_gene_log$log2_count) * 1.1,
             label = paste("p.adj =", format.pval(current_p_value_adj_log, digits = 3)))
  
  print(p) # Imprimir o gráfico para o PDF
  
  # Fechar o dispositivo PDF para este arquivo
  dev.off()
  
  print(paste("Gráfico de", current_gene, "salvo em", filename))
}


table(data_metrics$steroid)
