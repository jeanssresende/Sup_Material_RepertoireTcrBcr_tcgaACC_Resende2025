library(tidyr)
library(ggplot2)
library(dplyr)

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")

arquivos <- list.files("../7-diversity/metricsTrust4/results_pipeline_report/", 
                       pattern = "\\.tsv$", full.names = TRUE)

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_entropy <- data.frame(t(dados$Entropy))
  rownames(dados_entropy) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_entropy)
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
data_metrics$tumor_stage <- data_clones$tumor_stage


data <- data_metrics

rownames(data) <- data_clones$TCGA_barcode
head(data)

df <- data

# Remover valores "not reported"
df <- df %>% filter(tumor_stage != "not reported")

# Calcular N por estágio tumoral
count_data <- df %>%
  group_by(tumor_stage) %>%
  summarise(n = n())

# Aplicar Kruskal-Wallis nos dados brutos e armazenar p-valor
kruskal_results <- data.frame(Variável = character(), P_Valor = numeric())

for (var in c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")) {
  test <- kruskal.test(df[[var]] ~ df$tumor_stage)
  kruskal_results <- rbind(kruskal_results, data.frame(Variável = var, P_Valor = test$p.value))
}

# Transformação logarítmica APENAS para visualização
df_log <- df %>%
  mutate(across(c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"), ~ log10(. + 1)))

# Criar e salvar boxplots com escala logarítmica e valores p
for (var in c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")) {
  plot <- ggboxplot(df_log, x = "tumor_stage", y = var, 
                    color = "tumor_stage", palette = "jco") +
    stat_compare_means(method = "kruskal.test", label = "p.format") +  # Teste Kruskal-Wallis
    labs(title = paste("Boxplot de", var, "- Comparação por Estágio Tumoral"), 
         y = paste(var, "(log10)"), x = "Estágio Tumoral") +
    geom_text(data = count_data, aes(x = tumor_stage, y = max(df_log[[var]]) * 1.1, label = paste("N =", n)), 
              color = "black", size = 5)  # Adiciona N ao gráfico
  
  # Salvar gráfico
  ggsave(filename = paste0("boxplot_entropy_", var, "_estagios.pdf"), plot = plot, width = 7, height = 5)
  
  print(plot)  # Exibir o gráfico
}

# Exibir tabela com os valores p do Kruskal-Wallis nos dados brutos
print(kruskal_results)


library(dplyr)
library(rstatix)

# Aplicar teste de Dunn apenas nas variáveis com p < 0.05 no Kruskal-Wallis
significativas <- c("IGH", "IGK", "IGL", "TRB")  # Pegamos apenas as variáveis que deram p < 0.05

dunn_results <- list()

for (var in significativas) {
  dunn_results[[var]] <- df %>%
    rstatix::dunn_test(as.formula(paste(var, "~ tumor_stage")), p.adjust.method = "BH")
}

# Exibir os resultados ajustados
dunn_results
