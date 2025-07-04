# Carregar pacotes necessários
library(dplyr)
library(ggplot2)
library(tidyr)      # Para pivot_longer
library(ggpubr)     # Para stat_compare_means (útil para o rótulo)
library(rstatix)    # Para adjust_pvalue e estatísticas

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")

data <- data_clones[,c("IGH","IGK","IGL","TRA","TRB","TRG","TRD",
                       "steroid","tumor_stage")]

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
  ggsave(filename = paste0("boxplot_", var, "_estagios.pdf"), plot = plot, width = 7, height = 5)
  
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
