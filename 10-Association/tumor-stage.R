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

df <- df %>% filter(tumor_stage != "not reported")

# Criar os dois grupos
df <- df %>%
  mutate(tumor_group = case_when(
    tumor_stage %in% c("stage i", "stage ii") ~ "Grupo 1 (I+II)",
    tumor_stage %in% c("stage iii", "stage iv") ~ "Grupo 2 (III+IV)"
  ))

# Calcular N por grupo
count_data <- df %>%
  group_by(tumor_group) %>%
  summarise(n = n())

# Aplicar teste de Mann-Whitney e corrigir p-valor
adjusted_pvals <- data.frame(Variável = character(), P_Valor_Ajustado = numeric())

for (var in c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")) {
  test <- wilcox.test(df[[var]] ~ df$tumor_group)
  adjusted_pvals <- rbind(adjusted_pvals, data.frame(
    Variável = var, 
    P_Valor_Ajustado = p.adjust(test$p.value, method = "BH")
  ))
}

# Criar e salvar os boxplots com escala logarítmica, valores p ajustados e N
for (var in c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")) {
  plot <- ggboxplot(df, x = "tumor_group", y = var, 
                    color = "tumor_group", palette = "jco") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +  # Teste Mann-Whitney
    scale_y_log10() +  # Escala logarítmica
    labs(title = paste("Boxplot de", var, "- Comparação de Grupos"), 
         y = paste(var, "(log)"), x = "Grupo Tumoral") +
    geom_text(data = count_data, aes(x = tumor_group, y = max(df[[var]]) * 1.1, label = paste("N =", n)), 
              color = "black", size = 5)
  
  # Salvar gráfico
  ggsave(filename = paste0("boxplot_", var, "_grupos.png"), plot = plot, width = 7, height = 5)
  
  print(plot)  # Exibir o gráfico
}

# Exibir tabela com os valores p ajustados
print(adjusted_pvals)
