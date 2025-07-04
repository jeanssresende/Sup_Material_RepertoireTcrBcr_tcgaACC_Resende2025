library(ggplot2)
library(dplyr)
library(survival)
library(scales)
library(gridExtra)

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")
load("../6-dataExploration/tcgaACC_pre_processed.RData")

data <- data_clones[,c("TCGA_barcode","IGH","IGK","IGL","TRA","TRB","TRG","TRD",
                       "steroid","tumor_stage","vital_status")]

idx <- match(data$TCGA_barcode, tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

head(data)

# Carregar pacotes necessários
library(survival)
library(survminer)
library(dplyr)

df <- data

# Filtrar dados, garantindo apenas casos com sobrevida reportada
#df_survival <- df %>% filter(tumor_stage != "not reported")

df_survival <- df

# Criar objeto de sobrevivência
surv_object <- Surv(time = df_survival$OS.Time, event = df_survival$OS)

# Criar modelo Kaplan-Meier para comparar grupos de steroid
km_fit <- survfit(surv_object ~ steroid, data = df_survival)

# Aplicar teste log-rank para verificar diferenças estatísticas
log_rank_test <- survdiff(surv_object ~ steroid, data = df_survival)

# Criar gráfico Kaplan-Meier
plot_km <- ggsurvplot(km_fit, data = df_survival, 
                      pval = TRUE, conf.int = TRUE,
                      risk.table = TRUE, # Exibir tabela de risco
                      xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                      title = "Curva de Sobrevida Kaplan-Meier - Steroid",
                      legend.title = "Grupo Steroid",
                      legend.labs = c("Steroid High", "Steroid Low"),
                      palette = c("#E69F00", "#56B4E9"))

# Exibir gráfico
print(plot_km)

pdf("Kaplan_Meier_Steroid.pdf", width = 8, height = 6)
print(plot_km)
dev.off()

# Exibir teste log-rank
print(log_rank_test)

# Salvar o gráfico em PDF
ggsave("kaplan_meier_steroid.pdf", plot = plot_km$plot, width = 8, height = 6, device = "pdf")




# Filtrar apenas pacientes do estágio II
df_stage_ii <- df %>% filter(tumor_stage == "stage ii")

# Criar objeto de sobrevivência
surv_object_ii <- Surv(time = df_stage_ii$OS.Time, event = df_stage_ii$OS)

# Criar modelo Kaplan-Meier apenas para estágio II
km_fit_ii <- survfit(surv_object_ii ~ steroid, data = df_stage_ii)

# Aplicar teste log-rank para verificar diferenças estatísticas
log_rank_test_ii <- survdiff(surv_object_ii ~ steroid, data = df_stage_ii)

# Criar gráfico Kaplan-Meier
plot_km_ii <- ggsurvplot(km_fit_ii, data = df_stage_ii, 
                         pval = TRUE, conf.int = TRUE,
                         risk.table = TRUE,  # Exibir tabela de risco
                         xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                         title = "Curva de Sobrevida Kaplan-Meier - Estágio II",
                         legend.title = "Grupo Steroid",
                         legend.labs = c("Steroid High", "Steroid Low"),
                         palette = c("#E69F00", "#56B4E9"))




# Exibir gráfico
print(plot_km_ii)



# Salvar gráfico em PDF
ggsave("kaplan_meier_stage_ii.pdf", plot = plot_km_ii$plot, width = 8, height = 6, device = "pdf")

# Exibir teste log-rank
print(log_rank_test_ii)






library(survival)
library(survminer)
library(dplyr)

# Filtrar pacientes do estágio II
df_stage_ii <- df %>% filter(tumor_stage == "stage ii")

# Definir grupos de alta e baixa expressão (mediana)
df_stage_ii <- df_stage_ii %>%
  mutate(IGH_group = ifelse(IGH > median(IGH), "Alta", "Baixa"),
         IGK_group = ifelse(IGK > median(IGK), "Alta", "Baixa"),
         IGL_group = ifelse(IGL > median(IGL), "Alta", "Baixa"))

df_stage_ii <- df_stage_ii %>%
  mutate(TRA_group = ifelse(TRA > median(TRA), "Alta", "Baixa"),
         TRB_group = ifelse(TRB > median(TRB), "Alta", "Baixa"),
         TRG_group = ifelse(TRG > median(TRG), "Alta", "Baixa"),
         TRD_group = ifelse(TRD > median(TRD), "Alta", "Baixa"))

# Criar Kaplan-Meier para IGH
surv_object_IGH <- Surv(time = df_stage_ii$OS.Time, event = df_stage_ii$OS)
km_fit_IGH <- survfit(surv_object_IGH ~ IGH_group, data = df_stage_ii)
plot_km_IGH <- ggsurvplot(km_fit_IGH, data = df_stage_ii, pval = TRUE, conf.int = TRUE,
                          risk.table = TRUE, xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                          title = "Kaplan-Meier - IGH (Estágio II)", legend.title = "Grupo IGH",
                          legend.labs = c("Baixa Expressão", "Alta Expressão"), palette = c("#E69F00", "#56B4E9"))

# Exibir gráfico
print(plot_km_IGH)

# Aplicar teste log-rank
log_rank_test_IGH <- survdiff(surv_object_IGH ~ IGH_group, data = df_stage_ii)
print(log_rank_test_IGH)






# Criar Kaplan-Meier para IGH
surv_object_IGK <- Surv(time = df_stage_ii$OS.Time, event = df_stage_ii$OS)
km_fit_IGK <- survfit(surv_object_IGK ~ IGK_group, data = df_stage_ii)
plot_km_IGK <- ggsurvplot(km_fit_IGK, data = df_stage_ii, pval = TRUE, conf.int = TRUE,
                          risk.table = TRUE, xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                          title = "Kaplan-Meier - IGH (Estágio II)", legend.title = "Grupo IGK",
                          legend.labs = c("Baixa Expressão", "Alta Expressão"), palette = c("#E69F00", "#56B4E9"))

# Exibir gráfico
print(plot_km_IGK)

# Criar Kaplan-Meier para IGH
surv_object_IGL <- Surv(time = df_stage_ii$OS.Time, event = df_stage_ii$OS)
km_fit_IGL <- survfit(surv_object_IGL ~ IGL_group, data = df_stage_ii)
plot_km_IGL <- ggsurvplot(km_fit_IGL, data = df_stage_ii, pval = TRUE, conf.int = TRUE,
                          risk.table = TRUE, xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                          title = "Kaplan-Meier - IGH (Estágio II)", legend.title = "Grupo IGK",
                          legend.labs = c("Baixa Expressão", "Alta Expressão"), palette = c("#E69F00", "#56B4E9"))

# Exibir gráfico
print(plot_km_IGL)


# Criar Kaplan-Meier para IGH
surv_object_TRA <- Surv(time = df_stage_ii$OS.Time, event = df_stage_ii$OS)
km_fit_TRA <- survfit(surv_object_TRA ~ TRA_group, data = df_stage_ii)
plot_km_TRA <- ggsurvplot(km_fit_TRA, data = df_stage_ii, pval = TRUE, conf.int = TRUE,
                          risk.table = TRUE, xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                          title = "Kaplan-Meier - IGH (Estágio II)", legend.title = "Grupo IGK",
                          legend.labs = c("Baixa Expressão", "Alta Expressão"), palette = c("#E69F00", "#56B4E9"))

# Exibir gráfico
print(plot_km_TRA)


# Criar Kaplan-Meier para IGH
surv_object_TRB <- Surv(time = df_stage_ii$OS.Time, event = df_stage_ii$OS)
km_fit_TRB <- survfit(surv_object_TRB ~ TRB_group, data = df_stage_ii)
plot_km_TRB <- ggsurvplot(km_fit_TRB, data = df_stage_ii, pval = TRUE, conf.int = TRUE,
                          risk.table = TRUE, xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                          title = "Kaplan-Meier - IGH (Estágio II)", legend.title = "Grupo IGK",
                          legend.labs = c("Baixa Expressão", "Alta Expressão"), palette = c("#E69F00", "#56B4E9"))

# Exibir gráfico
print(plot_km_TRB)










# Criar modelo de Cox incluindo todas as cadeias
cox_model <- coxph(Surv(OS.Time, OS) ~ IGH + IGK + IGL + TRA + TRB + TRG + TRD, data = df_stage_ii)

# Resumo do modelo
summary(cox_model)


library(survival)

# Criar modelo de Cox incluindo steroid e todas as cadeias BCR/TCR
cox_model_steroid <- coxph(Surv(OS.Time, OS) ~ IGH + IGK + IGL + TRA + TRB + TRG + TRD + steroid, data = df)

# Resumo do modelo atualizado
summary(cox_model_steroid)


library(ggplot2)
library(dplyr)
library(survival)

# Filtrar apenas variáveis significativas (p < 0.05)
cox_summary <- summary(cox_model_steroid)
cox_results <- data.frame(
  Variável = rownames(cox_summary$coefficients),
  HR = exp(cox_summary$coefficients[, 1]),  # Hazard Ratio
  lower_CI = exp(cox_summary$conf.int[, 3]),  # Limite inferior do IC 95%
  upper_CI = exp(cox_summary$conf.int[, 4]),  # Limite superior do IC 95%
  p_valor = cox_summary$coefficients[, 5]  # Valor p
) #%>%
  #filter(p_valor < 0.05)  # Mantemos apenas variáveis significativas

# Criar gráfico de floresta com p-valores e eixo reduzido
ggplot(cox_results, aes(x = reorder(Variável, HR), y = HR)) +
  geom_point(aes(color = p_valor < 0.01), size = 3) +  # Pontos coloridos para destaque
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +  # Barras de erro
  geom_text(aes(label = paste0("p=", signif(p_valor, 3))), hjust = -0.5, size = 3) +  # Adicionar p-valor
  coord_flip(expand = FALSE) +  # Inverter e reduzir espaçamento
  #scale_y_log10(limits = c(0.5, 3)) +  # Escala logarítmica no eixo Y
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +  # Reduzir tamanho do texto do eixo Y
  labs(title = "Efeito das cadeias BCR/TCR e steroid na sobrevida",
       x = "Variáveis do modelo de Cox",
       y = "Hazard Ratio (HR) (escala log)") +
  scale_color_manual(values = c("red", "blue"), guide = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")  # Linha referência HR=1



# Criar modelo de Cox com interações entre steroid e as cadeias BCR/TCR
cox_model_interactions <- coxph(Surv(OS.Time, OS) ~ 
                                  IGH * steroid + IGK * steroid + IGL * steroid +
                                  TRA * steroid + TRB * steroid + TRG * steroid + TRD * steroid, 
                                data = df)

# Exibir resultados do modelo atualizado
summary(cox_model_interactions)


library(ggplot2)
library(dplyr)
library(survival)

# Filtrar apenas variáveis significativas (p < 0.05)
cox_summary <- summary(cox_model_interactions)
cox_results <- data.frame(
  Variável = rownames(cox_summary$coefficients),
  HR = exp(cox_summary$coefficients[, 1]),  # Hazard Ratio
  lower_CI = exp(cox_summary$conf.int[, 3]),  # Limite inferior do IC 95%
  upper_CI = exp(cox_summary$conf.int[, 4]),  # Limite superior do IC 95%
  p_valor = cox_summary$coefficients[, 5]  # Valor p
) %>%
  filter(p_valor < 0.05)  # Mantemos apenas variáveis significativas

# Criar gráfico de floresta com p-valores e eixo reduzido
ggplot(cox_results, aes(x = reorder(Variável, HR), y = HR)) +
  geom_point(aes(color = p_valor < 0.01), size = 3) +  # Pontos coloridos para destaque
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +  # Barras de erro
  geom_text(aes(label = paste0("p=", signif(p_valor, 3))), hjust = -0.5, size = 3) +  # Adicionar p-valor
  coord_flip(expand = FALSE) +  # Inverter e reduzir espaçamento
  scale_y_log10(limits = c(0.5, 3)) +  # Escala logarítmica no eixo Y
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +  # Reduzir tamanho do texto do eixo Y
  labs(title = "Efeito das cadeias BCR/TCR e steroid na sobrevida",
       x = "Variáveis do modelo de Cox",
       y = "Hazard Ratio (HR) (escala log)") +
  scale_color_manual(values = c("red", "blue"), guide = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")  # Linha referência HR=1


library(survival)

# Criar modelo de Cox considerando apenas as interações com steroid
cox_model_interactions_only <- coxph(Surv(OS.Time, OS) ~ 
                                       IGH:steroid + IGK:steroid + IGL:steroid +
                                       TRA:steroid + TRB:steroid + TRG:steroid + TRD:steroid, 
                                     data = df)

# Exibir resumo do modelo
summary(cox_model_interactions_only)


# Criar subconjuntos dos dados
df_low <- df %>% filter(steroid == "Steroid_Low")
df_high <- df %>% filter(steroid == "Steroid_High")

# Rodar modelo de Cox apenas para Steroid_Low
cox_model_low <- coxph(Surv(OS.Time, OS) ~ IGH + IGK + IGL + TRA + TRB + TRG + TRD, data = df_low)

# Rodar modelo de Cox apenas para Steroid_High
cox_model_high <- coxph(Surv(OS.Time, OS) ~ IGH + IGK + IGL + TRA + TRB + TRG + TRD, data = df_high)

# Exibir resultados
summary(cox_model_low)
summary(cox_model_high)

library(ggplot2)
library(dplyr)

# Extrair resultados do modelo de Cox
cox_summary <- summary(cox_model_interactions_only)  # Ajuste para seu modelo
cox_results <- data.frame(
  Variável = rownames(cox_summary$coefficients),
  HR = exp(cox_summary$coefficients[, 1]),  # Hazard Ratio
  lower_CI = exp(cox_summary$conf.int[, 3]),  # Limite inferior IC 95%
  upper_CI = exp(cox_summary$conf.int[, 4]),  # Limite superior IC 95%
  p_valor = cox_summary$coefficients[, 5]  # Valor p
) %>%
  arrange(desc(HR))  # Organizar variáveis pelo efeito

# Criar gráfico de floresta
ggplot(cox_results, aes(x = reorder(Variável, HR), y = HR)) +
  geom_point(aes(color = p_valor < 0.05), size = 3) +  # Pontos coloridos para destaque
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +  # Barras de erro IC 95%
  geom_text(aes(label = paste0("p=", signif(p_valor, 3))), hjust = -0.5, size = 3) +  # Adicionar p-valor
  coord_flip(expand = FALSE) +  # Ajustar posição das variáveis
  scale_y_log10(limits = c(0.5, 3)) +  # Escala log para melhorar visualização
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +  # Reduzir tamanho do texto
  labs(title = "Hazard Ratio das Interações - Modelo de Cox",
       x = "Variáveis",
       y = "Hazard Ratio (HR) (escala log)") +
  scale_color_manual(values = c("red", "blue"), guide = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")  # Linha referência HR=1




library(survival)
library(survminer)
library(dplyr)

# Criar colunas de alta/baixa expressão com base na mediana
df <- df %>%
  mutate(IGH_group = ifelse(IGH > median(IGH), "Alta", "Baixa"),
         IGK_group = ifelse(IGK > median(IGK), "Alta", "Baixa"),
         IGL_group = ifelse(IGL > median(IGL), "Alta", "Baixa"),
         TRA_group = ifelse(TRA > median(TRA), "Alta", "Baixa"),
         TRB_group = ifelse(TRB > median(TRB), "Alta", "Baixa"),
         TRG_group = ifelse(TRG > median(TRG), "Alta", "Baixa"),
         TRD_group = ifelse(TRD > median(TRD), "Alta", "Baixa"))

# Verificar distribuição dos grupos
table(df$IGH_group)


# Filtrar pacientes no grupo Steroid_High e Steroid_Low
df_high <- df %>% filter(steroid == "Steroid_High")
df_low <- df %>% filter(steroid == "Steroid_Low")

# Criar objetos de sobrevivência para cada grupo
surv_high <- Surv(df_high$OS.Time, df_high$OS)
surv_low <- Surv(df_low$OS.Time, df_low$OS)

# Ajustar modelos Kaplan-Meier
km_fit_high <- survfit(surv_high ~ TRG_group, data = df_high)
km_fit_low <- survfit(surv_low ~ TRG_group, data = df_low)

# Gráficos Kaplan-Meier comparando IGH nos grupos de steroid
plot_high <- ggsurvplot(km_fit_high, data = df_high, pval = TRUE, conf.int = TRUE,
                        title = "Kaplan-Meier IGH - Steroid High",
                        xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                        risk.table = TRUE, legend.title = "Grupo IGH", 
                        legend.labs = c("Baixa Expressão", "Alta Expressão"),
                        palette = c("#E69F00", "#56B4E9"))

plot_low <- ggsurvplot(km_fit_low, data = df_low, pval = TRUE, conf.int = TRUE,
                       title = "Kaplan-Meier IGH - Steroid Low",
                       xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                       risk.table = TRUE, legend.title = "Grupo IGH", 
                       legend.labs = c("Baixa Expressão", "Alta Expressão"),
                       palette = c("#E69F00", "#56B4E9"))

# Exibir gráficos
print(plot_high)
print(plot_low)




# Criar coluna de expressão alta/baixa baseada na mediana
df <- df %>%
  mutate(IGH_group = ifelse(IGH > median(IGH), "Alta", "Baixa"),
         steroid_IGH = paste(steroid, IGH_group, sep = "_"))#

# Gerar curvas Kaplan-Meier para interação entre steroid e IGH
km_fit_interaction <- survfit(Surv(OS.Time, OS) ~ steroid_IGH, data = df)

# Plotar Kaplan-Meier para interações
ggsurvplot(km_fit_interaction, data = df, pval = TRUE, conf.int = TRUE,
           title = "Kaplan-Meier - Interação IGH e Steroid na Sobrevida",
           xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
           risk.table = TRUE, legend.title = "Grupo",
           legend.labs = c("Steroid High - IGH Baixa", "Steroid High - IGH Alta",
                           "Steroid Low - IGH Baixa", "Steroid Low - IGH Alta"),
           palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))


# Lista de cadeias BCR/TCR
cadeias <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")

# Criar expressões alta/baixa para todas as cadeias
df <- df %>%
  mutate(across(all_of(cadeias), ~ ifelse(. > median(., na.rm = TRUE), "Alta", "Baixa"), 
                .names = "{col}_group"))


# Gerar e salvar gráficos
pdf("Kaplan_Meier_Cadeias_Steroid.pdf", width = 8, height = 6)  # Abrir arquivo PDF

for (cadeia in cadeias) {
  df <- df %>%
    mutate(Grupo = paste(steroid, get(paste0(cadeia, "_group")), sep = "_"))
  
  km_fit <- survfit(Surv(OS.Time, OS) ~ Grupo, data = df)
  
  plot_km <- ggsurvplot(km_fit, data = df, pval = TRUE, conf.int = TRUE,
                        title = paste("Kaplan-Meier - Interação", cadeia, "e Steroid"),
                        xlab = "Tempo (dias)", ylab = "Probabilidade de Sobrevida",
                        risk.table = TRUE, legend.title = "Grupo",
                        legend.labs = unique(df$Grupo),
                        palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"))
  
  print(plot_km)  # Salvar gráfico no PDF
}

dev.off() 
