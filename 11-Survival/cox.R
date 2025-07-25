library(survival)
library(survminer)
library(dplyr)

# Carregar os dados (assumindo que esta parte já está funcionando)
data_clones <- read.csv("../../6-dataExploration/df_data_clones.csv")
load("../../6-dataExploration/tcgaACC_pre_processed.RData")

arquivos <- list.files("../../7-diversity/metricsTrust4/results_pipeline_report/", 
                       pattern = "\\.tsv$", full.names = TRUE)

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados_temp$Abundance))
  rownames(dados_abundance) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_abundance)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))

# Aplicar a transformação logarítmica nos dados filtrados
dados_transformed_filtrados <- dados_filtrados %>%
  mutate(across(c(IGH, IGK, IGL, TRA, TRB, TRG, TRD), ~log2(.x + 1)))

# Agora rode a análise de Cox com os dados filtrados
fit_multivariate_filtrado <- coxph(Surv(OS.Time, OS) ~ TRA * TRB + steroid + 
                                     tumor_stage, data = dados_transformed_filtrados)
summary(fit_multivariate_filtrado)

cox.zph(fit_multivariate_filtrado) # Verificação da proporcionalidade dos riscos
ggforest(fit_multivariate_filtrado, data = dados_transformed_filtrados) # Forest plot


#------------------------------------------------------------------------------#
# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_entropy <- data.frame(t(dados_temp$Entropy))
  rownames(dados_entropy) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_entropy)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))

# Aplicar a transformação logarítmica nos dados filtrados
#dados_transformed_filtrados <- dados_filtrados %>%
#  mutate(across(c(IGH, IGK, IGL, TRA, TRB, TRG, TRD), ~log2(.x + 1)))

# Agora rode a análise de Cox com os dados filtrados
fit_multivariate_filtrado <- coxph(Surv(OS.Time, OS) ~ TRA * TRB + steroid + 
                                     tumor_stage, data = dados_filtrados)
summary(fit_multivariate_filtrado)

cox.zph(fit_multivariate_filtrado) # Verificação da proporcionalidade dos riscos
ggforest(fit_multivariate_filtrado, data = dados_filtrados) # Forest plot

#------------------------------------------------------------------------------#

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados_temp$Abundance))
  rownames(dados_abundance) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_abundance)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))

# Aplicar a transformação logarítmica nos dados filtrados
dados_transformed_filtrados <- dados_filtrados %>%
  mutate(across(c(IGH, IGK, IGL, TRA, TRB, TRG, TRD), ~log2(.x + 1)))

# --- Criar subgrupos de dados ---
dados_steroid_high <- dados_transformed_filtrados %>%
  filter(steroid == "Steroid_High")

dados_steroid_low <- dados_transformed_filtrados %>%
  filter(steroid == "Steroid_Low")

cat("Número de pacientes Steroid_High:", nrow(dados_steroid_high), "\n")
cat("Número de pacientes Steroid_Low:", nrow(dados_steroid_low), "\n")

# Verifique o número de eventos em cada grupo (crucial para o poder estatístico)
cat("Eventos em Steroid_High:", sum(dados_steroid_high$OS), "\n")
cat("Eventos em Steroid_Low:", sum(dados_steroid_low$OS), "\n")

# --- Rodar modelos para cada grupo ---

# Exemplo para Abundância de TRB no grupo Steroid_High
fit_high_abundance <- coxph(Surv(OS.Time, OS) ~ TRB + tumor_stage, data = dados_steroid_high)
summary(fit_high_abundance)
cox.zph(fit_high_abundance)
ggforest(fit_high_abundance, data = dados_steroid_high)

#------------------------------------------------------------------------------#

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_entropy <- data.frame(t(dados_temp$Entropy))
  rownames(dados_entropy) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_entropy)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))


# --- Criar subgrupos de dados ---
dados_steroid_high <- dados_filtrados %>%
  filter(steroid == "Steroid_High")

dados_steroid_low <- dados_filtrados %>%
  filter(steroid == "Steroid_Low")

cat("Número de pacientes Steroid_High:", nrow(dados_steroid_high), "\n")
cat("Número de pacientes Steroid_Low:", nrow(dados_steroid_low), "\n")

# Verifique o número de eventos em cada grupo (crucial para o poder estatístico)
cat("Eventos em Steroid_High:", sum(dados_steroid_high$OS), "\n")
cat("Eventos em Steroid_Low:", sum(dados_steroid_low$OS), "\n")

# --- Rodar modelos para cada grupo ---

# Exemplo para Abundância de TRB no grupo Steroid_High
fit_high_abundance <- coxph(Surv(OS.Time, OS) ~ TRB + tumor_stage, data = dados_steroid_high)
summary(fit_high_abundance)
cox.zph(fit_high_abundance)
ggforest(fit_high_abundance, data = dados_steroid_high)

#------------------------------------------------------------------------------#

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados_temp$Abundance))
  rownames(dados_abundance) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_abundance)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))

# Aplicar a transformação logarítmica nos dados filtrados
dados_transformed_filtrados <- dados_filtrados %>%
  mutate(across(c(IGH, IGK, IGL, TRA, TRB, TRG, TRD), ~log2(.x + 1)))

# --- Criar subgrupos de dados ---
dados_steroid_high <- dados_transformed_filtrados %>%
  filter(steroid == "Steroid_High")

dados_steroid_low <- dados_transformed_filtrados %>%
  filter(steroid == "Steroid_Low")

cat("Número de pacientes Steroid_High:", nrow(dados_steroid_high), "\n")
cat("Número de pacientes Steroid_Low:", nrow(dados_steroid_low), "\n")

# Verifique o número de eventos em cada grupo (crucial para o poder estatístico)
cat("Eventos em Steroid_High:", sum(dados_steroid_high$OS), "\n")
cat("Eventos em Steroid_Low:", sum(dados_steroid_low$OS), "\n")

# --- Rodar modelos para cada grupo ---

# Exemplo para Abundância de TRB no grupo Steroid_High
fit_high_abundance <- coxph(Surv(OS.Time, OS) ~ TRA * TRB + tumor_stage, data = dados_steroid_high)
summary(fit_high_abundance)
cox.zph(fit_high_abundance)
ggforest(fit_high_abundance, data = dados_steroid_high)

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_entropy <- data.frame(t(dados_temp$Entropy))
  rownames(dados_entropy) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_entropy)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))


# --- Criar subgrupos de dados ---
dados_steroid_high <- dados_filtrados %>%
  filter(steroid == "Steroid_High")

dados_steroid_low <- dados_filtrados %>%
  filter(steroid == "Steroid_Low")

cat("Número de pacientes Steroid_High:", nrow(dados_steroid_high), "\n")
cat("Número de pacientes Steroid_Low:", nrow(dados_steroid_low), "\n")

# Verifique o número de eventos em cada grupo (crucial para o poder estatístico)
cat("Eventos em Steroid_High:", sum(dados_steroid_high$OS), "\n")
cat("Eventos em Steroid_Low:", sum(dados_steroid_low$OS), "\n")

# --- Rodar modelos para cada grupo ---

# Exemplo para Abundância de TRB no grupo Steroid_High
fit_high_abundance <- coxph(Surv(OS.Time, OS) ~ TRB + tumor_stage, data = dados_steroid_high)
summary(fit_high_abundance)
cox.zph(fit_high_abundance)
ggforest(fit_high_abundance, data = dados_steroid_high)

#------------------------------------------------------------------------------#

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_entropy <- data.frame(t(dados_temp$Entropy))
  rownames(dados_entropy) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_entropy)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))


# --- Criar subgrupos de dados ---
dados_steroid_high <- dados_filtrados %>%
  filter(steroid == "Steroid_High")

dados_steroid_low <- dados_filtrados %>%
  filter(steroid == "Steroid_Low")

cat("Número de pacientes Steroid_High:", nrow(dados_steroid_high), "\n")
cat("Número de pacientes Steroid_Low:", nrow(dados_steroid_low), "\n")

# Verifique o número de eventos em cada grupo (crucial para o poder estatístico)
cat("Eventos em Steroid_High:", sum(dados_steroid_high$OS), "\n")
cat("Eventos em Steroid_Low:", sum(dados_steroid_low$OS), "\n")

# --- Rodar modelos para cada grupo ---

# Exemplo para Abundância de TRB no grupo Steroid_High
fit_high_abundance <- coxph(Surv(OS.Time, OS) ~ TRA * TRB + tumor_stage, data = dados_steroid_high)
summary(fit_high_abundance)
cox.zph(fit_high_abundance)
ggforest(fit_high_abundance, data = dados_steroid_high)

#------------------------------------------------------------------------------#
library(survival)
library(survminer)
library(dplyr)

# Carregar os dados (assumindo que esta parte já está funcionando)
data_clones <- read.csv("../../6-dataExploration/df_data_clones.csv")
load("../../6-dataExploration/tcgaACC_pre_processed.RData")

arquivos <- list.files("../../7-diversity/metricsTrust4/results_pipeline_report/", 
                       pattern = "\\.tsv$", full.names = TRUE)

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados_temp$Abundance))
  rownames(dados_abundance) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_abundance)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))

# Aplicar a transformação logarítmica nos dados filtrados
dados_transformed_filtrados <- dados_filtrados %>%
  mutate(across(c(IGH, IGK, IGL, TRA, TRB, TRG, TRD), ~log2(.x + 1)))


# Assumindo que dados_transformed_filtrados está carregado e tem as variáveis OS.Time, OS, steroid, tumor_stage

# 1. Curva de Kaplan-Meier para 'steroid'
fit_km_steroid <- surv_fit(Surv(OS.Time, OS) ~ steroid, data = dados_transformed_filtrados)
ggsurvplot(fit_km_steroid,
           data = dados_transformed_filtrados,
           pval = TRUE,             # Mostrar p-valor do teste Log-Rank
           conf.int = TRUE,         # Mostrar intervalo de confiança
           risk.table = TRUE,       # Mostrar tabela de risco (número de pacientes em risco ao longo do tempo)
           risk.table.col = "strata",
           surv.median.line = "hv", # Mostrar linha da mediana de sobrevida
           ggtheme = theme_bw(),    # Tema do gráfico
           palette = c("#E7B800", "#2E9FDF"), # Cores para os grupos (ex: High vs Low)
           title = "Curva de Kaplan-Meier por Status de Esteróide",
           xlab = "Tempo (meses)",
           legend.title = "Status de Esteróide"
)

# 2. Curva de Kaplan-Meier para 'tumor_stage'
fit_km_tumor_stage <- surv_fit(Surv(OS.Time, OS) ~ tumor_stage, data = dados_transformed_filtrados)
ggsurvplot(fit_km_tumor_stage,
           data = dados_transformed_filtrados,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           palette = "Set1", # Uma paleta de cores para múltiplos grupos
           title = "Curva de Kaplan-Meier por Estágio do Tumor",
           xlab = "Tempo (meses)",
           legend.title = "Estágio do Tumor"
)

# 3. Curva de Kaplan-Meier para uma cadeia (ex: TRA) - se quiser confirmar visualmente a falta de efeito
# Para fazer isso, você precisaria categorizar a variável contínua (TRA) em grupos (ex: acima/abaixo da mediana, ou tercis/quartis)
# Exemplo categorizando TRA pela mediana
 dados_transformed_filtrados$TRA_group <- ifelse(dados_transformed_filtrados$TRA > median(dados_transformed_filtrados$TRA), "High_TRA", "Low_TRA")
 fit_km_TRA <- surv_fit(Surv(OS.Time, OS) ~ TRA_group + steroid, data = dados_transformed_filtrados)
 ggsurvplot(fit_km_TRA,
            data = dados_transformed_filtrados,
            pval = TRUE,
            conf.int = TRUE,
            risk.table = TRUE,
            risk.table.col = "strata",
            surv.median.line = "hv",
            ggtheme = theme_bw(),
            title = "Curva de Kaplan-Meier por Abundância de TRA (Mediana)",
            xlab = "Tempo (meses)",
            legend.title = "Grupo de TRA"
 )
 

 
 dados_transformed_filtrados$TRB_group <- ifelse(dados_transformed_filtrados$TRB > median(dados_transformed_filtrados$TRB), "High_TRB", "Low_TRB")
 fit_km_TRB <- surv_fit(Surv(OS.Time, OS) ~ TRB_group + steroid, data = dados_transformed_filtrados)
 ggsurvplot(fit_km_TRB,
            data = dados_transformed_filtrados,
            pval = TRUE,
            conf.int = TRUE,
            risk.table = TRUE,
            risk.table.col = "strata",
            surv.median.line = "hv",
            ggtheme = theme_bw(),
            title = "Curva de Kaplan-Meier por Abundância de TRB (Mediana)",
            xlab = "Tempo (meses)",
            legend.title = "Grupo de TRB"
 )
 
 dados_transformed_filtrados$IGH_group <- ifelse(dados_transformed_filtrados$IGH > median(dados_transformed_filtrados$IGH), "High_IGH", "Low_IGH")
 fit_km_IGH <- surv_fit(Surv(OS.Time, OS) ~ IGH_group + steroid, data = dados_transformed_filtrados)
 ggsurvplot(fit_km_IGH,
            data = dados_transformed_filtrados,
            pval = TRUE,
            conf.int = TRUE,
            risk.table = TRUE,
            risk.table.col = "strata",
            surv.median.line = "hv",
            ggtheme = theme_bw(),
            title = "Curva de Kaplan-Meier por Abundância de IGH (Mediana)",
            xlab = "Tempo (meses)",
            legend.title = "Grupo de IGH"
 )

 dados_transformed_filtrados$IGK_group <- ifelse(dados_transformed_filtrados$IGK > median(dados_transformed_filtrados$IGK), "High_IGK", "Low_IGK")
 fit_km_IGK <- surv_fit(Surv(OS.Time, OS) ~ IGK_group + steroid, data = dados_transformed_filtrados)
 ggsurvplot(fit_km_IGK,
            data = dados_transformed_filtrados,
            pval = TRUE,
            conf.int = TRUE,
            risk.table = TRUE,
            risk.table.col = "strata",
            surv.median.line = "hv",
            ggtheme = theme_bw(),
            title = "Curva de Kaplan-Meier por Abundância de IGK (Mediana)",
            xlab = "Tempo (meses)",
            legend.title = "Grupo de IGK"
 )
 
 dados_transformed_filtrados$IGL_group <- ifelse(dados_transformed_filtrados$IGL > median(dados_transformed_filtrados$IGL), "High_IGL", "Low_IGL")
 fit_km_IGL <- surv_fit(Surv(OS.Time, OS) ~ IGL_group + steroid, data = dados_transformed_filtrados)
 ggsurvplot(fit_km_IGL,
            data = dados_transformed_filtrados,
            pval = TRUE,
            conf.int = TRUE,
            risk.table = TRUE,
            risk.table.col = "strata",
            surv.median.line = "hv",
            ggtheme = theme_bw(),
            title = "Curva de Kaplan-Meier por Abundância de IGL (Mediana)",
            xlab = "Tempo (meses)",
            legend.title = "Grupo de IGL"
 )
#------------------------------------------------------------------------------#
 library(survival)
 library(survminer)
 library(dplyr)
 
 # Carregar os dados (assumindo que esta parte já está funcionando)
 data_clones <- read.csv("../../6-dataExploration/df_data_clones.csv")
 load("../../6-dataExploration/tcgaACC_pre_processed.RData")
 
 arquivos <- list.files("../../7-diversity/metricsTrust4/results_pipeline_report/", 
                        pattern = "\\.tsv$", full.names = TRUE)
 
 
 # montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados_temp <- read.table(arquivo, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  dados_entropy <- data.frame(t(dados_temp$Entropy))
  rownames(dados_entropy) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_entropy)
}

colnames(data_metrics) <- dados_temp$chain

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

# **NOVA ETAPA:** Filtrar a categoria "not reported" de tumor_stage
dados_filtrados <- dados %>%
  filter(tumor_stage != "not reported") %>%
  # É uma boa prática re-dropar os níveis do fator após a filtragem
  mutate(tumor_stage = droplevels(tumor_stage))

dados_transformed_filtrados <- dados_filtrados
 
 

dados_transformed_filtrados$TRA_group <- ifelse(dados_transformed_filtrados$TRA > median(dados_transformed_filtrados$TRA), "High_TRA", "Low_TRA")
fit_km_TRA <- surv_fit(Surv(OS.Time, OS) ~ TRA_group + steroid, data = dados_transformed_filtrados)
ggsurvplot(fit_km_TRA,
           data = dados_transformed_filtrados,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           title = "Curva de Kaplan-Meier por Abundância de TRA (Mediana)",
           xlab = "Tempo (meses)",
           legend.title = "Grupo de TRA"
)



dados_transformed_filtrados$TRB_group <- ifelse(dados_transformed_filtrados$TRB > median(dados_transformed_filtrados$TRB), "High_TRB", "Low_TRB")
fit_km_TRB <- surv_fit(Surv(OS.Time, OS) ~ TRB_group + steroid, data = dados_transformed_filtrados)
ggsurvplot(fit_km_TRB,
           data = dados_transformed_filtrados,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           title = "Curva de Kaplan-Meier por Abundância de TRB (Mediana)",
           xlab = "Tempo (meses)",
           legend.title = "Grupo de TRB"
)

dados_transformed_filtrados$IGH_group <- ifelse(dados_transformed_filtrados$IGH > median(dados_transformed_filtrados$IGH), "High_IGH", "Low_IGH")
fit_km_IGH <- surv_fit(Surv(OS.Time, OS) ~ IGH_group + steroid, data = dados_transformed_filtrados)
ggsurvplot(fit_km_IGH,
           data = dados_transformed_filtrados,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           title = "Curva de Kaplan-Meier por Abundância de IGH (Mediana)",
           xlab = "Tempo (meses)",
           legend.title = "Grupo de IGH"
)

dados_transformed_filtrados$IGK_group <- ifelse(dados_transformed_filtrados$IGK > median(dados_transformed_filtrados$IGK), "High_IGK", "Low_IGK")
fit_km_IGK <- surv_fit(Surv(OS.Time, OS) ~ IGK_group + steroid, data = dados_transformed_filtrados)
ggsurvplot(fit_km_IGK,
           data = dados_transformed_filtrados,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           title = "Curva de Kaplan-Meier por Abundância de IGK (Mediana)",
           xlab = "Tempo (meses)",
           legend.title = "Grupo de IGK"
)

dados_transformed_filtrados$IGL_group <- ifelse(dados_transformed_filtrados$IGL > median(dados_transformed_filtrados$IGL), "High_IGL", "Low_IGL")
fit_km_IGL <- surv_fit(Surv(OS.Time, OS) ~ IGL_group + steroid, data = dados_transformed_filtrados)
ggsurvplot(fit_km_IGL,
           data = dados_transformed_filtrados,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           title = "Curva de Kaplan-Meier por Abundância de IGL (Mediana)",
           xlab = "Tempo (meses)",
           legend.title = "Grupo de IGL"
)
 
 
 
 
 
#------------------------------------------------------------------------------#
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 library(survival)
 library(survminer)
 library(dplyr)
 
 # Assumindo que dados_transformed_filtrados está carregado
 
 # 1. Filtrar os dados para incluir apenas pacientes no estágio II
 dados_stage_ii <- dados_transformed_filtrados %>%
   filter(tumor_stage == "stage ii")
 
 # Verificar o número de pacientes e eventos neste subgrupo
 # print(table(dados_stage_ii$steroid))
 # print(sum(dados_stage_ii$OS))
 
 # 2. Criar a curva de Kaplan-Meier para 'steroid' usando APENAS os dados do estágio II
 fit_km_steroid_stage_ii <- surv_fit(Surv(OS.Time, OS) ~ steroid, data = dados_stage_ii)
 
 ggsurvplot(fit_km_steroid_stage_ii,
            data = dados_stage_ii,
            pval = TRUE,             # Mostrar p-valor do teste Log-Rank
            conf.int = TRUE,         # Mostrar intervalo de confiança
            risk.table = TRUE,       # Mostrar tabela de risco
            risk.table.col = "strata",
            linetype = "strata",
            surv.median.line = "hv",
            ggtheme = theme_bw(),
            palette = c("#E7B800", "#2E9FDF"), # Cores para os grupos
            title = "Curva de Kaplan-Meier por Status de Esteróide (Apenas Estágio II)",
            xlab = "Tempo (meses)",
            legend.title = "Status de Esteróide"
 )
 
 # Opcional: Você também pode fazer isso para as abundâncias/entropias das cadeias dentro do Estágio II
 # Por exemplo, TRA dentro do Estágio II:
  dados_stage_ii$TRA_group <- ifelse(dados_stage_ii$TRA > median(dados_stage_ii$TRA), "High_TRA", "Low_TRA")
  fit_km_TRA_stage_ii <- surv_fit(Surv(OS.Time, OS) ~ TRA_group + steroid, data = dados_stage_ii)
  ggsurvplot(fit_km_TRA_stage_ii,
             data = dados_stage_ii,
             pval = TRUE,
             conf.int = TRUE,
             risk.table = TRUE,
             risk.table.col = "strata",
             linetype = "strata",
             surv.median.line = "hv",
             ggtheme = theme_bw(),
             title = "Curva de Kaplan-Meier por Abundância de TRA (Apenas Estágio II)",
             xlab = "Tempo (meses)",
             legend.title = "Grupo de TRA"
  )
  
  
  dados_stage_ii$TRB_group <- ifelse(dados_stage_ii$TRB > median(dados_stage_ii$TRB), "High_TRB", "Low_TRB")
  fit_km_TRB_stage_ii <- surv_fit(Surv(OS.Time, OS) ~ TRB_group + steroid, data = dados_stage_ii)
  ggsurvplot(fit_km_TRB_stage_ii,
             data = dados_stage_ii,
             pval = TRUE,
             conf.int = TRUE,
             risk.table = TRUE,
             risk.table.col = "strata",
             linetype = "strata",
             surv.median.line = "hv",
             ggtheme = theme_bw(),
             title = "Curva de Kaplan-Meier por Abundância de TRB (Apenas Estágio II)",
             xlab = "Tempo (meses)",
             legend.title = "Grupo de TRB"
  )
  
  dados_stage_ii$IGH_group <- ifelse(dados_stage_ii$IGH > median(dados_stage_ii$IGH), "High_IGH", "Low_IGH")
  fit_km_IGH_stage_ii <- surv_fit(Surv(OS.Time, OS) ~ IGH_group + steroid, data = dados_stage_ii)
  ggsurvplot(fit_km_IGH_stage_ii,
             data = dados_stage_ii,
             pval = TRUE,
             conf.int = TRUE,
             risk.table = TRUE,
             risk.table.col = "strata",
             linetype = "strata",
             surv.median.line = "hv",
             ggtheme = theme_bw(),
             title = "Curva de Kaplan-Meier por Abundância de IGH (Apenas Estágio II)",
             xlab = "Tempo (meses)",
             legend.title = "Grupo de IGH"
  )
  
  
  dados_stage_ii$IGK_group <- ifelse(dados_stage_ii$IGK > median(dados_stage_ii$IGK), "High_IGK", "Low_IGK")
  fit_km_IGK_stage_ii <- surv_fit(Surv(OS.Time, OS) ~ IGK_group + steroid, data = dados_stage_ii)
  ggsurvplot(fit_km_IGK_stage_ii,
             data = dados_stage_ii,
             pval = TRUE,
             conf.int = TRUE,
             risk.table = TRUE,
             risk.table.col = "strata",
             linetype = "strata",
             surv.median.line = "hv",
             ggtheme = theme_bw(),
             title = "Curva de Kaplan-Meier por Abundância de IGK (Apenas Estágio II)",
             xlab = "Tempo (meses)",
             legend.title = "Grupo de IGK"
  )
  
  
  dados_stage_ii$IGL_group <- ifelse(dados_stage_ii$IGL > median(dados_stage_ii$IGL), "High_IGL", "Low_IGL")
  fit_km_IGL_stage_ii <- surv_fit(Surv(OS.Time, OS) ~ IGL_group + steroid, data = dados_stage_ii)
  ggsurvplot(fit_km_IGL_stage_ii,
             data = dados_stage_ii,
             pval = TRUE,
             conf.int = TRUE,
             risk.table = TRUE,
             risk.table.col = "strata",
             linetype = "strata",
             surv.median.line = "hv",
             ggtheme = theme_bw(),
             title = "Curva de Kaplan-Meier por Abundância de IGL (Apenas Estágio II)",
             xlab = "Tempo (meses)",
             legend.title = "Grupo de IGL"
  )
  