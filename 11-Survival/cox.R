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
