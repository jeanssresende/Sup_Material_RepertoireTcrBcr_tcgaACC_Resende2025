library(survival)
library(survminer)
library(dplyr)

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")

load("../6-dataExploration/tcgaACC_pre_processed.RData")

arquivos <- list.files("../7-diversity/metricsTrust4/results_pipeline_report/", 
                       pattern = "\\.tsv$", full.names = TRUE)

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados$Abundance))
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
data_metrics$tumor_stage <- data_clones$tumor_stage


data <- data_metrics

rownames(data) <- data_clones$TCGA_barcode
head(data)

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

df <- df %>%
  filter(steroid == "Steroid_High")

# Lista de cadeias BCR/TCR
cadeias <- c("IGH", "IGK", "IGL", "TRA", "TRB")

# Criar expressões alta/baixa para todas as cadeias
df <- df %>%
  mutate(across(all_of(cadeias), ~ ifelse(. > median(., na.rm = TRUE), "Alta", "Baixa"), 
                .names = "{col}_group"))


# Gerar e salvar gráficos
pdf("sobrevida/Kaplan_Meier_Cadeias_Steroid_abundance_HSP.pdf", width = 8, height = 6)  # Abrir arquivo PDF

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




################################################################################

library(survival)
library(survminer)
library(dplyr)

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")

load("../6-dataExploration/tcgaACC_pre_processed.RData")

arquivos <- list.files("../7-diversity/metricsTrust4/results_pipeline_report/", 
                       pattern = "\\.tsv$", full.names = TRUE)

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados$Abundance))
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
data_metrics$tumor_stage <- data_clones$tumor_stage


data <- data_metrics

rownames(data) <- data_clones$TCGA_barcode
head(data)

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

df <- df %>%
  filter(steroid == "Steroid_Low")

# Lista de cadeias BCR/TCR
cadeias <- c("IGH", "IGK", "IGL", "TRA", "TRB")

# Criar expressões alta/baixa para todas as cadeias
df <- df %>%
  mutate(across(all_of(cadeias), ~ ifelse(. > median(., na.rm = TRUE), "Alta", "Baixa"), 
                .names = "{col}_group"))


# Gerar e salvar gráficos
pdf("sobrevida/Kaplan_Meier_Cadeias_Steroid_abundance_LSP.pdf", width = 8, height = 6)  # Abrir arquivo PDF

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


################################################################################


data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")

load("../6-dataExploration/tcgaACC_pre_processed.RData")

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

df <- df %>%
  filter(steroid == "Steroid_High")

# Lista de cadeias BCR/TCR
cadeias <- c("IGH", "IGK", "IGL", "TRA", "TRB")

# Criar expressões alta/baixa para todas as cadeias
df <- df %>%
  mutate(across(all_of(cadeias), ~ ifelse(. > median(., na.rm = TRUE), "Alta", "Baixa"), 
                .names = "{col}_group"))


# Gerar e salvar gráficos
pdf("sobrevida/Kaplan_Meier_Cadeias_Steroid_entropy_HSP.pdf", width = 8, height = 6)  # Abrir arquivo PDF

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


################################################################################

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")

load("../6-dataExploration/tcgaACC_pre_processed.RData")

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

idx <- match(rownames(data), tcgaACC$barcode)
data$OS <- tcgaACC$OS[idx]
data$OS.Time <- tcgaACC$OS.Time[idx]

df <- data

df <- df %>%
  filter(steroid == "Steroid_Low")

# Lista de cadeias BCR/TCR
cadeias <- c("IGH", "IGK", "IGL", "TRA", "TRB")

# Criar expressões alta/baixa para todas as cadeias
df <- df %>%
  mutate(across(all_of(cadeias), ~ ifelse(. > median(., na.rm = TRUE), "Alta", "Baixa"), 
                .names = "{col}_group"))


# Gerar e salvar gráficos
pdf("sobrevida/Kaplan_Meier_Cadeias_Steroid_entropy_LSP.pdf", width = 8, height = 6)  # Abrir arquivo PDF

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

p.adjust(c(0.66,0.32,0.16,
           0.66,0.41,
           0.091,0.42,0.51,
           0.069,0.057,
           0.055,0.31,0.27,
           0.41,0.088,
           0.039,0.71,0.68,
           0.069,0.33),
         method = "BH")

#  [1] 0.7100 0.5500 0.4000 
#      0.7100 0.5600 
#      0.2600 0.5600 0.6375 
#      0.2600 0.2600 
#      0.2600 0.5500 0.5500 
#      0.5600 0.2600 
#      0.2600 0.7100 0.7100 
#0.2600 0.5500