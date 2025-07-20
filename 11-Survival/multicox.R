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

dados <- df

# Converter variáveis para o tipo correto
dados$steroid <- as.factor(dados$steroid)
dados$tumor_stage <- as.factor(dados$tumor_stage)

dados_transformed <- dados %>%
  mutate(across(c(IGH, IGK, IGL, TRA, TRB, TRG, TRD), ~log2(.x + 1)))

# Análise de Cox Multivariada
fit_multivariate <- coxph(Surv(OS.Time, OS) ~ steroid + IGH, data = dados)
summary(fit_multivariate)





M_igh_igk_steroid <- coxph(Surv(OS.Time, OS) ~ IGH + IGK +steroid, data = dados)
summary(M_igh_igk_steroid)



# Lista de variáveis de repertório para testar
repertoire_vars_to_test <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")

# Armazenar os resultados
results_list <- list()

for(var in repertoire_vars_to_test){
  # Modelo univariado para cada variável
  fmla <- as.formula(paste("Surv(OS.Time, OS) ~ ", var, "+ steroid"))
  fit <- coxph(fmla, data = dados)
  results_list[[var]] <- summary(fit)
}

# Extrair p-valores para a variável de repertório de cada modelo
p_values <- sapply(results_list, function(x) x$coefficients[var, "Pr(>|z|)"])

# Aplicar a correção de Benjamini-Hochberg
p_values_adjusted <- p.adjust(p_values, method = "BH")

# Mostrar os resultados corrigidos
print(p_values_adjusted)







# 4. Correção para Múltiplas Comparações
# Como você faria múltiplos testes de Cox, você deve ajustar os p-valores.

# Lista de variáveis de repertório para testar
repertoire_vars_to_test <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")

# Crie um dataframe para armazenar os resultados
results_df <- data.frame(
  variable = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop para rodar os modelos e extrair os p-valores
for(var in repertoire_vars_to_test){
  # Crie a fórmula dinamicamente
  fmla <- as.formula(paste("Surv(OS.Time, OS) ~ Steroid_Phenotype + ", var))
  
  # Use tryCatch para lidar com modelos que podem falhar (e.g., com dados esparsos)
  fit_summary <- tryCatch({
    fit <- coxph(fmla, data = dados_transformed) # Usando os dados transformados
    summary(fit)
  }, error = function(e) {
    # Em caso de erro (ex: variância zero), retorne NA
    return(list(p_value = NA))
  })
  
  # Verifique se a análise foi bem-sucedida e se o p-valor existe
  if (!is.null(fit_summary$p_value)) {
    # Se fit_summary$p_value não for nulo, significa que houve um erro e retornamos NA
    p_val_var <- fit_summary$p_value
  } else {
    # Extrair o p-valor para a variável 'var' do sumário do modelo
    # o nome da linha do coeficiente da variável é 'var'
    p_val_var <- fit_summary$coefficients[var, "Pr(>|z|)"]
  }
  
  # Adicione o resultado ao dataframe
  results_df <- rbind(results_df, data.frame(variable = var, p_value = p_val_var))
}

# Remover linhas com p-valores NA se houver
results_df <- na.omit(results_df)

# Aplicar a correção de Benjamini-Hochberg (FDR)
results_df$p_adjusted <- p.adjust(results_df$p_value, method = "BH")

# Mostrar os resultados corrigidos
print(results_df)

