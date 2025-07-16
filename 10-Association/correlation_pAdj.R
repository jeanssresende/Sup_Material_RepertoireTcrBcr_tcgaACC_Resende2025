dir_metrics <- "../7-diversity/metricsTrust4/results_pipeline_report/"
load("../6-dataExploration/tcgaACC_pre_processed.RData")
data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")
load("data_purity.RData")
load("../4-IDs-metadata/metadata_2025.RData")


arquivos <- list.files(dir_metrics, pattern = "\\.tsv$", full.names = TRUE)

# montagem do objeto
data_abundance <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados$Abundance))
  rownames(dados_abundance) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_abundance <- rbind(data_abundance, dados_abundance)
}

colnames(data_abundance) <- dados$chain

# adicionando os casos que o TRUST4 nao encontrou TCR e BCR
mat_add <- data.frame(IGH = rep(0,10),
                      IGK = rep(0,10), 
                      IGL = rep(0,10), 
                      TRA = rep(0,10), 
                      TRB = rep(0,10), 
                      TRG = rep(0,10), 
                      TRD = rep(0,10))

#rownames(data_clones) <- data_clones$sample_id

rownames(mat_add) <- data_clones$sample_id[data_clones$sample_id %in% 
                                             rownames(data_abundance) == FALSE]

data_abundance <- rbind(data_abundance, mat_add)


# montagem do objeto
data_entropy <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_entropy <- data.frame(t(dados$Entropy))
  rownames(dados_entropy) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_entropy <- rbind(data_entropy, dados_entropy)
}

colnames(data_entropy) <- dados$chain

# adicionando os casos que o TRUST4 nao encontrou TCR e BCR
mat_add <- data.frame(IGH = rep(0,10),
                      IGK = rep(0,10), 
                      IGL = rep(0,10), 
                      TRA = rep(0,10), 
                      TRB = rep(0,10), 
                      TRG = rep(0,10), 
                      TRD = rep(0,10))

#rownames(data_clones) <- data_clones$sample_id

rownames(mat_add) <- data_clones$sample_id[data_clones$sample_id %in% 
                                             rownames(data_entropy) == FALSE]

data_entropy <- rbind(data_entropy, mat_add)

# metadata
data <- data.frame(
  barcode = tcgaACC$barcode,
  Th1.Cells = tcgaACC$Th1.Cells,
  Th2.Cells = tcgaACC$Th2.Cells,
  Th17.Cells = tcgaACC$Th17.Cells,
#  Wound.Healing = tcgaACC$Wound.Healing,
  IFN.gamma.Response = tcgaACC$IFN.gamma.Response,
  TGF.beta.Response = tcgaACC$TGF.beta.Response,
  Macrophage.Regulation = tcgaACC$Macrophage.Regulation,
  Lymphocyte.Infiltration = tcgaACC$Lymphocyte.Infiltration.Signature.Score#,
#  Nonsilent.Mutation.Rate = tcgaACC$Nonsilent.Mutation.Rate,
#  Silent.Mutation.Rate = tcgaACC$Silent.Mutation.Rate,
#  SNV.Neoantigens = tcgaACC$SNV.Neoantigens
)

#idx <- match(data$barcode, data_purity$barcode)
#data$paper_purity <- as.numeric(data_purity$paper_purity[idx]) 

# Abundancia
idx <- match(rownames(data_abundance), metadata$file_id)
rownames(data_abundance) <- metadata$TCGA_barcode[idx]

idx <- match(data$barcode, rownames(data_abundance))
data$IGH_Abundance <- data_abundance$IGH[idx]
data$IGK_Abundance <- data_abundance$IGK[idx]
data$IGL_Abundance <- data_abundance$IGL[idx]
data$TRA_Abundance <- data_abundance$TRA[idx]
data$TRB_Abundance <- data_abundance$TRB[idx]
data$TRG_Abundance <- data_abundance$TRG[idx]
data$TRD_Abundance <- data_abundance$TRD[idx]

# Entropia
idx <- match(rownames(data_entropy), metadata$file_id)
rownames(data_entropy) <- metadata$TCGA_barcode[idx]

idx <- match(data$barcode, rownames(data_entropy))
data$IGH_Entropy <- data_entropy$IGH[idx]
data$IGK_Entropy <- data_entropy$IGK[idx]
data$IGL_Entropy <- data_entropy$IGL[idx]
data$TRA_Entropy <- data_entropy$TRA[idx]
data$TRB_Entropy <- data_entropy$TRB[idx]
data$TRG_Entropy <- data_entropy$TRG[idx]
data$TRD_Entropy <- data_entropy$TRD[idx]

#barcode_bcrtcr <- coldataACC$barcode[coldataACC$sample_id %in%
#                                       gsub("_report", "", abundance$Sample)]

# # Abundance
# idx <- match(gsub("_report","",abundance$Sample), coldataACC$sample_id)
# abundance$barcode <- coldataACC$barcode[idx]
# 
# data <- data[data$barcode %in% barcode_bcrtcr,]
# 
# idx <- match(data$barcode, abundance$barcode)
# 
# data$IGH_Abundance <- abundance$IGH[idx]
# data$IGK_Abundance <- abundance$IGK[idx]
# data$IGL_Abundance <- abundance$IGL[idx]
# data$TRA_Abundance <- abundance$TRA[idx]
# data$TRB_Abundance <- abundance$TRB[idx]
# 
# # Entropy
# idx <- match(gsub("_report","",entropy$Sample), coldataACC$sample_id)
# entropy$barcode <- coldataACC$barcode[idx]
# 
# idx <- match(data$barcode, entropy$barcode)
# 
# data$IGH_Entropy <- entropy$IGH[idx]
# data$IGK_Entropy <- entropy$IGK[idx]
# data$IGL_Entropy <- entropy$IGL[idx]
# data$TRA_Entropy <- entropy$TRA[idx]
# data$TRB_Entropy <- entropy$TRB[idx]



shapiro.test(data$Th1.Cells) # normal
shapiro.test(data$Th2.Cells) # normal
shapiro.test(data$Th17.Cells) # nao normal
shapiro.test(data$Wound.Healing) # normal
shapiro.test(data$IFN.gamma.Response) # normal
shapiro.test(data$TGF.beta.Response) # nao normal
shapiro.test(data$Macrophage.Regulation) # normal
#...

library(ggplot2)
library(reshape2)  # Para transformar a matriz de correlação
library(ggpubr)    # Para adicionar estrelas de significância

# 1. Preparação dos Dados (Seu data.frame 'data' já está carregado)
# Verifique se as colunas existem e se são numéricas.  Adapte os nomes se necessário.

# Imprima os nomes das colunas do seu dataframe para verificar os nomes corretos.
print(colnames(data))

# Certifique-se de que todas as colunas que você vai usar na correlação são numéricas.
# Se alguma coluna não for, converta-a usando as.numeric().
# Exemplo (substitua "NomeDaColuna" pelo nome real):
# data$NomeDaColuna <- as.numeric(data$NomeDaColuna)

data <- textshape::column_to_rownames(data, "barcode")

barcodes_low <- metadata$TCGA_barcode[metadata$steroid == "Steroid_Low"]
barcodes_low <- barcodes_low[is.na(barcodes_low)==F]

data_low <- data[rownames(data) %in% barcodes_low,]

barcodes_high <- metadata$TCGA_barcode[metadata$steroid == "Steroid_High"]
barcodes_high <- barcodes_high[is.na(barcodes_high) == F]

data_high <- data[rownames(data) %in% barcodes_high,]

# 2. Cálculo da Correlação e dos Valores-P

# Função para calcular correlações e valores-p
cor.test.p <- function(x, y) {
  result <- cor.test(x, y, method = "spearman") #Usando Spearman
  data.frame(
    correlation = result$estimate,
    p.value = result$p.value
  )
}

# Obtenha todas as combinações de nomes de colunas
cols <- colnames(data)
combinations <- expand.grid(Var1 = cols, Var2 = cols, stringsAsFactors = FALSE)

cols_low <- colnames(data_low)
combinations_low <- expand.grid(Var1 = cols_low, Var2 = cols_low, stringsAsFactors = FALSE)

cols_high <- colnames(data_high)
combinations_high <- expand.grid(Var1 = cols_high, Var2 = cols_high, stringsAsFactors = FALSE)

# Aplique o teste de correlação a cada combinação
results <- apply(combinations, 1, function(row) {
  var1 <- row[["Var1"]]
  var2 <- row[["Var2"]]
  cor.test.p(data[[var1]], data[[var2]])
})

results_low <- apply(combinations_low, 1, function(row) {
  var1 <- row[["Var1"]]
  var2 <- row[["Var2"]]
  cor.test.p(data_low[[var1]], data_low[[var2]])
})

results_high <- apply(combinations_high, 1, function(row) {
  var1 <- row[["Var1"]]
  var2 <- row[["Var2"]]
  cor.test.p(data_high[[var1]], data_high[[var2]])
})


# Combine os resultados em um data.frame
correlation_data <- data.frame(
  Var1 = combinations$Var1,
  Var2 = combinations$Var2,
  correlation = unlist(lapply(results, function(x) x$correlation)),
  p.value = unlist(lapply(results, function(x) x$p.value))
)

correlation_data_low <- data.frame(
  Var1 = combinations$Var1,
  Var2 = combinations$Var2,
  correlation = unlist(lapply(results_low, function(x) x$correlation)),
  p.value = unlist(lapply(results_low, function(x) x$p.value))
)

correlation_data_high <- data.frame(
  Var1 = combinations$Var1,
  Var2 = combinations$Var2,
  correlation = unlist(lapply(results_high, function(x) x$correlation)),
  p.value = unlist(lapply(results_high, function(x) x$p.value))
)

correlation_data_low$p.value <- p.adjust(as.vector(correlation_data_low$p.value), method = "BH")
correlation_data_high$p.value <- p.adjust(as.vector(correlation_data_high$p.value), method = "BH")


# Crie rótulos de significância
correlation_data$significance <- ifelse(correlation_data$p.value < 0.05, "*", "")
correlation_data_low$significance <- ifelse(correlation_data_low$p.value < 0.05, "*", "")
correlation_data_high$significance <- ifelse(correlation_data_high$p.value < 0.05, "*", "")

# 2.1 Remover a linha Var1 == Var2

correlation_data <- correlation_data[correlation_data$Var1 != correlation_data$Var2,]
correlation_data_low <- correlation_data_low[correlation_data_low$Var1 != correlation_data_low$Var2,]
correlation_data_high <- correlation_data_high[correlation_data_high$Var1 != correlation_data_high$Var2,]

# 2.2 Remover as duplicadas (Var1 e Var2 invertidas)
correlation_data$pair <- apply(correlation_data[c("Var1", "Var2")], 1, function(x) paste(sort(x), collapse = "_"))
correlation_data <- correlation_data[!duplicated(correlation_data$pair),]

correlation_data_low$pair <- apply(correlation_data_low[c("Var1", "Var2")], 1, function(x) paste(sort(x), collapse = "_"))
correlation_data_low <- correlation_data_low[!duplicated(correlation_data_low$pair),]

correlation_data_high$pair <- apply(correlation_data_high[c("Var1", "Var2")], 1, function(x) paste(sort(x), collapse = "_"))
correlation_data_high <- correlation_data_high[!duplicated(correlation_data_high$pair),]


# 2.3 Preparar os dados para o gráfico
plot_data <- correlation_data
plot_data_low <- correlation_data_low
plot_data_high <- correlation_data_high

# 3. Ordenação das Variáveis (ADAPTE ESTA PARTE PARA A SUA ORDEM DESEJADA)
# Defina a ordem das variáveis como na Figura B (adapte os nomes se necessário!)
variaveis_y <- c("TRD_Entropy","TRG_Entropy","TRB_Entropy","TRA_Entropy",
                 "TRD_Abundance","TRG_Abundance","TRB_Abundance","TRA_Abundance",
                 "IGL_Entropy","IGK_Entropy","IGH_Entropy","IGL_Abundance",
                 "IGK_Abundance","IGH_Abundance")

#variaveis_y <- c("IGH_Abundance", "IGK_Abundance", "IGL_Abundance",
#                   "IGH_Entropy", "IGK_Entropy", "IGL_Entropy",
#                   "TRA_Abundance", "TRB_Abundance", "TRA_Entropy", "TRB_Entropy")

variaveis_x <- c("Lymphocyte.Infiltration", "Macrophage.Regulation",
                 "TGF.beta.Response", "IFN.gamma.Response", "Th1.Cells", "Th17.Cells",
                 "Th2.Cells")

#variaveis_x <- c("Lymphocyte.Infiltration", "Macrophage.Regulation",
#                 "TGF.beta.Response", "IFN.gamma.Response", "Th1.Cells", "Th17.Cells",
#                 "Th2.Cells", "Wound.Healing",
#                 "Nonsilent.Mutation.Rate", "Silent.Mutation.Rate",
#                 "SNV.Neoantigens", "paper_purity")

# Converta para fatores e defina a ordem
plot_data$Var1 <- factor(plot_data$Var1, levels = variaveis_y)
plot_data$Var2 <- factor(plot_data$Var2, levels = variaveis_x, ordered = TRUE)

plot_data_low$Var1 <- factor(plot_data_low$Var1, levels = variaveis_y)
plot_data_low$Var2 <- factor(plot_data_low$Var2, levels = variaveis_x, ordered = TRUE)

plot_data_high$Var1 <- factor(plot_data_high$Var1, levels = variaveis_y)
plot_data_high$Var2 <- factor(plot_data_high$Var2, levels = variaveis_x, ordered = TRUE)


plot_data <- plot_data[!is.na(plot_data$Var1),]
plot_data <- plot_data[!is.na(plot_data$Var2),]

plot_data_low <- plot_data_low[!is.na(plot_data_low$Var1),]
plot_data_low <- plot_data_low[!is.na(plot_data_low$Var2),]

plot_data_high <- plot_data_high[!is.na(plot_data_high$Var1),]
plot_data_high <- plot_data_high[!is.na(plot_data_high$Var2),]

# 4. Geração do Heatmap com ggplot2
heatmap_plot  <- ggplot(plot_data, aes(x = Var2, y = Var1, fill = correlation)) +
  geom_tile(color = "white") +  # Adicione bordas brancas às peças
  scale_fill_gradient2(low = "steelblue", high = "firebrick", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlação\nde Spearman") +
  theme_minimal() +  # Use um tema minimalista para uma aparência mais limpa
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 10, hjust = 1),
        axis.title.x = element_blank(), # Remova o título do eixo x
        axis.title.y = element_blank(),  # Remova o título do eixo y
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + # Remova características do painel
  geom_text(aes(label = significance), color = "black", size = 6) +  # Adicione as estrelas de significância
  coord_fixed()  # Garanta que as peças sejam quadradas

# 5. Salvar o gráfico em um arquivo
ggsave("heatmap_spearman_geral_20250715.pdf", plot = heatmap_plot, width = 6, height = 4, units = "in", dpi = 300)


heatmap_plot_low  <- ggplot(plot_data_low, aes(x = Var2, y = Var1, fill = correlation)) +
  geom_tile(color = "white") +  # Adicione bordas brancas às peças
  scale_fill_gradient2(low = "steelblue", high = "firebrick", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlação\nde Spearman") +
  theme_minimal() +  # Use um tema minimalista para uma aparência mais limpa
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 10, hjust = 1),
        axis.title.x = element_blank(), # Remova o título do eixo x
        axis.title.y = element_blank(),  # Remova o título do eixo y
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + # Remova características do painel
  geom_text(aes(label = significance), color = "black", size = 6) +  # Adicione as estrelas de significância
  coord_fixed()  # Garanta que as peças sejam quadradas

# 5. Salvar o gráfico em um arquivo
ggsave("heatmap_spearman_low_20250715.pdf", plot = heatmap_plot_low, width = 6, height = 4, units = "in", dpi = 300)


heatmap_plot_high  <- ggplot(plot_data_high, aes(x = Var2, y = Var1, fill = correlation)) +
  geom_tile(color = "white") +  # Adicione bordas brancas às peças
  scale_fill_gradient2(low = "steelblue", high = "firebrick", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlação\nde Spearman") +
  theme_minimal() +  # Use um tema minimalista para uma aparência mais limpa
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 10, hjust = 1),
        axis.title.x = element_blank(), # Remova o título do eixo x
        axis.title.y = element_blank(),  # Remova o título do eixo y
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + # Remova características do painel
  geom_text(aes(label = significance), color = "black", size = 6) +  # Adicione as estrelas de significância
  coord_fixed()  # Garanta que as peças sejam quadradas

# 5. Salvar o gráfico em um arquivo
ggsave("heatmap_spearman_high_20250715.pdf", plot = heatmap_plot_high, width = 6, height = 4, units = "in", dpi = 300)

