library(dplyr)
library(textshape)
library(IOBR)
library(readr)    
library(tibble)   
library(corrplot) 
library(ggplot2)  
library(tidyr)  

markers <- readxl::read_xlsx("markers_TB.xlsx")
head(markers)

length(table(markers$marker))

data_clones <- read.csv("../../6-dataExploration/df_data_clones.csv")

eset_acc <- read.delim("../lymphocyteFraction/TCGA-ACC.star_tpm.tsv")

## -- anotacao e normalizacao em TPM

# Remover a versão do Ensembl ID
eset_acc$Ensembl_ID <- sub("\\.[0-9]+$", "", eset_acc$Ensembl_ID)

# Agrupar pelo Ensembl_ID e calcular a média
eset_acc_agg <- eset_acc %>%
  group_by(Ensembl_ID) %>%
  summarise(across(everything(), median))

# Transformar Ensembl_ID em nomes de linha
eset_acc_final <- column_to_rownames(eset_acc_agg, "Ensembl_ID")

mat_tcgaACC <- eset_acc_final

marker_TB <- unique(markers$EnsemblGene)

data <- mat_tcgaACC[rownames(mat_tcgaACC) %in% marker_TB,]
data <- t(data)

data_clones_TB <- data_clones[,c(4:10,14,15)]
rownames(data_clones_TB) <- data_clones_TB$TCGA_barcode
data_clones_TB <- data_clones_TB[,-8]

rownames(data) <- gsub("\\.", "-", rownames(data))
rownames(data_clones_TB) <- substring(rownames(data_clones_TB),1,16)
data <- data[-64,]

idx <- match(rownames(data), rownames(data_clones_TB))
data_clones_TB <- data_clones_TB[idx,]

head(rownames(data_clones_TB))
head(rownames(data))

data_TB <- cbind(data, data_clones_TB)

head(data_TB)

# --- 1. Preparar o 'data_TB' ---
# Renomear a coluna de fenótipo para algo mais claro
data_TB <- data_TB %>%
  rename(Steroid_Phenotype = steroid)

# Identificar as colunas para marcadores e cadeias de TCR/BCR
marker_cols <- colnames(data_TB)[1:58] # Colunas 2 a 38 (ajuste se a 1a coluna for ID/Barcode)
repertoire_cols <- colnames(data_TB)[59:65] # Colunas 39 a 45

# --- 2. Separar os dados por Fenótipo Esteroidal ---
LSP_data <- data_TB %>% filter(Steroid_Phenotype == "Steroid_Low")
HSP_data <- data_TB %>% filter(Steroid_Phenotype == "Steroid_High")

# --- 3. Função Auxiliar para Calcular Correlação e P-valores (com correção) ---
# Esta função encapsula o cálculo de correlação e p-valores ajustados
# para um par de data frames (X = marcadores, Y = repertório)
calculate_correlations <- function(df, marker_vars, repertoire_vars) {
  # Seleciona apenas as colunas relevantes
  df_markers <- df %>% select(all_of(marker_vars))
  df_repertoire <- df %>% select(all_of(repertoire_vars))
  
  # Aplica transformação log2(x+1) para todas as colunas numéricas
  # (É importante fazer isso APENAS se os dados originais não estiverem transformados e tiverem zeros)
  # Verifique a natureza dos seus dados. Se forem já scores de deconvolução, talvez não precise.
  # Para expressão gênica crua e contagens de repertório, é geralmente recomendado.
  df_markers_transformed <- df_markers #%>% mutate(across(everything(), ~log2(.x + 1)))
  df_repertoire_transformed <- df_repertoire %>% mutate(across(everything(), ~log2(.x + 1)))
  
  # Calcula a matriz de correlação de Spearman
  cor_matrix <- cor(df_markers_transformed, df_repertoire_transformed, method = "spearman", use = "pairwise.complete.obs")
  
  # Calcula a matriz de p-valores
  p_matrix <- matrix(NA, nrow(cor_matrix), ncol(cor_matrix))
  rownames(p_matrix) <- rownames(cor_matrix)
  colnames(p_matrix) <- colnames(cor_matrix)
  
  for (i in 1:nrow(cor_matrix)) {
    for (j in 1:ncol(cor_matrix)) {
      test_result <- cor.test(df_markers_transformed[[rownames(cor_matrix)[i]]],
                              df_repertoire_transformed[[colnames(cor_matrix)[j]]],
                              method = "spearman", exact = FALSE) # exact=FALSE para amostras maiores
      p_matrix[i, j] <- test_result$p.value
    }
  }
  
  # Ajusta os p-valores para múltiplas comparações (Benjamini-Hochberg FDR)
  adj_p_matrix <- p.adjust(p_matrix, method = "BH")
  dim(adj_p_matrix) <- dim(p_matrix) # Reshape de volta para matriz
  rownames(adj_p_matrix) <- rownames(p_matrix)
  colnames(adj_p_matrix) <- colnames(p_matrix)
  
  return(list(cor = cor_matrix, p = p_matrix, p_adj = adj_p_matrix))
}


# --- 4. Calcular Correlações para LSP e HSP ---

cat("\nCalculando correlações para o grupo LSP...\n")
LSP_cor_results <- calculate_correlations(LSP_data, marker_cols, repertoire_cols)

cat("\nCalculando correlações para o grupo HSP...\n")
HSP_cor_results <- calculate_correlations(HSP_data, marker_cols, repertoire_cols)

# --- 5. Visualizar os Mapas de Calor ---

# Esquema de cores para o heatmap
my_colors <- colorRampPalette(c("blue", "white", "red"))(200)

# Mapa de Calor para LSP
cat("\nGerando heatmap para o grupo LSP...\n")
corrplot(LSP_cor_results$cor,
         p.mat = LSP_cor_results$p_adj,
         sig.level = 0.05,
         insig = "label_sig", 
         method = "color",
         type = "full",   
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.8,    
         cl.cex = 0.8,    
         col = my_colors,
         title = "Correlações (Spearman) em ACC-LSP: Marcadores vs Repertório",
         mar=c(0,0,3,0) 
)

# Mapa de Calor para HSP
cat("\nGerando heatmap para o grupo HSP...\n")
corrplot(HSP_cor_results$cor,
         p.mat = HSP_cor_results$p_adj,
         sig.level = 0.05,
         insig = "label_sig", 
         method = "color",
         type = "full",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.8,
         cl.cex = 0.8,
         col = my_colors,
         title = "Correlações (Spearman) em ACC-HSP: Marcadores vs Repertório",
         mar=c(0,0,3,0) 
)