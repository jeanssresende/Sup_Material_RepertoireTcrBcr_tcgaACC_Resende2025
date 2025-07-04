
library(ComplexHeatmap)

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")
dir_metrics <- "metricsTrust4/results_pipeline_report/"
load("../6-dataExploration/tcgaACC_pre_processed.RData")

# coldata
data <- read.csv("../6-dataExploration/df_data_clones.csv")
coldata <- data[, c("sample_id","TCGA_barcode","steroid","bcr_counts","tcr_counts",
                    "sum_bcr_tcr","reads")]
coldata <- coldata %>%
  arrange(desc(steroid), desc(sum_bcr_tcr))

idx <- match(coldata$TCGA_barcode, tcgaACC$barcode)

coldata$paper_C1A.C1B <- tcgaACC$paper_C1A.C1B[idx]
coldata$Cortisol <- tcgaACC$Cortisol[idx]
coldata$Immune.Subtype <- tcgaACC$Immune.Subtype[idx]
coldata$vital_status <- tcgaACC$vital_status[idx]
coldata$Lymphocyte.Inf.Score <- tcgaACC$Lymphocyte.Infiltration.Signature.Score[idx]
coldata$Leukocyte.Fraction <- tcgaACC$Leukocyte.Fraction[idx]
coldata$tumor_stage <- tcgaACC$tumor_stage[idx]

# -- cores dos metadados
col.bin.steroid <- c("Steroid_High"="#e41a1c", "Steroid_Low"="#377eb8")
col.imm_sub <- c("C1"="#4daf4a", "C2"="#984ea3", "C3"="#377eb8",
                 "C4"="#e41a1c", "C5"="#ff7f00", "C6"="#ffff33")
col.cAB <- c("C1A"="#e41a1c", "C1B"="#377eb8")
col.bin.vital_status <- c("Alive"="#525252", "Dead"="black")
col.stage <- c("stage i"= "#bdbdbd", "stage ii"="#969696",
               "stage iii"="#353535","stage iv"="#000000",
               "not reported"="#ffffff")
col.bin.cortisol <- c("Cortisol"="#e41a1c","No"="#377eb8")

col.seq.infilt <- colorRamp2(breaks = seq(0, 0.4, length.out=9),
                             colors = (brewer.pal(9,"Reds")))

col.seq.scores <- colorRamp2(breaks = seq(-4, 4, length.out=11),
                             colors = rev(brewer.pal(11,"RdGy")))

col.ha <- HeatmapAnnotation(steroid = coldata$steroid, #zheng
                            Cortisol = coldata$Cortisol, #zheng
                            C1A.C1B = coldata$paper_C1A.C1B, #zheng
                            vital_status = coldata$vital_status, #zheng
                            tumor_stage = coldata$tumor_stage, #zheng
                            Immune.Subtype = coldata$Immune.Subtype, #thorsson
                            Leukocyte.Fraction = coldata$Leukocyte.Fraction, #thorsson
                            Lymphocyte.Inf.Score = 
                              coldata$Lymphocyte.Inf.Score, #thorsson
                            col = list(steroid=col.bin.steroid,
                                       Cortisol = col.bin.cortisol,
                                       C1A.C1B = col.cAB,
                                       vital_status = col.bin.vital_status,
                                       tumor_stage = col.stage,
                                       Immune.Subtype = col.imm_sub,
                                       Leukocyte.Fraction = col.seq.infilt,
                                       Lymphocyte.Inf.Score = col.seq.scores))



arquivos <- list.files(dir_metrics, pattern = "\\.tsv$", full.names = TRUE)

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_abundance <- data.frame(t(dados$Abundance))
  rownames(dados_abundance) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_abundance)
}

colnames(data_metrics) <- dados$chain

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
data_metrics <- as.matrix(data_metrics)

# aplicando log2 e transformacao
data_metrics <- log2(data_metrics + 1)
data_metrics <- scale(data_metrics)

# estrutura
data_metrics <- t(data_metrics)

df <- as.data.frame(data_metrics)
df$cell_type <- c(rep("B",3), rep("T",4))

# setando escala
data_metrics[data_metrics > 2] <- 2
data_metrics[data_metrics < (-2)] <- (-2)


idx <- match(coldata$sample_id, colnames(data_metrics))

data_metrics <- data_metrics[,idx]

g_abundance <- Heatmap(data_metrics,
                      split = df$cell_type,
                      column_split = coldata$steroid,
                      top_annotation = col.ha,
                      show_column_names = FALSE,
                      cluster_columns = FALSE,
                      cluster_rows = F,
                      col=colorRamp2(breaks = seq(-2,2, length.out=9),
                                     colors = rev(brewer.pal(9,"RdYlBu"))),
                      row_title_gp = gpar(fontsize=10),
                      row_title_side = "left",
                      row_names_side = "right",
                      row_names_gp = gpar(fontsize=10))


# cpk

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_cpk <- data.frame(t(dados$CPK))
  rownames(dados_cpk) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_cpk)
}

colnames(data_metrics) <- dados$chain

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
data_metrics <- as.matrix(data_metrics)

data_metrics[is.na(data_metrics)] <- 0

# aplicando log2 e transformacao
data_metrics <- log2(data_metrics + 1)
data_metrics <- scale(data_metrics)

# estrutura
data_metrics <- t(data_metrics)

# setando escala
data_metrics[data_metrics > 2] <- 2
data_metrics[data_metrics < (-2)] <- (-2)


idx <- match(coldata$sample_id, colnames(data_metrics))

data_metrics <- data_metrics[,idx]

g_cpk <- Heatmap(data_metrics,
                       split = df$cell_type,
                       column_split = coldata$steroid,
                       #top_annotation = col.ha,
                       show_column_names = FALSE,
                       cluster_columns = FALSE,
                       cluster_rows = F,
                       col=colorRamp2(breaks = seq(-2,2, length.out=9),
                                      colors = rev(brewer.pal(9,"RdYlBu"))),
                       row_title_gp = gpar(fontsize=10),
                       row_title_side = "left",
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize=10))

# Entropy

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_entropy <- data.frame(t(dados$Entropy))
  rownames(dados_entropy) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_entropy)
}

colnames(data_metrics) <- dados$chain

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
data_metrics <- as.matrix(data_metrics)

data_metrics[is.na(data_metrics)] <- 0

# aplicando log2 e transformacao
data_metrics <- log2(data_metrics + 1)
data_metrics <- scale(data_metrics)

# estrutura
data_metrics <- t(data_metrics)

# setando escala
data_metrics[data_metrics > 2] <- 2
data_metrics[data_metrics < (-2)] <- (-2)


idx <- match(coldata$sample_id, colnames(data_metrics))

data_metrics <- data_metrics[,idx]

g_entropy <- Heatmap(data_metrics,
                 split = df$cell_type,
                 column_split = coldata$steroid,
                 #top_annotation = col.ha,
                 show_column_names = FALSE,
                 cluster_columns = FALSE,
                 cluster_rows = F,
                 col=colorRamp2(breaks = seq(-2,2, length.out=9),
                                colors = rev(brewer.pal(9,"RdYlBu"))),
                 row_title_gp = gpar(fontsize=10),
                 row_title_side = "left",
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize=10))


# Clonality

# montagem do objeto
data_metrics <- data.frame()

for (arquivo in arquivos) {
  dados <- read.table(arquivo, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  dados_clonality <- data.frame(t(dados$Clonality))
  rownames(dados_clonality) <- gsub("_clones.tsv", "", basename(arquivo))
  
  data_metrics <- rbind(data_metrics, dados_clonality)
}

colnames(data_metrics) <- dados$chain

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
data_metrics <- as.matrix(data_metrics)

data_metrics[is.na(data_metrics)] <- 0

# aplicando log2 e transformacao
data_metrics <- log2(data_metrics + 1)
data_metrics <- scale(data_metrics)

# estrutura
data_metrics <- t(data_metrics)

# setando escala
data_metrics[data_metrics > 2] <- 2
data_metrics[data_metrics < (-2)] <- (-2)


idx <- match(coldata$sample_id, colnames(data_metrics))

data_metrics <- data_metrics[,idx]

g_clonality <- Heatmap(data_metrics,
                     split = df$cell_type,
                     column_split = coldata$steroid,
                     #top_annotation = col.ha,
                     show_column_names = FALSE,
                     cluster_columns = FALSE,
                     cluster_rows = F,
                     col=colorRamp2(breaks = seq(-2,2, length.out=9),
                                    colors = rev(brewer.pal(9,"RdYlBu"))),
                     row_title_gp = gpar(fontsize=10),
                     row_title_side = "left",
                     row_names_side = "right",
                     row_names_gp = gpar(fontsize=10))

ht_list <- g_abundance %v% g_cpk %v% g_entropy %v% g_clonality


