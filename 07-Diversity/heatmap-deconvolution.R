library(dplyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)


load("../6-dataExploration/tcgaACC_pre_processed.RData")
data <- read.csv("../6-dataExploration/df_data_clones.csv")

# dados para o heatmap
mat_data <- data[,c("IGH","IGK","IGL","TRA","TRB","TRD","TRG")]
rownames(mat_data) <- data$TCGA_barcode

# coldata
coldata <- data[, c("TCGA_barcode","steroid","bcr_counts","tcr_counts",
                    "sum_bcr_tcr","reads")]

idx <- match(coldata$TCGA_barcode, tcgaACC$barcode)

coldata$Lymphocyte.Infiltration.Signature.Score <- 
  tcgaACC$Lymphocyte.Infiltration.Signature.Score[idx]

coldata$Leukocyte.Fraction <- tcgaACC$Leukocyte.Fraction[idx]


met_deconv_tcrbcr <- read.csv("lymphocyteFraction/metricas_deconv_tcrbcr.csv")
met_deconv_tcrbcr$X <- gsub("\\.", "-", met_deconv_tcrbcr$X)

idx <- match(substring(coldata$TCGA_barcode,1,16), met_deconv_tcrbcr$X)  

met_deconv_tcrbcr <- met_deconv_tcrbcr[idx,]
rownames(met_deconv_tcrbcr) <- met_deconv_tcrbcr$X
met_deconv_tcrbcr <- met_deconv_tcrbcr[,-1]

mat_data <- met_deconv_tcrbcr[,4:19]

# ordem
coldata <- coldata %>%
  arrange(desc(steroid), desc(sum_bcr_tcr))

idx <- match(substring(coldata$TCGA_barcode,1,16), rownames(mat_data))

mat_data <- mat_data[idx,]

data.scale <- scale(mat_data )

#colMeans(data.log2.scale)
#rowMeans(data.log2.scale, na.rm = T)

#sd(mat_data[1,])

data <- (data.scale)

#data <- as.matrix(data)

data[data > 2] <- 2
data[data < (-2)] <- (-2)

rownames(coldata) <- substring(coldata$TCGA_barcode,1,16)


# metadados
coldata$bcr_counts_log10 <- log10(coldata$bcr_counts+1)
coldata$tcr_counts_log10 <- log10(coldata$tcr_counts+1)
coldata$sum_bcr_tcr_log10 <- log10(coldata$sum_bcr_tcr+1)
coldata$reads_log10 <- log10(coldata$reads+1)



# -- cores dos metadados
col.bin.steroid <- c("Steroid_High"="#e41a1c", "Steroid_Low"="#377eb8")

col.div.scores <- colorRamp2(breaks = seq(-4, 2, length.out=11),
                             colors = rev(brewer.pal(11,"RdGy")))

col.seq.leukocyte <- colorRamp2(breaks = seq(0, 0.4, length.out=9),
                                colors = (brewer.pal(9,"Reds")))

col.seq_greys.read <- colorRamp2(breaks = seq(21000000, 83000000, length.out=9),
                                 colors = brewer.pal(9,"Greys"))

col.seq_greys <- colorRamp2(breaks = seq(0, 5, length.out=9),
                            colors = brewer.pal(9,"Greys"))

col.ha <- HeatmapAnnotation(steroid = coldata$steroid,
                            Leukocyte.Fraction = coldata$Leukocyte.Fraction,
                            Lymphocyte.Inf.Score = 
                              coldata$Lymphocyte.Infiltration.Signature.Score,
                            count_TCR_log10 = coldata$tcr_counts_log10,
                            count_BCR_log10 = coldata$bcr_counts_log10,
                            count_reads = coldata$reads,
                            col = list(steroid = col.bin.steroid,
                                       Leukocyte.Fraction = col.seq.leukocyte,
                                       Lymphocyte.Inf.Score = col.div.scores,
                                       count_TCR_log10 = col.seq_greys,
                                       count_BCR_log10 = col.seq_greys,
                                       count_reads = col.seq_greys.read))


# objeto contendo a informacao de TCR e BCR
#row.ha <- rowAnnotation(type=c(rep("B",3), rep("T",4)))
#df <- as.data.frame(data)
#df$cell_type <- c(rep("B",3), rep("T",4))                            

#data <- as.matrix(data)
data <- t(data)

head(rownames(coldata)) 
head(colnames(data))

g_xcell <- Heatmap(data,
                       #split = df$cell_type,
                       column_split = coldata$steroid,
                       #top_annotation = col.ha,
                       show_column_names = FALSE,
                       cluster_columns = FALSE,
                       cluster_rows = F,
                       col=colorRamp2(breaks = seq(-2,2, length.out=11),
                                      colors = rev(brewer.pal(11,"RdYlBu"))),
                       row_title_gp = gpar(fontsize=10),
                       row_title_side = "left",
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize=10))                          

# -- mcpcounter
mat_data <- met_deconv_tcrbcr[,20:23]

idx <- match(substring(coldata$TCGA_barcode,1,16), rownames(mat_data))

mat_data <- mat_data[idx,]

data.scale <- scale(mat_data )

data <- (data.scale)

#data <- as.matrix(data)

data[data > 2] <- 2
data[data < (-2)] <- (-2)

data <- t(data)

g_mcpcounter <- Heatmap(data,
                   #split = df$cell_type,
                   column_split = coldata$steroid,
                   top_annotation = col.ha,
                   show_column_names = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   col=colorRamp2(breaks = seq(-2,2, length.out=11),
                                  colors = rev(brewer.pal(11,"RdYlBu"))),
                   row_title_gp = gpar(fontsize=10),
                   row_title_side = "left",
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize=10))

# -- timer
mat_data <- met_deconv_tcrbcr[,24:26]

idx <- match(substring(coldata$TCGA_barcode,1,16), rownames(mat_data))

mat_data <- mat_data[idx,]

data.scale <- scale(mat_data )

data <- (data.scale)

#data <- as.matrix(data)

data[data > 2] <- 2
data[data < (-2)] <- (-2)

data <- t(data)

g_timer <- Heatmap(data,
                        #split = df$cell_type,
                        column_split = coldata$steroid,
                        #top_annotation = col.ha,
                        show_column_names = FALSE,
                        cluster_columns = FALSE,
                        cluster_rows = F,
                        col=colorRamp2(breaks = seq(-2,2, length.out=11),
                                       colors = rev(brewer.pal(11,"RdYlBu"))),
                        row_title_gp = gpar(fontsize=10),
                        row_title_side = "left",
                        row_names_side = "right",
                        row_names_gp = gpar(fontsize=10))

# -- epic
mat_data <- met_deconv_tcrbcr[,1:3]

idx <- match(substring(coldata$TCGA_barcode,1,16), rownames(mat_data))

mat_data <- mat_data[idx,]

data.scale <- scale(mat_data )

data <- (data.scale)

#data <- as.matrix(data)

data[data > 2] <- 2
data[data < (-2)] <- (-2)

data <- t(data)

head(rownames(coldata)) 
head(colnames(data))

g_epic <- Heatmap(data,
                   #split = df$cell_type,
                   column_split = coldata$steroid,
                   #top_annotation = col.ha,
                   show_column_names = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = F,
                   col=colorRamp2(breaks = seq(-2,2, length.out=11),
                                  colors = rev(brewer.pal(11,"RdYlBu"))),
                   row_title_gp = gpar(fontsize=10),
                   row_title_side = "left",
                   row_names_side = "right",
                   row_names_gp = gpar(fontsize=10))

# -- cibersort
mat_data <- met_deconv_tcrbcr[,27:35]

idx <- match(substring(coldata$TCGA_barcode,1,16), rownames(mat_data))

mat_data <- mat_data[idx,]

data.scale <- scale(mat_data )

data <- (data.scale)

#data <- as.matrix(data)

data[data > 2] <- 2
data[data < (-2)] <- (-2)

data <- t(data)

head(rownames(coldata)) 
head(colnames(data))

g_cibersort <- Heatmap(data,
                  #split = df$cell_type,
                  column_split = coldata$steroid,
                  #top_annotation = col.ha,
                  show_column_names = FALSE,
                  cluster_columns = FALSE,
                  cluster_rows = F,
                  col=colorRamp2(breaks = seq(-2,2, length.out=11),
                                 colors = rev(brewer.pal(11,"RdYlBu"))),
                  row_title_gp = gpar(fontsize=10),
                  row_title_side = "left",
                  row_names_side = "right",
                  row_names_gp = gpar(fontsize=10))

# -- quantiseq
mat_data <- met_deconv_tcrbcr[,36:39]

idx <- match(substring(coldata$TCGA_barcode,1,16), rownames(mat_data))

mat_data <- mat_data[idx,]

data.scale <- scale(mat_data )

data <- (data.scale)

#data <- as.matrix(data)

data[data > 2] <- 2
data[data < (-2)] <- (-2)

data <- t(data)

head(rownames(coldata)) 
head(colnames(data))

g_quantiseq <- Heatmap(data,
                       #split = df$cell_type,
                       column_split = coldata$steroid,
                       #top_annotation = col.ha,
                       show_column_names = FALSE,
                       cluster_columns = FALSE,
                       cluster_rows = F,
                       col=colorRamp2(breaks = seq(-2,2, length.out=11),
                                      colors = rev(brewer.pal(11,"RdYlBu"))),
                       row_title_gp = gpar(fontsize=10),
                       row_title_side = "left",
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize=10))

panel_1 <- g_mcpcounter %v% g_timer %v% g_xcell %v% g_epic %v% g_quantiseq %v% g_cibersort
















# -- cores dos metadados
col.bin.steroid <- c("Steroid_High"="#e41a1c", "Steroid_Low"="#377eb8")
col.bin.cortisol <- c("Cortisol"="#e41a1c","No"="#377eb8")
col.imm_sub <- c("C1"="#4daf4a", "C2"="#984ea3", "C3"="#377eb8",
                 "C4"="#e41a1c", "C5"="#ff7f00", "C6"="#ffff33")
col.cAB <- c("C1A"="#e41a1c", "C1B"="#377eb8")
col.bin.vital_status <- c("Alive"="#525252", "Dead"="black")

#col.stage <- c("stage i"= "#bdbdbd", "stage ii"="#969696",
#               "stage iii"="#353535","stage iv"="#000000",
#               "not reported"="#ffffff")
#col.oth_horm <- c("Mineralcorticoids"="#7fc97f", "Sexual"="#beaed4", "No"="#525252")

col.seq.leukocyte <- colorRamp2(breaks = seq(0, 0.4, length.out=9),
                                colors = (brewer.pal(9,"YlOrRd")))

col.div.scores <- colorRamp2(breaks = seq(-4, 4, length.out=11),
                             colors = rev(brewer.pal(11,"RdYlBu")))

col.seq_greys.read <- colorRamp2(breaks = seq(21000000, 83000000, length.out=9),
                                 colors = brewer.pal(9,"Greys"))

col.seq_greys <- colorRamp2(breaks = seq(0, 5, length.out=9),
                            colors = brewer.pal(9,"Greys"))

# objeto contendo a informacao de TCR e BCR
row.ha <- rowAnnotation(type=c(rep("B",3), rep("T",4)))
df <- as.data.frame(data)
df$cell_type <- c(rep("B",3), rep("T",4))






# coldata
coldata <- data[, c("TCGA_barcode", "steroid","cortisol.excess","immune.subtype",
                    "other.hormones","gender","vital_status","tumor_stage",
                    "bcr_counts","tcr_counts","sum_bcr_tcr","reads")]
