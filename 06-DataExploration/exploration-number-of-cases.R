df_data_clones <- read.csv("df_data_clones.csv")
df_data_reads <- read.csv("df_data_reads.csv")

head(df_data_reads)

sum(df_data_reads$bcr_counts)
sum(df_data_reads$tcr_counts)

table(df_data_clones$tumor_stage[df_data_clones$steroid=="Steroid_High"])
table(df_data_clones$tumor_stage[df_data_clones$steroid=="Steroid_Low"])

length(df_data_reads$TCGA_barcode[df_data_reads$bcr_counts == 0
                                  & df_data_reads$tcr_counts == 0])

bcr_tcr_zero <- df_data_reads[df_data_reads$bcr_counts == 0
                                      & df_data_reads$tcr_counts == 0,]

bcr_tcr <- df_data_reads[df_data_reads$TCGA_barcode %in% 
                           bcr_tcr_zero$TCGA_barcode == F,]

# quantos nao apresentaram sequencias de BCR e TCR?
length(bcr_tcr_zero$TCGA_barcode) # 10

# quantos casos apresentaram BCR e ou TCR?
length(bcr_tcr$TCGA_barcode) # 68

# quantos casos apresentaram apenas sequencias de BCR?
length(bcr_tcr$TCGA_barcode[bcr_tcr$bcr_counts != 0 & 
                              bcr_tcr$tcr_counts == 0]) # 12

# quantos casos apresentaram apenas sequencias de TCR?
length(bcr_tcr$TCGA_barcode[bcr_tcr$tcr_counts != 0 & 
                              bcr_tcr$bcr_counts == 0]) # 8

# quantos casos apresentaram sequencias de TCR e BCR?
length(bcr_tcr$TCGA_barcode[bcr_tcr$tcr_counts != 0 & 
                              bcr_tcr$bcr_counts != 0]) # 8


