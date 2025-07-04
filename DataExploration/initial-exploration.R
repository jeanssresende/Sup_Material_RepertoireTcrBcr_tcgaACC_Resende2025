library(immunarch)
library(dplyr)
library(purrr)
library(tidyr)

# - carregando os dados
immdata_reads <- repLoad("../3-extractionCDR3/cdr3/")
immdata_clones <- repLoad("../5-clonesIdentification/data/")

#sum(map_dbl(immdata_reads$data, ~sum(.x$Clones, na.rm = TRUE)))
#sum(map_dbl(immdata_clones$data, ~sum(.x$Clones, na.rm = TRUE)))

unique(gsub("_report","",names(immdata_reads$data)))
unique(gsub("_clones","",names(immdata_clones$data)))

immdata_reads$data$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`
immdata_clones$data$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_clones`



# Unindo todos os tibble de reads em um único dataframe
immdata_reads <- map(immdata_reads$data, as_tibble)
df_reads <- bind_rows(immdata_reads, .id = "sample_id")

# Unindo todos os tibble de clones em um único dataframe
immdata_clones <- map(immdata_clones$data, as_tibble)
df_clones <- bind_rows(immdata_clones, .id = "sample_id")

#unique(df_reads$sample_id)
#unique(df_clones$sample_id)

# filtrando duplicatas indesejadas
temp <- df_reads %>%
  group_by_all() %>%
  filter(n() > 1)

rm(temp)

# salvando as tabelas
#write.csv(df_reads, file = "df_reads.csv")
#write.csv(df_clones, file = "df_clones.csv")

# rmetadados
load("../4-IDs-metadata/metadata_2025.RData")

metadata$Sample <- paste(metadata$file_id,"_report", sep = "")

amostras_zeradas <- metadata$Sample[metadata$Sample %in% df_reads$sample_id == FALSE] 

df_amostras_zeradas <- data.frame(
  sample_id = amostras_zeradas,
  Clones = rep(0,length(amostras_zeradas)),
  Proportion = rep(0,length(amostras_zeradas)),
  CDR3.nt = "",
  CDR3.aa = "",
  V.name = "",
  D.name = "",
  J.name = "",
  V.end = "",
  D.start = "",
  D.end = "",
  J.start = "",
  VJ.ins = "",
  VD.ins = "",
  DJ.ins = "")

df_reads_all <- rbind(df_reads, df_amostras_zeradas)
df_clones_all <- rbind(df_clones, df_amostras_zeradas)

#df_reads_all <- df_reads_all[!is.na(df_reads_all$V.name),]
#df_clones_all <- df_clones_all[!is.na(df_clones_all$V.name),]

#unique(df_reads_all$sample_id) #79
#unique(df_clones_all$sample_id) #79

idx <- match(df_reads_all$sample_id, metadata$Sample)
df_reads_metadados <- cbind(df_reads_all, metadata[idx,])
 
idx <- match(gsub("_clones","",df_clones_all$sample_id),
             gsub("_report","",metadata$Sample))
df_clones_metadados <- cbind(df_clones_all, metadata[idx,])

# Criando a contagem baseada na soma da coluna `clones`
#df_counts_clones <- df_clones_metadados %>%
#  mutate(chain_type = case_when(
#    startsWith(V.name, "TRA") ~ "TRA",
#    startsWith(V.name, "TRB") ~ "TRB",
#    startsWith(V.name, "TRD") ~ "TRD",
#    startsWith(V.name, "TRG") ~ "TRG",
#    startsWith(V.name, "IGH") ~ "IGH",
#    startsWith(V.name, "IGK") ~ "IGK",
#    startsWith(V.name, "IGL") ~ "IGL"
#  )) %>%
#  group_by(sample_id, chain_type) %>%
#  summarise(total_clones = sum(Clones, na.rm = TRUE), .groups = "drop")  # Somando clones

# Criando a contagem baseada na soma da coluna `clones`, considerando J.name quando V.name estiver NA
df_counts_clones <- df_clones_metadados %>%
  mutate(chain_type = case_when(
    is.na(V.name) & startsWith(J.name, "TRA") ~ "TRA",
    is.na(V.name) & startsWith(J.name, "TRB") ~ "TRB",
    is.na(V.name) & startsWith(J.name, "TRD") ~ "TRD",
    is.na(V.name) & startsWith(J.name, "TRG") ~ "TRG",
    is.na(V.name) & startsWith(J.name, "IGH") ~ "IGH",
    is.na(V.name) & startsWith(J.name, "IGK") ~ "IGK",
    is.na(V.name) & startsWith(J.name, "IGL") ~ "IGL",
    startsWith(V.name, "TRA") ~ "TRA",
    startsWith(V.name, "TRB") ~ "TRB",
    startsWith(V.name, "TRD") ~ "TRD",
    startsWith(V.name, "TRG") ~ "TRG",
    startsWith(V.name, "IGH") ~ "IGH",
    startsWith(V.name, "IGK") ~ "IGK",
    startsWith(V.name, "IGL") ~ "IGL",
    TRUE ~ "Unclassified"  # Caso nenhuma regra se aplique
  )) %>%
  group_by(sample_id, chain_type) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE), .groups = "drop")

# Transformar para formato de tabela final
df_final_clones <- df_counts_clones %>%
  pivot_wider(names_from = chain_type, values_from = total_clones, values_fill = 0)


#df_counts_reads <- df_reads_metadados %>%
#  mutate(chain_type = case_when(
#    startsWith(V.name, "TRA") ~ "TRA",
#    startsWith(V.name, "TRB") ~ "TRB",
#    startsWith(V.name, "TRD") ~ "TRD",
#    startsWith(V.name, "TRG") ~ "TRG",
#    startsWith(V.name, "IGH") ~ "IGH",
#    startsWith(V.name, "IGK") ~ "IGK",
#    startsWith(V.name, "IGL") ~ "IGL"
#  )) %>%
#  group_by(sample_id, chain_type) %>%
#  summarise(total_clones = sum(Clones, na.rm = TRUE), .groups = "drop")  # Somando clones

# Criando a contagem baseada na soma da coluna `clones`, considerando J.name quando V.name estiver NA
df_counts_reads <- df_reads_metadados %>%
  mutate(chain_type = case_when(
    is.na(V.name) & startsWith(J.name, "TRA") ~ "TRA",
    is.na(V.name) & startsWith(J.name, "TRB") ~ "TRB",
    is.na(V.name) & startsWith(J.name, "TRD") ~ "TRD",
    is.na(V.name) & startsWith(J.name, "TRG") ~ "TRG",
    is.na(V.name) & startsWith(J.name, "IGH") ~ "IGH",
    is.na(V.name) & startsWith(J.name, "IGK") ~ "IGK",
    is.na(V.name) & startsWith(J.name, "IGL") ~ "IGL",
    startsWith(V.name, "TRA") ~ "TRA",
    startsWith(V.name, "TRB") ~ "TRB",
    startsWith(V.name, "TRD") ~ "TRD",
    startsWith(V.name, "TRG") ~ "TRG",
    startsWith(V.name, "IGH") ~ "IGH",
    startsWith(V.name, "IGK") ~ "IGK",
    startsWith(V.name, "IGL") ~ "IGL",
    TRUE ~ "Unclassified"  # Caso nenhuma regra se aplique
  )) %>%
  group_by(sample_id, chain_type) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE), .groups = "drop")

# Transformar para formato de tabela final
df_final_reads <- df_counts_reads %>%
  pivot_wider(names_from = chain_type, values_from = total_clones, values_fill = 0)

df_final_reads$sample_id <- gsub("_report","",df_final_reads$sample_id)
df_final_clones$sample_id <- gsub("_clones","",df_final_clones$sample_id)
df_final_clones$sample_id <- gsub("_report","",df_final_clones$sample_id)

idx <- match(df_final_reads$sample_id, gsub("_report","",metadata$Sample) )
df_data_reads <- cbind(df_final_reads, metadata[idx,])

idx <- match(df_final_clones$sample_id, gsub("_report","",metadata$Sample) )
df_data_clones <- cbind(df_final_clones, metadata[idx,])

df_data_reads <- df_data_reads[,c(1:4,7,5,8:22)]
df_data_clones <- df_data_clones[,c(1:4,7,5,8:22)]

# exploracao
# Filtrar apenas as amostras do grupo Steroid_High e que são TCR (V.name começa com "TR")
df_reads_metadados %>%
  filter(steroid == "Steroid_High", startsWith(V.name, "TR")) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE))

df_reads_metadados %>%
  filter(steroid == "Steroid_Low", startsWith(V.name, "TR")) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE))

df_reads_metadados %>%
  filter(steroid == "Steroid_High", startsWith(V.name, "IG")) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE))

df_reads_metadados %>%
  filter(steroid == "Steroid_Low", startsWith(V.name, "IG")) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE))

# Contar quantas sequências únicas de TCR existem no grupo Steroid_High
df_reads_metadados %>%
  filter(steroid == "Steroid_High", startsWith(V.name, "TR")) %>%
  summarise(unique_sequences = n_distinct(CDR3.nt))

df_reads_metadados %>%
  filter(steroid == "Steroid_Low", startsWith(V.name, "TR")) %>%
  summarise(unique_sequences = n_distinct(CDR3.nt))

df_reads_metadados %>%
  filter(steroid == "Steroid_High", startsWith(V.name, "IG")) %>%
  summarise(unique_sequences = n_distinct(CDR3.nt))

df_reads_metadados %>%
  filter(steroid == "Steroid_Low", startsWith(V.name, "IG")) %>%
  summarise(unique_sequences = n_distinct(CDR3.nt))




# Filtrar apenas as amostras do grupo Steroid_High e que são TCR (V.name começa com "TR")
df_clones_metadados %>%
  filter(steroid == "Steroid_High", startsWith(V.name, "TR")) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE))

df_clones_metadados %>%
  filter(steroid == "Steroid_Low", startsWith(V.name, "TR")) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE))

df_clones_metadados %>%
  filter(steroid == "Steroid_High", startsWith(V.name, "IG")) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE))

df_clones_metadados %>%
  filter(steroid == "Steroid_Low", startsWith(V.name, "IG")) %>%
  summarise(total_clones = sum(Clones, na.rm = TRUE))

# Contar quantas sequências únicas de TCR existem no grupo Steroid_High
df_clones_metadados %>%
  filter(steroid == "Steroid_High", startsWith(V.name, "TR")) %>%
  summarise(unique_sequences = n_distinct(CDR3.nt))

df_clones_metadados %>%
  filter(steroid == "Steroid_Low", startsWith(V.name, "TR")) %>%
  summarise(unique_sequences = n_distinct(CDR3.nt))

df_clones_metadados %>%
  filter(steroid == "Steroid_High", startsWith(V.name, "IG")) %>%
  summarise(unique_sequences = n_distinct(CDR3.nt))

df_clones_metadados %>%
  filter(steroid == "Steroid_Low", startsWith(V.name, "IG")) %>%
  summarise(unique_sequences = n_distinct(CDR3.nt))



# -- table 1
table(df_data_clones$steroid)
table(df_data_clones$vital_status == "Dead" & df_data_clones$steroid == "Steroid_High")
table(df_data_clones$vital_status == "Dead" & df_data_clones$steroid == "Steroid_Low")

table(df_data_clones$steroid == "Steroid_High" & df_data_clones$gender == "male")
table(df_data_clones$steroid == "Steroid_High" & df_data_clones$gender == "male" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$gender == "male")
table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$gender == "male" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_High" & df_data_clones$gender == "female")
table(df_data_clones$steroid == "Steroid_High" & df_data_clones$gender == "female" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$gender == "female")
table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$gender == "female" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_High" & df_data_clones$tumor_stage == "stage i")
table(df_data_clones$steroid == "Steroid_High" & df_data_clones$tumor_stage == "stage i" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$tumor_stage == "stage i")
table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$tumor_stage == "stage i" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_High" & df_data_clones$tumor_stage == "stage ii")
table(df_data_clones$steroid == "Steroid_High" & df_data_clones$tumor_stage == "stage ii" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$tumor_stage == "stage ii")
table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$tumor_stage == "stage ii" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_High" & df_data_clones$tumor_stage == "stage iii")
table(df_data_clones$steroid == "Steroid_High" & df_data_clones$tumor_stage == "stage iii" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$tumor_stage == "stage iii")
table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$tumor_stage == "stage iii" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_High" & df_data_clones$tumor_stage == "stage iv")
table(df_data_clones$steroid == "Steroid_High" & df_data_clones$tumor_stage == "stage iv" & df_data_clones$vital_status == "Dead")

table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$tumor_stage == "stage iv")
table(df_data_clones$steroid == "Steroid_Low" & df_data_clones$tumor_stage == "stage iv" & df_data_clones$vital_status == "Dead")

df_data_clones <- df_data_clones[!is.na(df_data_clones$steroid),]
df_data_reads <- df_data_reads[!is.na(df_data_reads$steroid),]

load("tcgaACC_pre_processed.RData")

idx <- match(df_data_clones$TCGA_barcode, tcgaACC$barcode)
df_data_clones$age_at_diagnosis <-  tcgaACC$age_at_diagnosis[idx]

idx <- match(df_data_reads$TCGA_barcode, tcgaACC$barcode)
df_data_reads$age_at_diagnosis <-  tcgaACC$age_at_diagnosis[idx]

round(mean(df_data_clones$age_at_diagnosis[df_data_clones$steroid == "Steroid_High" &
                                             df_data_clones$gender == "male"]),1)
round(sd(df_data_clones$age_at_diagnosis[df_data_clones$steroid == "Steroid_High" &
                                             df_data_clones$gender == "male"]),1)

round(mean(df_data_clones$age_at_diagnosis[df_data_clones$steroid == "Steroid_Low" &
                                             df_data_clones$gender == "male"]),1)
round(sd(df_data_clones$age_at_diagnosis[df_data_clones$steroid == "Steroid_Low" &
                                           df_data_clones$gender == "male"]),1)

round(mean(df_data_clones$age_at_diagnosis[df_data_clones$steroid == "Steroid_High" &
                                             df_data_clones$gender == "female"]),1)
round(sd(df_data_clones$age_at_diagnosis[df_data_clones$steroid == "Steroid_High" &
                                           df_data_clones$gender == "female"]),1)

round(mean(df_data_clones$age_at_diagnosis[df_data_clones$steroid == "Steroid_Low" &
                                             df_data_clones$gender == "female"]),1)
round(sd(df_data_clones$age_at_diagnosis[df_data_clones$steroid == "Steroid_Low" &
                                           df_data_clones$gender == "female"]),1)

# -- Table S1 - Sequencing Summary
sum(df_data_reads$IGH[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGH != 0])
median(df_data_reads$IGH[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGH != 0])
range(df_data_reads$IGH[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGH != 0])
length(df_data_reads$IGH[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGH != 0])

sum(df_data_reads$IGH[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGH != 0])
median(df_data_reads$IGH[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGH != 0])
range(df_data_reads$IGH[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGH != 0])
length(df_data_reads$IGH[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGH != 0])

sum(df_data_reads$IGK[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGK != 0])
median(df_data_reads$IGK[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGK != 0])
range(df_data_reads$IGK[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGK != 0])
length(df_data_reads$IGK[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGK != 0])

sum(df_data_reads$IGK[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGK != 0])
median(df_data_reads$IGK[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGK != 0])
range(df_data_reads$IGK[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGK != 0])
length(df_data_reads$IGK[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGK != 0])

sum(df_data_reads$IGL[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGL != 0])
median(df_data_reads$IGL[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGL != 0])
range(df_data_reads$IGL[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGL != 0])
length(df_data_reads$IGL[df_data_reads$steroid == "Steroid_High" & df_data_reads$IGL != 0])

sum(df_data_reads$IGL[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGL != 0])
median(df_data_reads$IGL[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGL != 0])
range(df_data_reads$IGL[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGL != 0])
length(df_data_reads$IGL[df_data_reads$steroid == "Steroid_Low" & df_data_reads$IGL != 0])

sum(df_data_reads$TRA[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRA != 0])
median(df_data_reads$TRA[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRA != 0])
range(df_data_reads$TRA[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRA != 0])
length(df_data_reads$TRA[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRA != 0])

sum(df_data_reads$TRA[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRA != 0])
median(df_data_reads$TRA[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRA != 0])
range(df_data_reads$TRA[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRA != 0])
length(df_data_reads$TRA[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRA != 0])

sum(df_data_reads$TRB[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRB != 0])
median(df_data_reads$TRB[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRB != 0])
range(df_data_reads$TRB[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRB != 0])
length(df_data_reads$TRB[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRB != 0])

sum(df_data_reads$TRB[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRB != 0])
median(df_data_reads$TRB[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRB != 0])
range(df_data_reads$TRB[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRB != 0])
length(df_data_reads$TRB[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRB != 0])

sum(df_data_reads$TRD[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRD != 0])
median(df_data_reads$TRD[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRD != 0])
range(df_data_reads$TRD[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRD != 0])
length(df_data_reads$TRD[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRD != 0])

sum(df_data_reads$TRD[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRD != 0])
median(df_data_reads$TRD[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRD != 0])
range(df_data_reads$TRD[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRD != 0])
length(df_data_reads$TRD[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRD != 0])

sum(df_data_reads$TRG[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRG != 0])
median(df_data_reads$TRG[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRG != 0])
range(df_data_reads$TRG[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRG != 0])
length(df_data_reads$TRG[df_data_reads$steroid == "Steroid_High" & df_data_reads$TRG != 0])

sum(df_data_reads$TRG[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRG != 0])
median(df_data_reads$TRG[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRG != 0])
range(df_data_reads$TRG[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRG != 0])
length(df_data_reads$TRG[df_data_reads$steroid == "Steroid_Low" & df_data_reads$TRG != 0])




sum(df_data_clones$IGH[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGH != 0])
median(df_data_clones$IGH[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGH != 0])
range(df_data_clones$IGH[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGH != 0])
length(df_data_clones$IGH[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGH != 0])

sum(df_data_clones$IGH[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGH != 0])
median(df_data_clones$IGH[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGH != 0])
range(df_data_clones$IGH[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGH != 0])
length(df_data_clones$IGH[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGH != 0])

sum(df_data_clones$IGK[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGK != 0])
median(df_data_clones$IGK[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGK != 0])
range(df_data_clones$IGK[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGK != 0])
length(df_data_clones$IGK[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGK != 0])

sum(df_data_clones$IGK[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGK != 0])
median(df_data_clones$IGK[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGK != 0])
range(df_data_clones$IGK[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGK != 0])
length(df_data_clones$IGK[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGK != 0])

sum(df_data_clones$IGL[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGL != 0])
median(df_data_clones$IGL[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGL != 0])
range(df_data_clones$IGL[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGL != 0])
length(df_data_clones$IGL[df_data_clones$steroid == "Steroid_High" & df_data_clones$IGL != 0])

sum(df_data_clones$IGL[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGL != 0])
median(df_data_clones$IGL[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGL != 0])
range(df_data_clones$IGL[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGL != 0])
length(df_data_clones$IGL[df_data_clones$steroid == "Steroid_Low" & df_data_clones$IGL != 0])

sum(df_data_clones$TRA[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRA != 0])
median(df_data_clones$TRA[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRA != 0])
range(df_data_clones$TRA[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRA != 0])
length(df_data_clones$TRA[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRA != 0])

sum(df_data_clones$TRA[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRA != 0])
median(df_data_clones$TRA[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRA != 0])
range(df_data_clones$TRA[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRA != 0])
length(df_data_clones$TRA[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRA != 0])

sum(df_data_clones$TRB[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRB != 0])
median(df_data_clones$TRB[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRB != 0])
range(df_data_clones$TRB[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRB != 0])
length(df_data_clones$TRB[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRB != 0])

sum(df_data_clones$TRB[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRB != 0])
median(df_data_clones$TRB[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRB != 0])
range(df_data_clones$TRB[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRB != 0])
length(df_data_clones$TRB[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRB != 0])

sum(df_data_clones$TRD[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRD != 0])
median(df_data_clones$TRD[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRD != 0])
range(df_data_clones$TRD[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRD != 0])
length(df_data_clones$TRD[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRD != 0])

sum(df_data_clones$TRD[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRD != 0])
median(df_data_clones$TRD[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRD != 0])
range(df_data_clones$TRD[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRD != 0])
length(df_data_clones$TRD[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRD != 0])

sum(df_data_clones$TRG[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRG != 0])
median(df_data_clones$TRG[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRG != 0])
range(df_data_clones$TRG[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRG != 0])
length(df_data_clones$TRG[df_data_clones$steroid == "Steroid_High" & df_data_clones$TRG != 0])

sum(df_data_clones$TRG[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRG != 0])
median(df_data_clones$TRG[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRG != 0])
range(df_data_clones$TRG[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRG != 0])
length(df_data_clones$TRG[df_data_clones$steroid == "Steroid_Low" & df_data_clones$TRG != 0])

#write.csv(df_data_reads, file = "df_data_reads.csv")
#write.csv(df_data_clones, file = "df_data_clones.csv")
