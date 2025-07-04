library(jsonlite)

# Carregue o JSON SEM simplificar vetores automaticamente (ajuste o nome do arquivo)
metadata_list <- fromJSON("../1-obtainingData/metadata.cart.2022-04-12.json", simplifyVector = FALSE) #

# Extrair o mapeamento 
mapping_list <- lapply(metadata_list, function(entry) {
  # Adicionar uma verificação para garantir que 'entry' é uma lista e tem os campos esperados
  if (is.list(entry) && !is.null(entry$file_id) && !is.null(entry$cases) && 
      length(entry$cases) > 0 && !is.null(entry$cases[[1]]$submitter_id)) {
    data.frame(
      file_id = entry$file_id,
      file_name = entry$file_name,
      TCGA_Case_ID = entry$cases[[1]]$submitter_id, # Usar [[1]] pois 'cases' é uma lista dentro de 'entry'
      # Posso adicionar outros campos aqui, se necessário:
      TCGA_Sample_ID = entry$cases[[1]]$samples[[1]]$submitter_id,
      TCGA_Aliquot_ID = entry$associated_entities[[1]]$entity_submitter_id
    )
  } else {
    # Retorna NULL se a estrutura não for a esperada, para evitar erros
    NULL
  }
})

# Remover entradas NULL que podem ter sido retornadas pela verificação
mapping_list <- Filter(Negate(is.null), mapping_list)

# Combinar a lista de data frames em um único data frame
if (length(mapping_list) > 0) {
  mapping_df <- do.call(rbind, mapping_list)
  
  # Remover duplicatas, se houver
  mapping_df <- unique(mapping_df)
  
  # Exibir as primeiras linhas para verificação
  print(head(mapping_df))
  
} else {
  print("Nenhuma informação válida de mapeamento foi extraída do arquivo JSON.")
  mapping_df <- NULL # Define como NULL se nada foi extraído
}

# regex

sample_names <- gsub("_report.tsv","", list.files("../3-extractionCDR3/cdr3/")) 

mapping_df$file_name[1]

sample_names_meta <- sub("^[^.]*\\.[^.]*\\.", "", mapping_df$file_name)
sample_names_meta <- sub("\\.tar\\.gz$", "", sample_names_meta)

mapping_df$sample_names_meta <- sample_names_meta

# 
load("coldataACC.RData")

metadata <- data.frame(
  TCGA_file_id = mapping_df$file_id,
  TCGA_sample_names_meta = mapping_df$sample_names_meta,
  TCGA_file_name = mapping_df$file_name,
  TCGA_barcode = mapping_df$TCGA_Aliquot_ID)

head(metadata)

idx <- match(metadata$TCGA_barcode, coldataACC$barcode)

metadata$steroid <- coldataACC$steroid[idx]
metadata$cortisol.excess <- coldataACC$cortisol.excess[idx]
metadata$immune.subtype <- coldataACC$immune.subtype[idx]
metadata$other.hormones <- coldataACC$other.hormones[idx]
metadata$tumor_stage <- coldataACC$tumor_stage[idx]
metadata$vital_status <- coldataACC$vital_status[idx]
metadata$gender <- coldataACC$gender[idx]

head(metadata)

idx <- match(metadata$TCGA_barcode, coldataACC$barcode)

metadata$file_id <- coldataACC$sample_id[idx]

metadata$file_id[is.na(metadata$file_id)] <- metadata$TCGA_file_id[is.na(metadata$file_id)] 

metadata$file_id[metadata$file_id == "fe298483-c42c-4753-babb-29a9cb1b6312"] <-
  "94a6881c-34ca-4f00-8fac-640f97a56c01"

metadata$file_id[metadata$file_id == "ae4bfe63-f0c0-4624-ad97-4f7ff13e2384"] <-
  "126e4c65-dd29-4040-a256-8d5f2feac3dc"

metadata$file_id[metadata$file_id == "ff34117b-5ee4-4c96-a8a2-2ce6a654cbe5"] <-
  "06cc0bc6-1b36-4318-99ee-dfc83cc134c0"

save(metadata, file = "metadata_2025.RData")


