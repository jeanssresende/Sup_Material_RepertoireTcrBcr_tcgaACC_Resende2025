# -- calculo da distancia e agrupamento -- #####################################

# pacotes necessarios
library(immunarch)
library(dplyr)

# -- preparando o ambiente de trabalho -- ######################################
immdata <- repLoad("../3-extractionCDR3/cdr3/")
load("../4-IDs-metadata/metadata_2025.RData")

immdata$meta <- metadata[metadata$sample_id %in% 
                           gsub("_report", "", immdata$meta$Sample),]

# colune de ligacao do meta
metadata$Sample <- paste(metadata$file_id,"_report", sep = "")  

# metadata$Sample[metadata$Sample %in% names(immdata$data) == F] 
immdata$meta <- metadata[metadata$Sample %in% names(immdata$data),] 

# -- separacao dos repertorios de TCR e BCR -- #################################
immdata_bcr <- repFilter(immdata,
                         .method = "by.clonotype",
                         .query = list(V.name = include("IG")),
                         .match = "substring")

immdata_tcr <- repFilter(immdata,
                         .method = "by.clonotype",
                         .query = list(V.name = include("TR")),
                         .match = "substring")

# -- funcao para remover linhas com NA na coluna J.name -- #####################
remove_NA <- function(df){
  df %>%
    filter(!is.na(J.name))
}

# -- funcao para separar os NA e acrescentar a tabela final depois -- ##########
capturar_NA <- function(df){
  df %>%
    filter(is.na(J.name))
}

# -- NAs em J.name para acrescentar -- #########################################
TCR_NAs_acrescentar <- lapply(immdata_tcr$data, capturar_NA)
BCR_NAs_acrescentar <- lapply(immdata_bcr$data, capturar_NA)

# -- removendo linhas com valores NA na coluna J.name -- #######################
immdata_tcr$data <- lapply(immdata_tcr$data, remove_NA)
immdata_bcr$data <- lapply(immdata_bcr$data, remove_NA)

# -- calculo da distancia entre as sequencias -- ###############################
dist_BCR <- seqDist(immdata_bcr$data, 
                    .group_by = c("V.name", "J.name"),
                    .group_by_seqLength = TRUE,
                    .col = "CDR3.nt")

clust_BCR <- seqCluster(immdata_bcr$data, dist_BCR, .perc_similarity = 0.9)


dist_TCR <- seqDist(immdata_tcr$data[-12], 
                    .group_by = c("V.name", "J.name"),
                    .group_by_seqLength = TRUE,
                    .col = "CDR3.nt")

clust_TCR <- seqCluster(immdata_tcr$data[-12], dist_TCR, .perc_similarity = 0.95)


temp <- immdata_tcr$data[12]

temp$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`$J.name
nchar(temp$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`$CDR3.nt)

temp$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`$Cluster <- 
  c("TRAV9-2/TRAJ17_length_39","TRAV9-2/TRAJ17_length_39")

names(temp)

clust_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report` <- 
  as.data.frame(temp$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`) 


clust_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_TAGCTT_L006_report`

clust_TCR$`130723_UNC9-SN296_0386_BC2E4WACXX_ACTTGA_L003_report`

# -- funcao para capturar as linhas de cada dataframe que contem NAs em Cluster
capturar_NA_cluster <- function(df){
  df %>%
    filter(is.na(Cluster))
}

# -- NAs em cluster para acrescentar -- ########################################
TCR_NAs_acrescentar_cluster <- lapply(clust_TCR, capturar_NA_cluster)
BCR_NAs_acrescentar_cluster <- lapply(clust_BCR, capturar_NA_cluster)

# -- funcao para remover linhas com NA na coluna Cluster -- ####################
remove_NA_cluster <- function(df){
  df %>%
    filter(!is.na(Cluster))
}

# -- removendo NAs para acrescentar depois -- ##################################
clust_TCR <- lapply(clust_TCR, remove_NA_cluster)
clust_BCR <- lapply(clust_BCR, remove_NA_cluster)

# -- funcao para agrupar os dataframes -- ######################################
agrupar_por_cluster <- function(df){
  df %>%
    group_by(Cluster) %>%
    mutate(Clones = sum(Clones),
           Proportion = sum(Proportion)) %>%
    ungroup() %>%
    distinct(Cluster, .keep_all = TRUE) %>%
    as.data.frame()
}

# -- aplicando a funcao na lista de dataframes -- ##############################
clust_TCR_agrupados <- lapply(clust_TCR, agrupar_por_cluster)
clust_BCR_agrupados <- lapply(clust_BCR, agrupar_por_cluster)

# -- adicionando coluna cluster -- #############################################
add_col <- function(lista){
  for (i in seq_along(lista)) {
    lista[[i]]$Cluster <- NA
  }
  return(lista)
}

BCR_NAs_acrescentar <- add_col(BCR_NAs_acrescentar)
TCR_NAs_acrescentar <- add_col(TCR_NAs_acrescentar)

juntar_df <- function(df1, df2, df3){
  return(rbind(df1,df2,df3))
}

clust_TCR_agrupados <- Map(juntar_df, clust_TCR_agrupados, 
                           TCR_NAs_acrescentar, TCR_NAs_acrescentar_cluster)

clust_BCR_agrupados <- Map(juntar_df, clust_BCR_agrupados,
                           BCR_NAs_acrescentar, BCR_NAs_acrescentar_cluster)

# -- funcao para juntar listas -- ##############################################
juntar_listas <- function(lista1, lista2){
  nomes_comuns <- intersect(names(lista1), names(lista2))
  lista_final <- list()
  
  for (nome in nomes_comuns) {
    df1 <- lista1[[nome]]
    df2 <- lista2[[nome]]
    
    lista_final[[nome]] <- rbind(df1,df2)
  }
  
  nomes_lista1 <- setdiff(names(lista1), nomes_comuns)
  for (nome in nomes_lista1) {
    lista_final[[nome]] <- lista1[[nome]]
  }
  
  nomes_lista2 <- setdiff(names(lista2), nomes_comuns)
  for (nome in nomes_lista2) {
    lista_final[[nome]] <- lista2[[nome]]
  }
  
  return(lista_final)
}

lista_tcr_bcr <- juntar_listas(clust_TCR_agrupados, clust_BCR_agrupados)

lista_tcr_bcr$`130723_UNC9-SN296_0386_BC2E4WACXX_ACTTGA_L003_report`

# -- funcao para substituir NA por "." -- ######################################
substituir_na <- function(df){
  df[is.na(df)] <- "."
  return(df)
}

for (i in seq_along(lista_tcr_bcr)) {
  lista_tcr_bcr[[i]] <- substituir_na(lista_tcr_bcr[[i]])
}

lista_tcr_bcr$`130723_UNC9-SN296_0386_BC2E4WACXX_ACTTGA_L003_report`

# -- funcao para processar o dataframe -- ######################################
processar_dataframe <- function(df){
  
  # remover as colunas
  df <- df %>%
    select(-c(V.end, D.start, D.end, J.start, VJ.ins, VD.ins, DJ.ins, Cluster))
  
  # alterar o nome das colunas
  names(df) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J")
  
  # adicionar a coluna "C"
  df$C <- ifelse(substr(df$V, 1,3) == "IGL","IGLC",
                 ifelse(substr(df$V, 1,3) == "IGK","IGKC",
                        ifelse(substr(df$V, 1,3) == "IGH","IGHC",
                               ifelse(substr(df$V, 1,3) == "TRA","TRAC",
                                      ifelse(substr(df$V, 1,3) == "TRB","TRBC",
                                             ifelse(substr(df$V, 1,3) == "TRG",
                                                    "TRGC","TRDC"))))))
  return(df)
}

# -- funcao para processar cada dataframe -- ###################################
processar_lista <- function(lista){
  nomes <- names(lista)
  for (i in seq_along(lista)) {
    lista[[i]] <- processar_dataframe(lista[[i]])
  }
  names(lista) <- nomes
  return(lista)
}

# -- aplicando as operacoes na lista de dataframes -- ##########################
lista_tcr_bcr_formatada <- processar_lista(lista_tcr_bcr)

lista_tcr_bcr_formatada$`130723_UNC9-SN296_0386_BC2E4WACXX_ACTTGA_L003_report`

# -- funcao para atualizar os valores da coluna frequency -- ###################
atualizar_frequency<- function(df){
  df %>%
    mutate(
      frequency = ifelse(substr(V,1,2) == "IG",
                         `count` / sum(filter(df, substr(V,1,2) == "IG")$`count`),
                         `count` / sum(filter(df, substr(V,1,2) == "TR")$`count`))
    )
}

# --funcao para atualizar cada dataframe na lista -- ###########################
atualizar_lista <- function(lista){
  for (nome_df in names(lista)) {
    lista[[nome_df]] <- atualizar_frequency(lista[[nome_df]])
  }
  
  return(lista)
}

# -- aplica as operacoes na lista de dataframes -- #############################
lista_tcr_bcr_formatada_frequency <- atualizar_lista(lista_tcr_bcr_formatada)

# - funcao para salvar cada dataframe individualmente no formato desejado
salvar_dataframes_formatado <- function(lista, dir_destino){
  for (nome_df in names(lista)) {
    arquivo <- paste0(dir_destino, "/", gsub("_report", "", nome_df), "_clones.tsv")
    df_para_salvar <- lista[[nome_df]]
    colnames(df_para_salvar) <- c("#count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "C")
    write.table(df_para_salvar, arquivo, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

# - salvar os dataframes como .tsv no formato desejado
salvar_dataframes_formatado(lista_tcr_bcr_formatada_frequency, dir_destino = "data")

save(lista_tcr_bcr_formatada_frequency,
     file = "lista_tcrbcr_formatadaFrequency.RData")
