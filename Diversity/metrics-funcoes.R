# reports
outputTrust4_report <- "../../5-clonesIdentification/data/"
files <- list.files(outputTrust4_report)
writeLines(files, "samplesNames.txt") # gerando samplesNames.txt
file.exists("samplesNames.txt")

# executar o metrics-trust4.sh

## -- sequencias brutas
dir_sequencias <- "../../5-clonesIdentification/data/"

arquivos_sequencias <- list.files(dir_sequencias)

tcrbcr_sequencias <- list()

# loop para ler cada arquivo e armazenar os dados na lista
for (arquivo in arquivos_sequencias) {
  nome <- gsub("\\.tsv$", "", arquivo)
  dados <- read.delim(file.path(dir_sequencias , arquivo), sep = "\t")
  tcrbcr_sequencias[[nome]] <- dados
}

# abundance
calc.abundance <- function(count){
  return(sum(count))
}

# richness
calc.richness <- function(count){
  return(length(count))
}

# cpk
calc.cpk <- function(count){
  return(length(count) / sum(count) * 1000)
}

# entropy
calc.entropy <- function(count){
  j <- 0
  for (i in seq_along(count)) {
    t <- (-count[i])/sum(count) * log(count[i]/sum(count))
    j <- j + t
  }
  return(j)
}

# clonality
calc.clonality <- function(count){
  return(1 - calc.entropy(count) / log(length(count)))
}

calc_tcrbcr_sequencias <- list()

for (nome in names(tcrbcr_sequencias)) {
  
  calc_tcrbcr_sequencias[[nome]] <- data.frame(
    chain = c("IGH","IGK","IGL","TRA","TRB","TRG","TRD"),
    Abundance = rep(NA, 7),
    Richness = rep(NA, 7),
    CPK = rep(NA, 7),
    Entropy = rep(NA, 7),
    Clonality = rep(NA, 7)
  )
} 

for (nome in names(tcrbcr_sequencias)) {
  df <- tcrbcr_sequencias[[nome]]
  df2 <- calc_tcrbcr_sequencias[[nome]]
  
  # abundance
  df2[df2$chain == "IGH", "Abundance"] <- 
    calc.abundance(df$X.count[substr(df$V, 1,3) == "IGH"])
  
  df2[df2$chain == "IGK", "Abundance"] <- 
    calc.abundance(df$X.count[substr(df$V, 1,3) == "IGK"])
  
  df2[df2$chain == "IGL", "Abundance"] <- 
    calc.abundance(df$X.count[substr(df$V, 1,3) == "IGL"])
  
  df2[df2$chain == "TRA", "Abundance"] <- 
    calc.abundance(df$X.count[substr(df$V, 1,3) == "TRA"])
  
  df2[df2$chain == "TRB", "Abundance"] <- 
    calc.abundance(df$X.count[substr(df$V, 1,3) == "TRB"])
  
  df2[df2$chain == "TRG", "Abundance"] <- 
    calc.abundance(df$X.count[substr(df$V, 1,3) == "TRG"])
  
  df2[df2$chain == "TRD", "Abundance"] <- 
    calc.abundance(df$X.count[substr(df$V, 1,3) == "TRD"])
  
  
  # richness
  df2[df2$chain == "IGH", "Richness"] <-
    calc.richness(df$X.count[substr(df$V, 1,3) == "IGH"])
  
  df2[df2$chain == "IGK", "Richness"] <-
    calc.richness(df$X.count[substr(df$V, 1,3) == "IGK"])
  
  df2[df2$chain == "IGL", "Richness"] <-
    calc.richness(df$X.count[substr(df$V, 1,3) == "IGL"])
  
  df2[df2$chain == "TRA", "Richness"] <-
    calc.richness(df$X.count[substr(df$V, 1,3) == "TRA"])
  
  df2[df2$chain == "TRB", "Richness"] <-
    calc.richness(df$X.count[substr(df$V, 1,3) == "TRB"])
  
  df2[df2$chain == "TRG", "Richness"] <-
    calc.richness(df$X.count[substr(df$V, 1,3) == "TRG"])
  
  df2[df2$chain == "TRD", "Richness"] <-
    calc.richness(df$X.count[substr(df$V, 1,3) == "TRD"])
  
  # cpk
  df2[df2$chain == "IGH", "CPK"] <-
    calc.cpk(df$X.count[substr(df$V, 1,3) == "IGH"])
  
  df2[df2$chain == "IGK", "CPK"] <-
    calc.cpk(df$X.count[substr(df$V, 1,3) == "IGK"])
  
  df2[df2$chain == "IGL", "CPK"] <-
    calc.cpk(df$X.count[substr(df$V, 1,3) == "IGL"])
  
  df2[df2$chain == "TRA", "CPK"] <-
    calc.cpk(df$X.count[substr(df$V, 1,3) == "TRA"])
  
  df2[df2$chain == "TRB", "CPK"] <-
    calc.cpk(df$X.count[substr(df$V, 1,3) == "TRB"])
  
  df2[df2$chain == "TRG", "CPK"] <-
    calc.cpk(df$X.count[substr(df$V, 1,3) == "TRG"])
  
  df2[df2$chain == "TRD", "CPK"] <-
    calc.cpk(df$X.count[substr(df$V, 1,3) == "TRD"])
  
  # entropy
  df2[df2$chain == "IGH", "Entropy"] <-
    calc.entropy(df$X.count[substr(df$V, 1,3) == "IGH"])
  
  df2[df2$chain == "IGK", "Entropy"] <-
    calc.entropy(df$X.count[substr(df$V, 1,3) == "IGK"])
  
  df2[df2$chain == "IGL", "Entropy"] <-
    calc.entropy(df$X.count[substr(df$V, 1,3) == "IGL"])
  
  df2[df2$chain == "TRA", "Entropy"] <-
    calc.entropy(df$X.count[substr(df$V, 1,3) == "TRA"])
  
  df2[df2$chain == "TRB", "Entropy"] <-
    calc.entropy(df$X.count[substr(df$V, 1,3) == "TRB"])
  
  df2[df2$chain == "TRG", "Entropy"] <-
    calc.entropy(df$X.count[substr(df$V, 1,3) == "TRG"])
  
  df2[df2$chain == "TRD", "Entropy"] <-
    calc.entropy(df$X.count[substr(df$V, 1,3) == "TRD"])
  
  # clonality
  df2[df2$chain == "IGH", "Clonality"] <-
    calc.clonality(df$X.count[substr(df$V, 1,3) == "IGH"])
  
  df2[df2$chain == "IGK", "Clonality"] <-
    calc.clonality(df$X.count[substr(df$V, 1,3) == "IGK"])
  
  df2[df2$chain == "IGL", "Clonality"] <-
    calc.clonality(df$X.count[substr(df$V, 1,3) == "IGL"])
  
  df2[df2$chain == "TRA", "Clonality"] <-
    calc.clonality(df$X.count[substr(df$V, 1,3) == "TRA"])
  
  df2[df2$chain == "TRB", "Clonality"] <-
    calc.clonality(df$X.count[substr(df$V, 1,3) == "TRB"])
  
  df2[df2$chain == "TRG", "Clonality"] <-
    calc.clonality(df$X.count[substr(df$V, 1,3) == "TRG"])
  
  df2[df2$chain == "TRD", "Clonality"] <-
    calc.clonality(df$X.count[substr(df$V, 1,3) == "TRD"])
  
  
  calc_tcrbcr_sequencias[[nome]] <- df2
}

calc_tcrbcr_sequencias$`06cc0bc6-1b36-4318-99ee-dfc83cc134c0_clones`

# funcao para salvar cada dataframe individualmente
salvar_dataframes <- function(lista, dir_destino){
  for (nome_df in names(lista)) {
    arquivo <- paste0(dir_destino,"/", nome_df, ".tsv")
    write.table(lista[[nome_df]], arquivo, sep = "\t", row.names = FALSE)
  }
}

# salva os dataframes como .tsv
salvar_dataframes(calc_tcrbcr_sequencias, dir_destino = "results_pipeline_report")

# salva as listas  
save(calc_tcrbcr_sequencias, file = "calc_tcrbcr_sequencias.RData")
