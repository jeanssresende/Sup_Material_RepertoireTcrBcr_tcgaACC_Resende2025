library(UCSCXenaTools)
library(tidyverse)
library(IOBR)

eset_acc <- read.table("TCGA-ACC.star_counts.tsv", header = T)

## -- anotacao e normalizacao em TPM

# Remover a versão do Ensembl ID
eset_acc$Ensembl_ID <- sub("\\.[0-9]+$", "", eset_acc$Ensembl_ID)

# Agrupar pelo Ensembl_ID e calcular a média
eset_acc_agg <- eset_acc %>%
  group_by(Ensembl_ID) %>%
  summarise(across(everything(), median))

# Transformar Ensembl_ID em nomes de linha
eset_acc_final <- column_to_rownames(eset_acc_agg, var = "Ensembl_ID")

eset_acc_final <- (2^eset_acc_final)+1 # revertendo para a contagem original (pois estao em log2)
eset_acc_final <- count2tpm(countMat = eset_acc_final, idType = "Ensembl", org = "hsa") # normalizacao em TPM

epic <- deconvo_tme(eset = log2(eset_acc_final+1), method = "epic", arrays = FALSE)

xcell <- deconvo_tme(eset = log2(eset_acc_final+1), method = "xcell", arrays = FALSE)

mcpcounter <- deconvo_tme(eset = log2(eset_acc_final+1), method = "mcpcounter", arrays = FALSE)

estimate <- deconvo_tme(eset = log2(eset_acc_final+1), method = "estimate", arrays = FALSE)

timer <- deconvo_tme(eset = log2(eset_acc_final+1),
                     method = "timer",
                     group_list = rep("acc",ncol(eset_acc_final)),
                     arrays = FALSE)

ips <- deconvo_tme(eset = log2(eset_acc_final+1), method = "ips", plot = FALSE)

cibersort <- deconvo_tme(eset = log2(eset_acc_final+1),
                         method = "cibersort",
                         arrays = FALSE,
                         perm = 200)

quantiseq <- deconvo_tme(eset = log2(eset_acc_final+1),
                         tumor = TRUE,
                         arrays = FALSE,
                         scale_mrna = TRUE,
                         method = "quantiseq")

write.csv(epic, file = "epic.csv")
write.csv(xcell, file = "xcell.csv")
write.csv(mcpcounter, file = "mcpcounter.csv")
write.csv(estimate, file = "estimate.csv")
write.csv(timer, file = "timer.csv")
write.csv(ips, file = "ips.csv")
write.csv(cibersort, file = "cibersort.csv")
write.csv(quantiseq, file = "quantiseq.csv")

colnames(quantiseq)

epic_tcrbcr <- epic[,c("Bcells_EPIC","CD4_Tcells_EPIC","CD8_Tcells_EPIC")]

xcell_tcrbcr <- xcell[,c("B-cells_xCell","Memory_B-cells_xCell",
                         "naive_B-cells_xCell","pro_B-cells_xCell",
                         "CD4+_memory_T-cells_xCell","CD4+_naive_T-cells_xCell",
                         "CD4+_T-cells_xCell","CD4+_Tcm_xCell","CD4+_Tem_xCell",
                         "CD8+_naive_T-cells_xCell","CD8+_T-cells_xCell",
                         "CD8+_Tcm_xCell","CD8+_Tem_xCell","Th1_cells_xCell",
                         "Th2_cells_xCell","Tregs_xCell")]

mcpcounter_tcrbcr <- mcpcounter[,c("T_cells_MCPcounter","CD8_T_cells_MCPcounter",
                                   "B_lineage_MCPcounter","Cytotoxic_lymphocytes_MCPcounter")]

timer_tcrbcr <- timer[,c("B_cell_TIMER","T_cell_CD4_TIMER","T_cell_CD8_TIMER")]

cibersort_tcrbcr <- cibersort[,c("B_cells_naive_CIBERSORT","B_cells_memory_CIBERSORT",
                                 "T_cells_CD8_CIBERSORT","T_cells_CD4_naive_CIBERSORT",
                                 "T_cells_CD4_memory_resting_CIBERSORT",
                                 "T_cells_CD4_memory_activated_CIBERSORT",
                                 "T_cells_follicular_helper_CIBERSORT",
                                 "T_cells_regulatory_(Tregs)_CIBERSORT",
                                 "T_cells_gamma_delta_CIBERSORT")]

quantiseq_tcrbcr <- quantiseq[,c("B_cells_quantiseq","T_cells_CD4_quantiseq",
                                 "T_cells_CD8_quantiseq","Tregs_quantiseq")]

metricas_deconv_tcrbcr <- cbind(epic_tcrbcr, xcell_tcrbcr, mcpcounter_tcrbcr,
                                timer_tcrbcr, cibersort_tcrbcr, quantiseq_tcrbcr)

rownames(metricas_deconv_tcrbcr) <- epic$ID

write.csv(metricas_deconv_tcrbcr, file = "metricas_deconv_tcrbcr.csv")
