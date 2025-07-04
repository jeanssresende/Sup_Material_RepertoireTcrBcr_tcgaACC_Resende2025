library(ggplot2)

data_clones <- read.csv("../../6-dataExploration/df_data_clones.csv")
mcpcounter <- read.csv("mcpcounter.csv")

data <- mcpcounter

data$ID <- gsub("//.","-", data$ID)

data_clones$tcga_cod <- gsub("-",".",data_clones$TCGA_barcode) 
data_clones$tcga_cod <- substr(data_clones$tcga_cod,1,16)

data <- data[data$ID %in% data_clones$tcga_cod,
             c("ID","T_cells_MCPcounter","CD8_T_cells_MCPcounter",
               "Cytotoxic_lymphocytes_MCPcounter", "B_lineage_MCPcounter")]

idx <- match(data$ID, data_clones$tcga_cod)

data$Steroid_phenotype <- data_clones$steroid[idx]

shapiro.test(data$T_cells_MCPcounter)
shapiro.test(data$CD8_T_cells_MCPcounter)
shapiro.test(data$Cytotoxic_lymphocytes_MCPcounter)
shapiro.test(data$B_lineage_MCPcounter)

p_T_cells_MCPcounter <- wilcox.test(T_cells_MCPcounter ~ Steroid_phenotype, data = data)
p_CD8_T_cells_MCPcounter <- wilcox.test(CD8_T_cells_MCPcounter ~ Steroid_phenotype, data = data)
p_Cytotoxic_lymphocytes_MCPcounter <- wilcox.test(Cytotoxic_lymphocytes_MCPcounter ~ Steroid_phenotype, data = data)
p_B_lineage_MCPcounter <- wilcox.test(B_lineage_MCPcounter ~ Steroid_phenotype, data = data)


ggplot(data, aes(x = Steroid_phenotype, y = log2(T_cells_MCPcounter + 1), fill = Steroid_phenotype)) +
  geom_boxplot() +
  labs(title = "Boxplot de log2(T_cells_MCPcounter + 1)",
       x = "Steroid",
       y = "log2(Score + 1)") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(log2(data$T_cells_MCPcounter + 1)) * 1.1,
           label = paste("p =", format.pval(p_T_cells_MCPcounter$p.value, digits = 3)))

ggplot(data, aes(x = Steroid_phenotype, y = log2(CD8_T_cells_MCPcounter + 1), fill = Steroid_phenotype)) +
  geom_boxplot() +
  labs(title = "Boxplot de log2(CD8_T_cells_MCPcounter + 1)",
       x = "Steroid",
       y = "log2(Score + 1)") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(log2(data$CD8_T_cells_MCPcounter + 1)) * 1.1,
           label = paste("p =", format.pval(p_CD8_T_cells_MCPcounter$p.value, digits = 3)))

ggplot(data, aes(x = Steroid_phenotype, y = log2(Cytotoxic_lymphocytes_MCPcounter + 1), fill = Steroid_phenotype)) +
  geom_boxplot() +
  labs(title = "Boxplot de log2(Cytotoxic_lymphocytes_MCPcounter + 1)",
       x = "Steroid",
       y = "log2(Score + 1)") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(log2(data$Cytotoxic_lymphocytes_MCPcounter + 1)) * 1.1,
           label = paste("p =", format.pval(p_Cytotoxic_lymphocytes_MCPcounter$p.value, digits = 3)))

ggplot(data, aes(x = Steroid_phenotype, y = log2(B_lineage_MCPcounter + 1), fill = Steroid_phenotype)) +
  geom_boxplot() +
  labs(title = "Boxplot de log2(B_lineage_MCPcounter + 1)",
       x = "Steroid",
       y = "log2(Score + 1)") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(log2(data$B_lineage_MCPcounter + 1)) * 1.1,
           label = paste("p =", format.pval(p_B_lineage_MCPcounter$p.value, digits = 3)))



