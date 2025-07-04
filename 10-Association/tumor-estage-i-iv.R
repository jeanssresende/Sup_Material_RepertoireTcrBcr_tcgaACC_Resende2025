# Carregar pacotes necessários
library(dplyr)
library(ggplot2)
library(tidyr)      # Para pivot_longer
library(ggpubr)     # Para stat_compare_means (útil para o rótulo)
library(rstatix)    # Para adjust_pvalue e estatísticas

data_clones <- read.csv("../6-dataExploration/df_data_clones.csv")

data <- data_clones[,c("IGH","IGK","IGL","TRA","TRB","TRG","TRD",
                       "steroid","tumor_stage")]

rownames(data) <- data_clones$TCGA_barcode
head(data)









dados_cadeias <- data


# Definir as colunas das cadeias BCR e TCR
bcr_chains <- c("IGH", "IGK", "IGL")
tcr_chains <- c("TRA", "TRB", "TRG", "TRD")
all_chains_cols <- c(bcr_chains, tcr_chains)

# Converter para formato longo
data_long <- dados_cadeias %>%
  pivot_longer(
    cols = all_chains_cols,
    names_to = "chain",
    values_to = "abundance"
  ) %>%
  mutate(
    chain_type = case_when(
      chain %in% bcr_chains ~ "BCR",
      chain %in% tcr_chains ~ "TCR",
      TRUE ~ NA_character_ # Para qualquer outra cadeia que possa aparecer
    ),
    # Adicionar 1 e aplicar log10 para melhor visualização
    abundance_log10 = log10(abundance + 1)
  )

# Tratar e reordenar 'tumor_stage'
data_long$tumor_stage <- factor(
  data_long$tumor_stage,
  levels = c("stage i", "stage ii", "stage iii", "stage iv", "not reported") # Defina a ordem desejada
)

# Filtrar `not reported` se não for para incluir nas comparações de estágio, mas sim no plot
# Se você quer incluir "not reported" nas comparações (mesmo que com NA no p-valor), mantenha.
# Para este exemplo, vou manter, mas você pode filtrar se quiser apenas estágios numéricos.
# data_long_filtered <- data_long %>% filter(tumor_stage != "not reported")


# --- 2. Função para gerar Boxplot e informações do Teste ---
# Esta função AGORA NÃO SALVA O GRÁFICO DIRETAMENTE, nem adiciona o p-valor
# Ela retorna o objeto ggplot e as informações do teste (p-valor bruto, método)
gerar_boxplot_com_teste <- function(data_subset_for_plot, current_chain_name, current_tumor_stage_name, variavel_agrupadora) {
  
  # Nome da variável de abundância transformada (sempre abundance_log10)
  y_var_name <- "abundance_log10" 
  
  # Validar se o subset de dados é suficiente para plotagem e teste
  if (nrow(data_subset_for_plot) == 0) {
    warning(paste0("No data for Chain: ", current_chain_name, " - Stage: ", current_tumor_stage_name, ". Skipping."))
    return(list(plot = NULL, p_value_bruto = NA, test_method = "N/A", error_msg = "No data"))
  }
  
  # Obter os nomes dos grupos de esteroides presentes neste subset
  groups_present <- unique(data_subset_for_plot[[variavel_agrupadora]])
  
  if (length(groups_present) < 2) {
    warning(paste0("Only one or no steroid group in Chain: ", current_chain_name, " - Stage: ", current_tumor_stage_name, ". Skipping comparison."))
    return(list(plot = NULL, p_value_bruto = NA, test_method = "N/A", error_msg = "Single group"))
  }
  
  # Extrair os dados para cada grupo de esteroides
  data_group1 <- data_subset_for_plot %>% filter(.data[[variavel_agrupadora]] == groups_present[1]) %>% pull(!!y_var_name)
  data_group2 <- data_subset_for_plot %>% filter(.data[[variavel_agrupadora]] == groups_present[2]) %>% pull(!!y_var_name)
  
  # --- Realizar o Teste de Mann-Whitney U (Wilcoxon) ---
  # Sempre usaremos Wilcoxon conforme sua preferência, sem teste de normalidade prévio.
  current_p_value_bruto <- NA
  metodo_teste <- "wilcox.test"
  
  # Usar tryCatch para lidar com casos onde o teste não pode ser computado (ex: todos os valores idênticos, muito poucos dados)
  if (length(data_group1) >= 2 && length(data_group2) >= 2) { # Wilcoxon precisa de pelo menos 2 observações por grupo
    test_res <- tryCatch({
      wilcox.test(data_group1, data_group2)
    }, error = function(e) {
      warning(paste0("Erro no Wilcoxon para Chain: ", current_chain_name, " - Stage: ", current_tumor_stage_name, " (", groups_present[1], " vs ", groups_present[2], "): ", e$message))
      return(list(p.value = NA)) # Retorna NA se houver erro
    }, finally = {
      # Captura o warning 'não é possível computar o valor de p exato com o de desempate'
      # Isso é normal e não impede o cálculo da aproximação
    })
    current_p_value_bruto <- test_res$p.value
  } else {
    warning(paste0("Dados insuficientes para Wilcoxon (min 2 obs por grupo) para Chain: ", current_chain_name, " - Stage: ", current_tumor_stage_name, ". Retornando NA para p-valor bruto."))
  }
  
  # --- Criar o Boxplot ---
  # Nomes dos grupos para o caption
  n_grupo1 <- length(data_group1)
  n_grupo2 <- length(data_group2)
  
  boxplot_obj <- ggplot(data_subset_for_plot, aes(x = .data[[variavel_agrupadora]], y = .data[[y_var_name]], fill = .data[[variavel_agrupadora]])) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_point(aes(color = .data[[variavel_agrupadora]]), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1), alpha = 0.6, shape = 16) +
    # Cores fixas para esteroides para consistência em todos os plots
    scale_fill_manual(values = c("Steroid_Low" = "lightblue", "Steroid_High" = "salmon")) +
    scale_color_manual(values = c("Steroid_Low" = "darkblue", "Steroid_High" = "darkred")) +
    labs(
      title = paste0("log10(Abundância + 1) de ", current_chain_name, "\nStage: ", current_tumor_stage_name),
      x = "Status Esteroides",
      y = "log10(Abundância + 1)",
      fill = "Status Esteroides",
      color = "Status Esteroides",
      caption = paste0("N(", groups_present[1], ") = ", n_grupo1, ", N(", groups_present[2], ") = ", n_grupo2, "\nTeste: ", metodo_teste)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.x = element_blank(), # Esconde o rótulo do eixo X se já está no título
      legend.position = "none" # A legenda de fill/color já está no facet grid
    )
  
  return(list(plot = boxplot_obj, p_value_bruto = current_p_value_bruto, test_method = metodo_teste, error_msg = NA_character_))
}


# --- 3. Execução Principal: Coletar, Ajustar P-valores, Plotar ---

# Preparar uma lista para armazenar todos os resultados de plotagem (plot, p-valor bruto, etc.)
all_individual_plot_results <- list()

# Iterar sobre cada cadeia e cada estágio para gerar os boxplots e p-valores brutos
for (chain_name in all_chains_cols) { # Itera por IGH, IGK, IGL, TRA, TRB, TRG, TRD
  for (stage_name in levels(data_long$tumor_stage)) { # Itera por stage i, ii, iii, iv, not reported
    
    # Filtrar os dados para a combinação atual de cadeia e estágio
    # Certifique-se de filtrar pelo nome da cadeia e estágio
    df_subset_current <- data_long %>%
      filter(chain == chain_name, tumor_stage == stage_name)
    
    # Chamar a função para gerar o boxplot e obter o p-valor bruto
    plot_key <- paste0(chain_name, "_", stage_name)
    
    plot_info <- gerar_boxplot_com_teste(
      data_subset_for_plot = df_subset_current,
      current_chain_name = chain_name,
      current_tumor_stage_name = stage_name,
      variavel_agrupadora = "steroid"
    )
    
    # Armazenar o resultado, se não for NULL (para casos de dados insuficientes)
    if (!is.null(plot_info$plot)) {
      all_individual_plot_results[[plot_key]] <- list(
        plot_obj = plot_info$plot,
        p_value_bruto = plot_info$p_value_bruto,
        test_method = plot_info$test_method,
        chain = chain_name,
        tumor_stage = stage_name,
        error_msg = plot_info$error_msg # Armazena qualquer erro/aviso
      )
    } else {
      message(paste0("Skipping plot for ", chain_name, " - Stage: ", stage_name, ": ", plot_info$error_msg))
    }
  }
}

# --- Coletar todos os p-valores brutos para ajuste ---
p_values_to_adjust <- data.frame(
  chain = character(),
  tumor_stage = character(),
  p_value = numeric(),
  test_method = character(),
  stringsAsFactors = FALSE
)

for (plot_key in names(all_individual_plot_results)) {
  info <- all_individual_plot_results[[plot_key]]
  if (!is.na(info$p_value_bruto)) { # Apenas p-valores válidos
    p_values_to_adjust <- bind_rows(p_values_to_adjust, data.frame(
      chain = info$chain,
      tumor_stage = info$tumor_stage,
      p_value = info$p_value_bruto,
      test_method = info$test_method,
      stringsAsFactors = FALSE
    ))
  }
}

adjusted_p_values_table <- NULL
if (nrow(p_values_to_adjust) > 0) {
  adjusted_p_values_table <- p_values_to_adjust %>%
    rstatix::adjust_pvalue(method = "BH") # Correção Benjamini-Hochberg (FDR)
} else {
  message("Nenhum p-valor bruto válido para ajuste. Verifique seus dados de entrada.")
}

# --- Adicionar p-valores ajustados aos gráficos e salvar/exibir ---
if (!is.null(adjusted_p_values_table) && nrow(adjusted_p_values_table) > 0) {
  # Lista para armazenar os gráficos finais com p-valores adicionados
  final_plots_list <- list()
  
  for (plot_key in names(all_individual_plot_results)) {
    current_plot_info <- all_individual_plot_results[[plot_key]]
    current_plot_obj <- current_plot_info$plot_obj
    
    # Encontrar o p-valor ajustado correspondente
    adj_info_row <- adjusted_p_values_table %>%
      filter(chain == current_plot_info$chain,
             tumor_stage == current_plot_info$tumor_stage)
    
    p_adj_val_display <- NA # Valor padrão para exibição
    if (nrow(adj_info_row) > 0) {
      p_adj_val_display <- adj_info_row$p.adj[1]
    }
    
    # Formatar o p-valor para exibição (ex: <0.001)
    formatted_p_value <- format.pval(p_adj_val_display, digits = 3, eps = 0.001) # Mostra <0.001 para p muito pequenos
    
    # Adicionar o p-valor ajustado ao gráfico
    plot_build <- ggplot_build(current_plot_obj)
    y_max_coord <- tryCatch(max(plot_build$data[[1]]$y, na.rm = TRUE),
                            error = function(e) {
                              warning(paste("Could not determine y_max for", plot_key, ". Using default based on plot data."))
                              # Fallback to max of actual data
                              return(max(current_plot_obj$data$abundance_log10, na.rm = TRUE) * 1.1)
                            })
    if (is.infinite(y_max_coord) || is.na(y_max_coord) || y_max_coord == 0) y_max_coord <- 1 # Default if max is 0 or invalid
    
    plot_with_adj_p <- current_plot_obj +
      annotate("text",
               x = 1.5, # Posição X (meio entre os dois grupos de esteroides)
               y = y_max_coord * 1.15, # Posição Y (um pouco acima do boxplot)
               label = paste0("p.adj = ", formatted_p_value),
               color = "black", # Mudei para preto para contraste
               size = 3
      )
    
    # Adicionar o objeto plot_with_adj_p à lista de plots finais
    final_plots_list[[plot_key]] <- plot_with_adj_p
    
    # Salvar o gráfico individual em PDF
    filename_to_save <- paste0("boxplot_", current_plot_info$chain, "_stage_", current_plot_info$tumor_stage, ".pdf")
    ggsave(filename = filename_to_save, plot = plot_with_adj_p, width = 4, height = 4) # Tamanho menor para plots individuais
  }
  
  # --- Opcional: Combinar todos os plots em um único PDF grande (com patchwork) ---
  # Se você tiver muitos plots, pode ser melhor tê-los em PDFs separados
  # ou usar um pacote como patchwork para combiná-los em um grid.
  # install.packages("patchwork") # Se ainda não tiver
  # library(patchwork)
  
  # # Exemplo de como você poderia combinar os plots BCR
  # bcr_plots <- final_plots_list[grep("BCR", names(final_plots_list))]
  # # Organize os plots em um grid (ajuste ncol/nrow conforme necessário)
  # if (length(bcr_plots) > 0) {
  #   combined_bcr_plot <- wrap_plots(bcr_plots, ncol = 3) +
  #     plot_annotation(title = 'BCR Chains Abundance by Stage and Steroid Status (Adjusted p-values)') &
  #     theme(plot.title = element_text(hjust = 0.5))
  #   ggsave("combined_bcr_boxplots_adjusted_p.pdf", plot = combined_bcr_plot, width = 12, height = 8)
  # }
  
  # # Exemplo de como você poderia combinar os plots TCR
  # tcr_plots <- final_plots_list[grep("TCR", names(final_plots_list))]
  # if (length(tcr_plots) > 0) {
  #   combined_tcr_plot <- wrap_plots(tcr_plots, ncol = 4) +
  #     plot_annotation(title = 'TCR Chains Abundance by Stage and Steroid Status (Adjusted p-values)') &
  #     theme(plot.title = element_text(hjust = 0.5))
  #   ggsave("combined_tcr_boxplots_adjusted_p.pdf", plot = combined_tcr_plot, width = 16, height = 8)
  # }
  
  # Exibir a tabela de p-valores ajustados no final
  message("\n--- Tabela de P-valores Ajustados ---")
  print(adjusted_p_values_table)
  
} else {
  message("Não há p-valores ajustados para plotar. Verifique os dados de entrada ou se todos os testes falharam.")
}
