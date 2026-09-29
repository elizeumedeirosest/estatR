linhas <- readLines('C:/Users/Elize/OneDrive/Documentos/RStudio/estatR/R/hipotese.R', encoding='UTF-8')

# 1. Substituir fill='#D90429' por fill='red1' nas áreas de densidade (t e prop)
linhas <- gsub('fill = \"#D90429\", alpha = 0.35', 'fill = \"red1\", alpha = 0.35', linhas)

# 2. Remover grade dupla de todos, deixando apenas no boxplot
# Em teste_t_uma_amostra, teste_proporcao, teste_t_pareado:
# Vamos achar todos os meu_tema(grade = "dupla") e ver qual substituir.
# Em vez de regex complexa, vamos percorrer e substituir nas funcoes certas.

in_boxplot <- FALSE
for (i in 1:length(linhas)) {
  if (grepl('teste_t_duas_amostras <-', linhas[i])) in_boxplot <- TRUE
  if (grepl('teste_t_pareado <-', linhas[i])) in_boxplot <- FALSE
  
  if (grepl('meu_tema\\(grade = \"dupla\"\\)', linhas[i])) {
    if (!in_boxplot) {
      linhas[i] <- gsub('meu_tema\\(grade = \"dupla\"\\)', 'meu_tema()', linhas[i])
    }
  }
  
  # 3. Elevar a posição do t e Z no grafico
  if (grepl('y = max\\(df_curv\\\\) \\* 0\\.85', linhas[i])) {
    linhas[i] <- gsub('\\* 0\\.85', '* 1.05', linhas[i])
  }
}

# 4. Modificar o bloco do grafico em teste_t_duas_amostras (boxplot)
# Precisamos adicionar o bracket e remover o subtitle
idx_start_box_plot <- grep('p_plot <- ggplot2::ggplot\\(\\) \\+ eval\\(chamada_box\\) \\+', linhas)
if (length(idx_start_box_plot) > 0) {
  idx <- idx_start_box_plot[1]
  
  # Procurar onde ele define o p_sub antes do p_plot
  idx_p_sub <- grep('p_sub <- sprintf\\(\"p %s', linhas)[1]
  
  # Construir o novo codigo do grafico do boxplot
  novo_grafico <- c(
    '        max_y <- max(df_plot[[x_nome]], na.rm = TRUE)',
    '        min_y <- min(df_plot[[x_nome]], na.rm = TRUE)',
    '        amp <- max_y - min_y',
    '        y_bar <- max_y + amp * 0.06',
    '        y_tick <- max_y + amp * 0.03',
    '        y_text <- max_y + amp * 0.10',
    '',
    '        old_w <- getOption(\"warn\"); options(warn = -1)',
    '        p_plot <- ggplot2::ggplot() + eval(chamada_box) +',
    '          ggplot2::annotate(\"segment\", x = 1, xend = 2, y = y_bar, yend = y_bar, color = \"black\", linewidth = 0.5) +',
    '          ggplot2::annotate(\"segment\", x = 1, xend = 1, y = y_tick, yend = y_bar, color = \"black\", linewidth = 0.5) +',
    '          ggplot2::annotate(\"segment\", x = 2, xend = 2, y = y_tick, yend = y_bar, color = \"black\", linewidth = 0.5) +',
    '          ggplot2::annotate(\"text\", x = 1.5, y = y_text, label = p_sub, size = 4.5, fontface = \"bold\") +',
    '          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +',
    '          ggplot2::labs(',
    '            title    = sprintf(\"Compara\u00e7\u00e3o: %s por %s\", var_nome, grupo_nome_plot),',
    '            x        = grupo_nome_plot,',
    '            y        = var_nome,',
    '            caption  = \"estatR\"',
    '          ) +',
    '          meu_tema(grade = \"dupla\")'
  )
  
  # Achar onde termina o bloco labs atual
  idx_end_labs <- idx + grep('meu_tema\\(grade = \"dupla\"\\)', linhas[idx:length(linhas)])[1] - 1
  
  # Substituir
  linhas <- c(linhas[1:(idx-2)], novo_grafico, linhas[(idx_end_labs+1):length(linhas)])
}

# 5. Dar um expand no eixo Y para as curvas de densidade (t e prop) para o texto caber
# Procurar p_plot <- p_plot + em teste_t e teste_prop e injetar scale_y_continuous
idx_expand <- grep('ggplot2::geom_vline\\(xintercept = ', linhas)
for (i in rev(idx_expand)) {
  linhas <- c(linhas[1:(i-1)], 
              '          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +', 
              linhas[i:length(linhas)])
}

writeLines(linhas, 'C:/Users/Elize/OneDrive/Documentos/RStudio/estatR/R/hipotese.R', useBytes=TRUE)
