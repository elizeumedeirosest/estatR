linhas <- readLines('C:/Users/Elize/OneDrive/Documentos/RStudio/estatR/R/hipotese.R', encoding='UTF-8')

# 1. Substituir cores e remover grade onde nao eh boxplot
linhas <- gsub('fill = "#D90429", alpha = 0.35', 'fill = "red1", alpha = 0.35', linhas)

in_boxplot <- FALSE
for (i in 1:length(linhas)) {
  if (grepl('teste_t_duas_amostras <-', linhas[i])) in_boxplot <- TRUE
  if (grepl('teste_t_pareado <-', linhas[i])) in_boxplot <- FALSE
  
  if (grepl('meu_tema\\(grade = "dupla"\\)', linhas[i]) && !in_boxplot) {
    linhas[i] <- gsub('meu_tema\\(grade = "dupla"\\)', 'meu_tema()', linhas[i])
  }
  
  if (grepl('y = max\\(df_curv\\$y\\) \\* 0\\.85', linhas[i])) {
    linhas[i] <- gsub('\\* 0\\.85', '* 1.05', linhas[i])
  }
}

# 2. Modificar o boxplot (teste_t_duas_amostras)
idx_box <- grep('p_plot <- ggplot2::ggplot\\(\\) \\+ eval\\(chamada_box\\)', linhas)
if (length(idx_box) > 0) {
  idx <- idx_box[1]
  
  novo_grafico <- c(
    '        max_y <- max(df_plot[[x_nome]], na.rm = TRUE)',
    '        min_y <- min(df_plot[[x_nome]], na.rm = TRUE)',
    '        amp <- max_y - min_y',
    '        y_bar <- max_y + amp * 0.08',
    '        y_tick <- max_y + amp * 0.05',
    '        y_text <- max_y + amp * 0.12',
    '',
    '        old_w <- getOption("warn"); options(warn = -1)',
    '        p_plot <- ggplot2::ggplot() + eval(chamada_box) +',
    '          ggplot2::annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, color = "black", linewidth = 0.6) +',
    '          ggplot2::annotate("segment", x = 1, xend = 1, y = y_tick, yend = y_bar, color = "black", linewidth = 0.6) +',
    '          ggplot2::annotate("segment", x = 2, xend = 2, y = y_tick, yend = y_bar, color = "black", linewidth = 0.6) +',
    '          ggplot2::annotate("text", x = 1.5, y = y_text, label = p_sub, size = 4.5, fontface = "bold") +',
    '          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.20))) +',
    '          ggplot2::labs(',
    '            title    = sprintf("Compara\u00e7\u00e3o: %s por %s", var_nome, grupo_nome_plot),',
    '            x        = grupo_nome_plot,',
    '            y        = var_nome,',
    '            caption  = "estatR"',
    '          ) +',
    '          meu_tema(grade = "dupla")'
  )
  
  idx_end <- idx + grep('meu_tema\\(grade = "dupla"\\)', linhas[idx:length(linhas)])[1] - 1
  linhas <- c(linhas[1:(idx-2)], novo_grafico, linhas[(idx_end+1):length(linhas)])
}

# 3. Dar espaço no topo dos graficos de densidade para nao cortar o T / Z
idx_expand <- grep('ggplot2::geom_vline\\(xintercept =', linhas)
for (i in rev(idx_expand)) {
  linhas <- c(linhas[1:(i-1)], 
              '          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.20))) +', 
              linhas[i:length(linhas)])
}

writeLines(linhas, 'C:/Users/Elize/OneDrive/Documentos/RStudio/estatR/R/hipotese.R', useBytes=TRUE)
