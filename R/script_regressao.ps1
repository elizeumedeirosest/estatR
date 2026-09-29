# Le o arquivo original
linhas <- readLines('C:/Users/Elize/OneDrive/Documentos/RStudio/estatR/R/regressao.R', encoding = 'UTF-8')

# Identifica o inicio de regressao_linear e metricas
ini_reg <- grep('^regressao_linear <- function', linhas) - 7 # Para pegar os comentarios
ini_met <- grep('^#\\' @title Métricas do Modelo', linhas)

# Pega apenas de metricas pra baixo
linhas_restantes <- linhas[ini_met:length(linhas)]

# Cria as novas funcoes
novas_funcoes <- "
#' @title Regressão Linear Simples e Múltipla
#' @description Ajusta um modelo de regressão linear e retorna uma saída formatada com coeficientes, ajuste e interpretação.
#' @param formula Fórmula do modelo (ex: y ~ x)
#' @param dados Data frame com os dados
#' @param grafico Lógico. Se TRUE, gera um gráfico de dispersão para regressão simples.
#' @param ... Outros argumentos passados para o gráfico.
#' @export
regressao_linear <- function(formula, dados, grafico = TRUE, ...) {
  modelo <- lm(formula, data = dados)
  resumo <- summary(modelo)
  coefs <- resumo\\
  r2 <- resumo\\.squared
  r2_adj <- resumo\\.r.squared
  rmse <- sqrt(mean(modelo\\^2))
  f_stat <- resumo\\
  
  if (!is.null(f_stat)) {
    p_global <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  } else {
    p_global <- NA
  }
  
  tem_vif <- FALSE
  vif_vals <- NULL
  if (length(modelo\\) > 2) {
    vif_vals <- rep(NA, length(modelo\\) - 1)
    nomes_x <- all.vars(formula)[-1]
    try({
      X <- model.matrix(modelo)[, -1, drop=FALSE]
      for (i in 1:ncol(X)) {
        if (ncol(X) > 1) {
          mod_vif <- lm(X[, i] ~ X[, -i])
          r2_vif <- summary(mod_vif)\\.squared
          vif_vals[i] <- 1 / (1 - r2_vif)
        }
      }
      tem_vif <- TRUE
    }, silent = TRUE)
  }
  
  ic_vals <- confint(modelo, level = 0.95)
  
  formata_p <- function(p) {
    if (is.na(p)) return(\"NA\")
    if (p < 0.001) return(\"< 0.001\")
    return(formatC(p, format = \"f\", digits = 3, decimal.mark = \",\"))
  }
  
  get_asterisk <- function(p) {
    if (is.na(p)) return(\"\")
    if (p < 0.01) return(\"***\")
    if (p < 0.05) return(\"**\")
    if (p < 0.10) return(\"*\")
    return(\"\")
  }
  
  pad <- function(s, w, align = \"left\") {
    s <- as.character(s)
    n <- nchar(s)
    if (n >= w) return(s)
    sp <- w - n
    if (align == \"left\") paste0(s, strrep(\" \", sp))
    else if (align == \"right\") paste0(strrep(\" \", sp), s)
    else paste0(strrep(\" \", floor(sp/2)), s, strrep(\" \", ceiling(sp/2)))
  }
  
  vars <- all.vars(formula)
  y_name <- vars[1]
  x_names <- vars[-1]
  
  cat(sprintf(\"\\n── RESULTADOS DA REGRESSÃO %s ──\\n\", ifelse(length(x_names) > 1, \"MÚLTIPLA\", \"SIMPLES\")))
  cat(sprintf(\"  Equação do modelo: %s\\n\", deparse(formula)))
  cat(sprintf(\"  Observações: %d\\n\\n\", nrow(modelo\\)))
  
  cat(\"  ▶ QUADRO ANOVA\\n\")
  tabela_anova <- anova(modelo)
  
  w_a <- c(fv=20, gl=12, sq=12, qm=12, f=9, p=12)
  hdr_anova <- paste(pad(\"Fonte de Variação\", w_a[\"fv\"], \"left\"),
                     pad(\"Graus Lib.\", w_a[\"gl\"], \"center\"),
                     pad(\"Soma Quad.\", w_a[\"sq\"], \"center\"),
                     pad(\"Quad. Médio\", w_a[\"qm\"], \"center\"),
                     pad(\"F\", w_a[\"f\"], \"center\"),
                     pad(\"p-valor\", w_a[\"p\"], \"center\"))
  sep_anova <- paste(rep(\"─\", nchar(hdr_anova) + 2), collapse=\"\")
  
  cat(\"  \", sep_anova, \"\\n\", sep=\"\")
  cat(\"  \", hdr_anova, \"\\n\", sep=\"\")
  cat(\"  \", sep_anova, \"\\n\", sep=\"\")
  
  rn_anova <- rownames(tabela_anova)
  idx_res <- which(rownames(tabela_anova) == \"Residuals\")
  df_res <- tabela_anova[idx_res, \"Df\"]
  sq_res <- tabela_anova[idx_res, \"Sum Sq\"]
  qm_res <- tabela_anova[idx_res, \"Mean Sq\"]
  
  df_mod <- sum(tabela_anova[-idx_res, \"Df\"])
  sq_mod <- sum(tabela_anova[-idx_res, \"Sum Sq\"])
  qm_mod <- sq_mod / df_mod
  f_val_mod <- f_stat[1]
  p_val_mod <- p_global
  
  df_tot <- df_mod + df_res
  sq_tot <- sq_mod + sq_res
  
  cat(\"  \", paste(
    pad(\"Modelo (Regressão)\", w_a[\"fv\"], \"left\"),
    pad(df_mod, w_a[\"gl\"], \"center\"),
    pad(sprintf(\"%.2f\", sq_mod), w_a[\"sq\"], \"center\"),
    pad(sprintf(\"%.2f\", qm_mod), w_a[\"qm\"], \"center\"),
    pad(if(!is.na(f_val_mod)) sprintf(\"%.2f\", f_val_mod) else \"\", w_a[\"f\"], \"center\"),
    pad(paste(formata_p(p_val_mod), get_asterisk(p_val_mod)), w_a[\"p\"], \"center\")
  ), \"\\n\", sep=\"\")
  
  cat(\"  \", paste(
    pad(\"Resíduos (Erro)\", w_a[\"fv\"], \"left\"),
    pad(df_res, w_a[\"gl\"], \"center\"),
    pad(sprintf(\"%.2f\", sq_res), w_a[\"sq\"], \"center\"),
    pad(sprintf(\"%.2f\", qm_res), w_a[\"qm\"], \"center\"),
    pad(\"\", w_a[\"f\"], \"center\"),
    pad(\"\", w_a[\"p\"], \"center\")
  ), \"\\n\", sep=\"\")
  
  cat(\"  \", paste(
    pad(\"Total\", w_a[\"fv\"], \"left\"),
    pad(df_tot, w_a[\"gl\"], \"center\"),
    pad(sprintf(\"%.2f\", sq_tot), w_a[\"sq\"], \"center\"),
    pad(\"\", w_a[\"qm\"], \"center\"),
    pad(\"\", w_a[\"f\"], \"center\"),
    pad(\"\", w_a[\"p\"], \"center\")
  ), \"\\n\", sep=\"\")
  
  cat(\"  \", sep_anova, \"\\n\", sep=\"\")
  cat(\"  Notas:\\n\")
  cat(\"  Significância: *** p < 0.01   ** p < 0.05   * p < 0.10\\n\\n\")
  
  cat(\"  ▶ COEFICIENTES DO MODELO\\n\")
  
  w_c <- c(var=18, est=12, err=12, ic=16, p=14, vif=8)
  if(!tem_vif) w_c[\"vif\"] <- 0
  
  hdr_coef <- paste(pad(\"Variável\", w_c[\"var\"], \"left\"),
                    pad(\"Estimativa\", w_c[\"est\"], \"center\"),
                    pad(\"Erro Pad.\", w_c[\"err\"], \"center\"),
                    pad(\"IC (95%)\", w_c[\"ic\"], \"center\"),
                    pad(\"p-valor\", w_c[\"p\"], \"center\"),
                    if(tem_vif) pad(\"VIF\", w_c[\"vif\"], \"center\") else \"\")
  sep_coef <- paste(rep(\"─\", nchar(hdr_coef) + 2), collapse=\"\")
  
  cat(\"  \", sep_coef, \"\\n\", sep=\"\")
  cat(\"  \", hdr_coef, \"\\n\", sep=\"\")
  cat(\"  \", sep_coef, \"\\n\", sep=\"\")
  
  rn <- rownames(coefs)
  rn[rn == \"(Intercept)\"] <- \"(Intercepto)\"
  
  for (i in 1:nrow(coefs)) {
    ic_str <- sprintf(\"[%s; %s]\", formatC(ic_vals[i, 1], format=\"f\", digits=2), formatC(ic_vals[i, 2], format=\"f\", digits=2))
    p_ast <- paste(formata_p(coefs[i, 4]), get_asterisk(coefs[i, 4]))
    
    vif_str <- \"\"
    if(tem_vif) {
      if(i == 1) vif_str <- \"-\"
      else {
        v <- vif_vals[i - 1]
        vif_str <- if(is.na(v)) \"-\" else sprintf(\"%.2f\", v)
      }
    }
    
    linha_c <- paste(
      pad(rn[i], w_c[\"var\"], \"left\"),
      pad(sprintf(\"%.3f\", coefs[i, 1]), w_c[\"est\"], \"center\"),
      pad(sprintf(\"%.3f\", coefs[i, 2]), w_c[\"err\"], \"center\"),
      pad(ic_str, w_c[\"ic\"], \"center\"),
      pad(p_ast, w_c[\"p\"], \"center\"),
      if(tem_vif) pad(vif_str, w_c[\"vif\"], \"center\") else \"\"
    )
    cat(\"  \", trimws(linha_c, \"right\"), \"\\n\", sep=\"\")
  }
  cat(\"  \", sep_coef, \"\\n\", sep=\"\")
  cat(\"  Notas:\\n\")
  cat(\"  Significância: *** p < 0.01   ** p < 0.05   * p < 0.10\\n\")
  if(tem_vif) {
    cat(\"  VIF: Valores > 5 indicam multicolinearidade moderada; > 10 grave.\\n\")
  }
  cat(\"\\n\")
  
  cat(\"  ▶ AJUSTE DO MODELO\\n\")
  w_aj <- c(med = 22, val = 30)
  hdr_aj <- paste(pad(\"Medidas de Ajuste\", w_aj[\"med\"], \"left\"), pad(\"Valor\", w_aj[\"val\"], \"left\"))
  sep_aj <- paste(rep(\"─\", nchar(hdr_aj) + 2), collapse=\"\")
  
  cat(\"  \", sep_aj, \"\\n\", sep=\"\")
  cat(\"  \", hdr_aj, \"\\n\", sep=\"\")
  cat(\"  \", sep_aj, \"\\n\", sep=\"\")
  
  cat(\"  \", paste(pad(\"R² (Explicação):\", w_aj[\"med\"], \"left\"), pad(paste0(formatC(r2 * 100, format=\"f\", digits=2), \"%\"), w_aj[\"val\"], \"left\")), \"\\n\", sep=\"\")
  cat(\"  \", paste(pad(\"R² ajustado:\", w_aj[\"med\"], \"left\"), pad(paste0(formatC(r2_adj * 100, format=\"f\", digits=2), \"%\"), w_aj[\"val\"], \"left\")), \"\\n\", sep=\"\")
  cat(\"  \", paste(pad(\"RMSE:\", w_aj[\"med\"], \"left\"), pad(sprintf(\"%.2f\", rmse), w_aj[\"val\"], \"left\")), \"\\n\", sep=\"\")
  
  f_txt <- if(!is.na(p_global)) {
    sprintf(\"%.2f (p %s %s)\", f_stat[1], ifelse(p_global < 0.001, \"< 0.001\", paste0(\"= \", formatC(p_global, format=\"f\", digits=3, decimal.mark=\",\"))), get_asterisk(p_global))
  } else { \"NA\" }
  cat(\"  \", paste(pad(\"F-estatística:\", w_aj[\"med\"], \"left\"), pad(f_txt, w_aj[\"val\"], \"left\")), \"\\n\", sep=\"\")
  cat(\"  \", sep_aj, \"\\n\\n\", sep=\"\")
  
  cat(\"  ▶ INTERPRETAÇÃO\\n\")
  cat(sprintf(\"  O modelo de regressão linear foi ajustado para explicar %s a partir de %s.\\n\",
              y_name, paste(x_names, collapse = \" e \")))
  
  valido_txt <- if(!is.na(p_global) && p_global < 0.05) \"é globalmente válido\" else \"não é globalmente válido\"
  
  cat(sprintf(\"\\n  O modelo %s (p %s) e explica aproximadamente %s\\n  da variabilidade observada em %s.\\n\\n\",
              valido_txt,
              ifelse(is.na(p_global), \"NA\", ifelse(p_global < 0.001, \"< 0.001\", paste0(\"= \", formatC(p_global, format=\"f\", digits=3, decimal.mark=\",\")))),
              paste0(formatC(r2 * 100, format = \"f\", digits = 1, decimal.mark = \",\"), \"%\"), y_name))
              
  for (i in 2:nrow(coefs)) {
    p_val <- coefs[i, 4]
    if (p_val < 0.05) {
      direcao <- ifelse(coefs[i, 1] > 0, \"positiva\", \"negativa\")
      mantendo <- ifelse(length(x_names) > 1, \", mantendo os demais constantes\", \"\")
      p_texto <- ifelse(p_val < 0.001, \"p < 0.001\",
                        sprintf(\"p = %s\", formatC(p_val, format=\"f\", digits=3, decimal.mark=\",\")))
      cat(sprintf(\"  O preditor '%s' apresentou associação %s estatisticamente significativa\\n  (%s) com %s%s.\\n\\n\",
                  rn[i], direcao, p_texto, y_name, mantendo))
    }
  }
  
  for (i in 2:nrow(coefs)) {
    beta <- coefs[i, 1]
    nome_pred <- rn[i]
    
    if (length(x_names) > 1) {
      outros <- rn[2:nrow(coefs)]
      outros <- outros[outros != nome_pred]
      ceteris <- sprintf(\"Mantendo %s constante(s), um\", paste(outros, collapse = \" e \"))
    } else {
      ceteris <- \"Um\"
    }
    
    if (beta > 0) {
      cat(sprintf(\"  %s aumento de uma unidade em %s está associado\\n  a um aumento estimado de %.2f unidades em %s.\\n\\n\",
                  ceteris, nome_pred, abs(beta), y_name))
    } else {
      cat(sprintf(\"  %s aumento de uma unidade em %s está associado\\n  a uma redução estimada de %.2f unidades em %s.\\n\\n\",
                  ceteris, nome_pred, abs(beta), y_name))
    }
  }
  
  cat(strrep(\"─\", nchar(sep_anova)), \"\\n\")
  
  if (grafico) {
    if (length(x_names) == 1) {
      tryCatch({
        if (exists(\"grafico_de_dispersao\") && exists(\"meu_tema\")) {
          intercepto <- coefs[1, 1]
          inclinacao <- coefs[2, 1]
          sinal <- ifelse(inclinacao >= 0, \"+\", \"-\")
          equacao_texto <- sprintf(\"%s = %.2f %s %.2f * %s\\nR² = %.1f%%\",
                                   y_name, intercepto, sinal, abs(inclinacao), x_names[1], r2 * 100)
          
          if (inclinacao < 0) {
            pos_x <- Inf;  pos_y <- Inf; anc_h <- 1.1;  anc_v <- 1.5
          } else {
            pos_x <- -Inf; pos_y <- Inf; anc_h <- -0.1; anc_v <- 1.5
          }
          
          chamada_grafico <- bquote(
            grafico_de_dispersao(data = dados, x = .(as.name(x_names[1])), y = .(as.name(y_name)),
                                 reta = TRUE, ic = TRUE, outlier = TRUE, deteccao_outliers = \"residuos\",
                                 paleta = 1)
          )
          
          old_w <- getOption(\"warn\")
          options(warn = -1)
          
          p_plot <- ggplot2::ggplot() + eval(chamada_grafico) +
            ggplot2::annotate(\"text\", x = pos_x, y = pos_y, label = equacao_texto,
                              hjust = anc_h, vjust = anc_v, size = 6.5,
                              fontface = \"bold\", color = \"#333333\") +
            ggplot2::labs(
              title = sprintf(\"Regressão Linear: %s vs %s\", y_name, x_names[1]),
              caption = \"estatR\"
            ) +
            meu_tema(grade = \"dupla\")
          
          suppressMessages(print(p_plot))
          options(warn = old_w)
        }
      }, error = function(e) {})
    }
  }
  
  invisible(modelo)
}

"

# Funcao de predicao
funcao_pred <- "
#' @title Predição com o Modelo de Regressão
#' @description Gera previsões usando um modelo linear ajustado, com intervalos de confiança ou predição.
#' @param modelo Objeto do tipo lm (retornado pela função regressao_linear).
#' @param novos_dados Data frame contendo os novos dados para predição.
#' @param intervalo Tipo de intervalo: \"confianca\", \"predicao\" ou \"nenhum\". Padrão é \"confianca\".
#' @param nivel Nível de confiança (padrão 0.95).
#' @param decimais Casas decimais para exibição (padrão 2).
#' @param grafico Lógico. Se TRUE e for regressão simples, plota o resultado.
#' @return Retorna invisivelmente um data.frame com as previsões.
#' @export
predicao <- function(modelo, novos_dados, intervalo = c(\"confianca\", \"predicao\", \"nenhum\"), nivel = 0.95, decimais = 2, grafico = TRUE) {
  intervalo <- match.arg(intervalo)
  
  if (!inherits(modelo, \"lm\")) stop(\"O objeto 'modelo' deve ser da classe 'lm'.\")
  if (!is.data.frame(novos_dados)) stop(\"'novos_dados' deve ser um data.frame.\")
  
  if (intervalo == \"nenhum\") {
    preds <- predict(modelo, newdata = novos_dados)
    df_res <- data.frame(Previsao = preds)
  } else {
    tipo_int <- ifelse(intervalo == \"confianca\", \"confidence\", \"prediction\")
    preds <- predict(modelo, newdata = novos_dados, interval = tipo_int, level = nivel)
    df_res <- as.data.frame(preds)
    colnames(df_res) <- c(\"Previsao\", \"LI\", \"LS\")
  }
  
  y_name <- as.character(formula(modelo)[[2]])
  x_names <- all.vars(formula(modelo))[-1]
  
  df_x <- novos_dados[, x_names, drop = FALSE]
  df_final <- cbind(df_x, df_res)
  
  cat(sprintf(\"\\n── PREVISÕES DO MODELO (Variável: %s) ──\\n\", y_name))
  
  if (intervalo != \"nenhum\") {
    nome_int <- ifelse(intervalo == \"confianca\", \"Confiança\", \"Predição\")
    cat(sprintf(\"  Tipo de Margem: Intervalo de %s (%.0f%%)\\n\\n\", nome_int, nivel * 100))
  } else {
    cat(\"  Margens de erro (intervalos) não solicitadas.\\n\\n\")
  }
  
  pad <- function(s, w, align = \"left\") {
    s <- as.character(s)
    n <- nchar(s)
    if (n >= w) return(s)
    sp <- w - n
    if (align == \"left\") paste0(s, strrep(\" \", sp))
    else if (align == \"right\") paste0(strrep(\" \", sp), s)
    else paste0(strrep(\" \", floor(sp/2)), s, strrep(\" \", ceiling(sp/2)))
  }
  
  w_x <- pmax(sapply(x_names, nchar), sapply(df_x, function(col) max(nchar(as.character(col))))) + 2
  w_x <- pmax(w_x, 8)
  w_p <- 14
  w_l <- 20
  
  hdr_x <- paste(mapply(pad, x_names, w_x, \"left\"), collapse = \"\")
  
  if (intervalo == \"nenhum\") {
    hdr <- paste0(hdr_x, pad(\"Previsão (Y)\", w_p, \"center\"))
  } else {
    hdr <- paste0(hdr_x, pad(\"Previsão (Y)\", w_p, \"center\"), pad(\"Lim Inferior (LI)\", w_l, \"center\"), pad(\"Lim Superior (LS)\", w_l, \"center\"))
  }
  
  sep_hdr <- paste(rep(\"─\", nchar(hdr) + 2), collapse=\"\")
  
  cat(\"  \", sep_hdr, \"\\n\", sep=\"\")
  cat(\"  \", hdr, \"\\n\", sep=\"\")
  cat(\"  \", sep_hdr, \"\\n\", sep=\"\")
  
  for (i in 1:nrow(df_final)) {
    linha_x <- paste(mapply(function(nm, w) pad(df_final[i, nm], w, \"left\"), x_names, w_x), collapse = \"\")
    
    if (intervalo == \"nenhum\") {
      linha_y <- pad(sprintf(paste0(\"%.\", decimais, \"f\"), df_final[i, \"Previsao\"]), w_p, \"center\")
      cat(\"  \", linha_x, linha_y, \"\\n\", sep=\"\")
    } else {
      linha_y <- paste0(
        pad(sprintf(paste0(\"%.\", decimais, \"f\"), df_final[i, \"Previsao\"]), w_p, \"center\"),
        pad(sprintf(paste0(\"%.\", decimais, \"f\"), df_final[i, \"LI\"]), w_l, \"center\"),
        pad(sprintf(paste0(\"%.\", decimais, \"f\"), df_final[i, \"LS\"]), w_l, \"center\")
      )
      cat(\"  \", linha_x, linha_y, \"\\n\", sep=\"\")
    }
  }
  cat(\"  \", sep_hdr, \"\\n\\n\", sep=\"\")
  
  if (grafico && length(x_names) == 1 && intervalo != \"nenhum\") {
    tryCatch({
      if (exists(\"grafico_de_dispersao\") && exists(\"meu_tema\")) {
        dados_orig <- modelo\\
        
        chamada_grafico <- bquote(
          grafico_de_dispersao(data = dados_orig, x = .(as.name(x_names[1])), y = .(as.name(y_name)),
                               reta = TRUE, ic = TRUE, nivel_ic = .(nivel), paleta = 1)
        )
        
        old_w <- getOption(\"warn\")
        options(warn = -1)
        
        p_plot <- ggplot2::ggplot() + eval(chamada_grafico) +
          ggplot2::geom_point(data = df_final, ggplot2::aes(x = .data[[x_names[1]]], y = Previsao),
                              color = \"#D90429\", size = 4) +
          ggplot2::geom_errorbar(data = df_final, ggplot2::aes(x = .data[[x_names[1]]], ymin = LI, ymax = LS),
                                 color = \"#D90429\", width = 0.2, linewidth = 1) +
          ggplot2::labs(
            title = sprintf(\"Predição: %s vs %s\", y_name, x_names[1]),
            subtitle = sprintf(\"Intervalo de %s (%.0f%%)\", nome_int, nivel * 100),
            caption = \"estatR\"
          ) +
          meu_tema(grade = \"dupla\")
        
        suppressMessages(print(p_plot))
        options(warn = old_w)
      }
    }, error = function(e) {})
  }
  
  invisible(df_final)
}
"

arquivo_final <- paste(c(novas_funcoes, linhas_restantes, funcao_pred), collapse = "\n")
writeLines(arquivo_final, 'C:/Users/Elize/OneDrive/Documentos/RStudio/estatR/R/regressao.R', useBytes = TRUE)
