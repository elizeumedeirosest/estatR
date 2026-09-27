#' @title Regressão Linear Simples e Múltipla
#' @description Ajusta um modelo de regressão linear e retorna uma saída formatada com coeficientes, ajuste e interpretação.
#' @param formula Fórmula do modelo (ex: y ~ x)
#' @param dados Data frame com os dados
#' @param grafico Lógico. Se TRUE, gera um gráfico de dispersão para regressão simples.
#' @param ... Outros argumentos passados para o gráfico.
#' @export
regressao_linear <- function(formula, dados, grafico = TRUE, ...) {
  # 1. Ajustar o modelo
  modelo <- lm(formula, data = dados)
  
  # 2. Extrair dados
  resumo <- summary(modelo)
  coefs <- resumo$coefficients
  r2 <- resumo$r.squared
  r2_adj <- resumo$adj.r.squared
  rmse <- sqrt(mean(modelo$residuals^2))
  f_stat <- resumo$fstatistic
  
  # Calcular p-valor global
  if(!is.null(f_stat)) {
    p_global <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  } else {
    p_global <- NA
  }
  
  # 3. Formatar saída do console
  vars <- all.vars(formula)
  y_name <- vars[1]
  x_names <- vars[-1]
  
  # Helper interno para formatar p-valor
  formata_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("<0,001")
    return(formatC(p, format = "f", digits = 3, decimal.mark = ","))
  }
  
  # Helper interno para centralizar strings (textos e números)
  centraliza <- function(val, w) {
    s <- trimws(as.character(val))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  
  linha_divisoria <- paste0("\n", strrep("\u2500", 54), "\n\n")
  
  cat("\nREGRESSÃO LINEAR\n\n")
  cat(sprintf("Variável resposta: %s\n", y_name))
  cat(sprintf("Preditor(es): %s\n", paste(x_names, collapse = ", ")))
  cat(sprintf("Observações: %d\n", nrow(modelo$model)))
  
  cat(linha_divisoria)
  cat("COEFICIENTES\n\n")
  
  # Formatar tabela de coeficientes
  rn <- rownames(coefs)
  rn[rn == "(Intercept)"] <- "Intercepto"
  
  cat(sprintf("%-15s %s %s %s %s\n", "Termo", centraliza("Estimativa", 12), centraliza("EP", 8), centraliza("t", 8), centraliza("p", 10)))
  for (i in 1:nrow(coefs)) {
    cat(sprintf("%-15s %s %s %s %s\n",
                rn[i],
                centraliza(sprintf("%.2f", coefs[i, 1]), 12),
                centraliza(sprintf("%.2f", coefs[i, 2]), 8),
                centraliza(sprintf("%.2f", coefs[i, 3]), 8),
                centraliza(formata_p(coefs[i, 4]), 10)))
  }
  
  cat(linha_divisoria)
  cat("TABELA ANOVA\n\n")
  
  tabela_anova <- anova(modelo)
  cat(sprintf("%-15s %s %s %s %s %s\n", "FV", centraliza("GL", 6), centraliza("SQ", 10), centraliza("QM", 10), centraliza("F", 8), centraliza("p", 10)))
  rn_anova <- rownames(tabela_anova)
  rn_anova[rn_anova == "Residuals"] <- "Resíduos"
  
  for (i in 1:nrow(tabela_anova)) {
    if (is.na(tabela_anova[i, "F value"])) {
       cat(sprintf("%-15s %s %s %s %s %s\n",
                   rn_anova[i], 
                   centraliza(tabela_anova[i, "Df"], 6), 
                   centraliza(sprintf("%.2f", tabela_anova[i, "Sum Sq"]), 10),
                   centraliza(sprintf("%.2f", tabela_anova[i, "Mean Sq"]), 10), 
                   centraliza("", 8), 
                   centraliza("", 10)))
    } else {
       cat(sprintf("%-15s %s %s %s %s %s\n",
                   rn_anova[i], 
                   centraliza(tabela_anova[i, "Df"], 6), 
                   centraliza(sprintf("%.2f", tabela_anova[i, "Sum Sq"]), 10),
                   centraliza(sprintf("%.2f", tabela_anova[i, "Mean Sq"]), 10), 
                   centraliza(sprintf("%.2f", tabela_anova[i, "F value"]), 8), 
                   centraliza(formata_p(tabela_anova[i, "Pr(>F)"]), 10)))
    }
  }

  cat(linha_divisoria)
  cat("AJUSTE DO MODELO\n\n")
  cat(sprintf("%-20s %10s\n", "R\u00b2", formatC(r2, format = "f", digits = 3, decimal.mark = ",")))
  cat(sprintf("%-20s %10s\n", "R\u00b2 ajustado", formatC(r2_adj, format = "f", digits = 3, decimal.mark = ",")))
  cat(sprintf("%-20s %10.2f\n", "RMSE", rmse))
  if(!is.na(p_global)) {
    cat(sprintf("%-20s %10s\n", "F-estatística",
                paste0(formatC(f_stat[1], format="f", digits=2), " (p ",
                       ifelse(p_global < 0.001, "< 0,001",
                              paste0("= ", formatC(p_global, format="f", digits=3, decimal.mark=","))), ")")))
  }
  
  cat(linha_divisoria)
  cat("INTERPRETAÇÃO\n\n")
  
  cat(sprintf("O modelo de regressão linear foi ajustado para explicar\n%s a partir de %s.\n\n",
              y_name, paste(x_names, collapse = " e ")))
  
  cat(sprintf("O modelo explica aproximadamente %s da variabilidade\nobservada em %s.\n\n",
              paste0(formatC(r2 * 100, format = "f", digits = 1, decimal.mark = ","), "%"), y_name))
  
  # Interpretação dos preditores significativos
  for (i in 2:nrow(coefs)) {
    p_val <- coefs[i, 4]
    if (p_val < 0.05) {
      direcao <- ifelse(coefs[i, 1] > 0, "positiva", "negativa")
      mantendo <- ifelse(length(x_names) > 1, ", mantendo os demais constantes", "")
      p_texto <- ifelse(p_val < 0.001, "p < 0,001",
                        sprintf("p = %s", formatC(p_val, format="f", digits=3, decimal.mark=",")))
      cat(sprintf("O preditor '%s' apresentou associação %s estatisticamente\nsignificativa (%s) com %s%s.\n\n",
                  rn[i], direcao, p_texto, y_name, mantendo))
    }
  }
  
  # Interpretação prática: efeito por unidade de cada preditor
  for (i in 2:nrow(coefs)) {
    beta <- coefs[i, 1]
    nome_pred <- rn[i]
    
    if (length(x_names) > 1) {
      # Múltipla: listar os outros preditores mantidos constantes
      outros <- rn[2:nrow(coefs)]
      outros <- outros[outros != nome_pred]
      ceteris <- sprintf("Mantendo %s constante(s), ", paste(outros, collapse = " e "))
    } else {
      ceteris <- ""
    }
    
    if (beta > 0) {
      cat(sprintf("%sum aumento de uma unidade em %s está associado\na um aumento estimado de %.2f unidades em %s.\n\n",
                  ceteris, nome_pred, abs(beta), y_name))
    } else {
      cat(sprintf("%sum aumento de uma unidade em %s está associado\na uma redução estimada de %.2f unidades em %s.\n\n",
                  ceteris, nome_pred, abs(beta), y_name))
    }
  }
  
  # Gerar Gráfico de Dispersão se for regressão simples (1 preditor)
  if (grafico) {
    if (length(x_names) == 1) {
      tryCatch({
        if (exists("grafico_de_dispersao") && exists("meu_tema")) {
          # Extraindo os valores para a equação
          intercepto <- coefs[1, 1]
          inclinacao <- coefs[2, 1]
          sinal <- ifelse(inclinacao >= 0, "+", "-")
          equacao_texto <- sprintf("%s = %.2f %s %.2f * %s\nR\u00b2 = %.1f%%",
                                   y_name, intercepto, sinal, abs(inclinacao), x_names[1], r2 * 100)
          
          # Posição Inteligente: reta cai -> topo-direito | reta sobe -> topo-esquerdo
          if (inclinacao < 0) {
            pos_x <- Inf;  pos_y <- Inf; anc_h <- 1.1;  anc_v <- 1.5
          } else {
            pos_x <- -Inf; pos_y <- Inf; anc_h <- -0.1; anc_v <- 1.5
          }
          
          # Avalia strings como nomes de colunas (symbols) para a função do metaR
          chamada_grafico <- bquote(
            grafico_de_dispersao(data = dados, x = .(as.name(x_names[1])), y = .(as.name(y_name)),
                                 reta = TRUE, ic = TRUE, outlier = TRUE, deteccao_outliers = "residuos",
                                 paleta = 1)
          )
          
          # Silenciador de warnings apenas durante o print
          old_w <- getOption("warn")
          options(warn = -1)
          
          p_plot <- ggplot2::ggplot() + eval(chamada_grafico) +
            ggplot2::annotate("text", x = pos_x, y = pos_y, label = equacao_texto,
                              hjust = anc_h, vjust = anc_v, size = 6.5,
                              fontface = "bold", color = "#333333") +
            ggplot2::labs(
              title = sprintf("Regressão Linear: %s vs %s", y_name, x_names[1]),
              caption = "estatR"
            ) +
            meu_tema(grade = "dupla")
          
          suppressMessages(print(p_plot))
          options(warn = old_w)
        } else {
          message("\n[Aviso] Funções 'grafico_de_dispersao' ou 'meu_tema' não encontradas. Verifique se o metaR está carregado.")
        }
      }, error = function(e) {
        message("\n[Aviso] Não foi possível gerar o gráfico automático: ", e$message)
      })
    }
  }
  
  invisible(modelo)
}

#' @title Métricas do Modelo
#' @description Exibe as principais métricas de avaliação e erro para um modelo de regressão.
#' @param modelo Objeto do tipo lm
#' @export
metricas <- function(modelo) {
  resumo <- summary(modelo)
  r2     <- resumo$r.squared
  r2_adj <- resumo$adj.r.squared
  
  res <- modelo$residuals
  n   <- length(res)
  k   <- length(modelo$coefficients)
  
  rmse     <- sqrt(mean(res^2))
  mae      <- mean(abs(res))
  aic_val  <- AIC(modelo)
  bic_val  <- BIC(modelo)
  aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
  
  # Larguras fixas para garantir alinhamento perfeito
  lw <- 20   # largura da coluna do rótulo
  vw <- 10   # largura da coluna do valor
  linha <- paste0(strrep("\u2500", lw + vw + 1), "\n")
  
  f3 <- function(x) formatC(x, format = "f", digits = 3, decimal.mark = ",")
  f2 <- function(x) formatC(x, format = "f", digits = 2, decimal.mark = ",")
  
  cat("\nMÉTRICAS E CRITÉRIOS DE SELEÇÃO\n\n")
  cat(linha)
  cat(sprintf("%-*s %*s\n", lw, "Métrica", vw, "Valor"))
  cat(linha)
  cat(sprintf("%-*s %*s\n", lw, "R\u00b2",           vw, f3(r2)))
  cat(sprintf("%-*s %*s\n", lw, "R\u00b2 ajustado",  vw, f3(r2_adj)))
  cat(sprintf("%-*s %*s\n", lw, "RMSE",          vw, f3(rmse)))
  cat(sprintf("%-*s %*s\n", lw, "MAE",           vw, f3(mae)))
  cat(sprintf("%-*s %*s\n", lw, "AIC",           vw, f2(aic_val)))
  cat(sprintf("%-*s %*s\n", lw, "AICc",          vw, f2(aicc_val)))
  cat(sprintf("%-*s %*s\n", lw, "BIC",           vw, f2(bic_val)))
  cat(linha, "\n", sep = "")
  
  invisible(data.frame(
    Metrica = c("R2", "R2_ajustado", "RMSE", "MAE", "AIC", "AICc", "BIC"),
    Valor   = c(r2, r2_adj, rmse, mae, aic_val, aicc_val, bic_val)
  ))
}

#' @title Análise Residual
#' @description Realiza os testes estatísticos de diagnóstico dos resíduos e plota o painel gráfico internamente.
#' @param modelo Objeto do tipo lm
#' @param grafico Lógico. Se TRUE, exibe o painel 2x2 de diagnóstico dos resíduos.
#' @export
analise_residual <- function(modelo, grafico = TRUE) {
  res <- modelo$residuals
  
  # 1. Normalidade (Shapiro-Wilk)
  st     <- shapiro.test(res)
  p_norm <- st$p.value
  w_norm <- st$statistic
  
  # 2. Homocedasticidade (Breusch-Pagan, base R)
  res_sq       <- res^2
  dados_modelo <- modelo$model
  bp_mod       <- lm(res_sq ~ ., data = dados_modelo[, -1, drop = FALSE])
  r2_bp        <- summary(bp_mod)$r.squared
  n            <- length(res)
  bp_stat      <- n * r2_bp
  df_bp        <- length(modelo$coefficients) - 1
  p_bp         <- pchisq(bp_stat, df_bp, lower.tail = FALSE)
  
  formata_p <- function(p) {
    if (p < 0.001) return("<0,001")
    return(formatC(p, format = "f", digits = 3, decimal.mark = ","))
  }
  
  linha <- paste0(strrep("\u2500", 54), "\n")
  
  cat("\nDIAGNÓSTICO DOS RESÍDUOS\n\n")
  cat(linha)
  cat(sprintf("%-28s %-14s %s\n", "Teste", "Estatística", "p-valor"))
  cat(linha)
  cat(sprintf("%-28s W = %-10.3f %s\n", "Normalidade (Shapiro-Wilk)", w_norm, formata_p(p_norm)))
  cat(sprintf("%-28s BP = %-9.3f %s\n",  "Homocedasticidade (B-P)",   bp_stat, formata_p(p_bp)))
  cat(linha)
  
  cat("\nINTERPRETAÇÃO DO DIAGNÓSTICO\n\n")
  
  if (p_norm < 0.05) {
    cat(sprintf("[!] : O teste de normalidade rejeitou a hipótese nula\n    (p = %s). Os resíduos não seguem uma distribuição normal,\n    o que pode afetar a confiabilidade dos intervalos de confiança.\n\n",
                formata_p(p_norm)))
  } else {
    cat(sprintf("[\u2713] : O teste não encontrou evidências para rejeitar a normalidade\n    dos resíduos (p = %s).\n\n",
                formata_p(p_norm)))
  }
  
  if (p_bp < 0.05) {
    cat(sprintf("[!] : O teste detectou heterocedasticidade (p = %s),\n    indicando que a variância dos resíduos não é constante.\n\n",
                formata_p(p_bp)))
  } else {
    cat(sprintf("[\u2713] : O teste indicou homocedasticidade (p = %s),\n    ou seja, a variância dos resíduos é considerada constante\n    ao longo dos valores ajustados.\n\n",
                formata_p(p_bp)))
  }
  
  if (grafico) {
    tryCatch({
      df_res <- data.frame(
        ajustados = fitted(modelo),
        residuos  = res,
        res_pad   = rstandard(modelo),
        leverage  = hatvalues(modelo)
      )
      
      tema_painel <- if (exists("meu_tema")) meu_tema(grade = "dupla") else ggplot2::theme_minimal(base_size = 11)
      
      p1 <- ggplot2::ggplot(df_res, ggplot2::aes(x = ajustados, y = residuos)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
        ggplot2::geom_point(color = "#555555", alpha = 0.8) +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429",
                             linewidth = 0.7, formula = y ~ x) +
        ggplot2::labs(title = "Resíduos vs Ajustados",
                      x = "Valores Ajustados", y = "Resíduos") +
        tema_painel
      
      p2 <- ggplot2::ggplot(df_res, ggplot2::aes(sample = res_pad)) +
        ggplot2::stat_qq(color = "#555555", alpha = 0.8) +
        ggplot2::stat_qq_line(color = "#D90429", linewidth = 0.7) +
        ggplot2::labs(title = "QQ-Plot Normal",
                      x = "Quantis Teóricos", y = "Quantis Amostrais") +
        tema_painel
      
      p3 <- ggplot2::ggplot(df_res, ggplot2::aes(x = ajustados, y = sqrt(abs(res_pad)))) +
        ggplot2::geom_point(color = "#555555", alpha = 0.8) +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429",
                             linewidth = 0.7, formula = y ~ x) +
        ggplot2::labs(title = "Escala-Localização",
                      x = "Valores Ajustados", y = "\u221a|Res. Padronizados|") +
        tema_painel
      
      p4 <- ggplot2::ggplot(df_res, ggplot2::aes(x = res_pad)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 10, fill = "#555555", color = "black", alpha = 0.8) +
        ggplot2::geom_density(color = "#D90429", linewidth = 0.8) +
        ggplot2::labs(title = "Histograma dos Resíduos",
                      x = "Resíduos Padronizados", y = "Densidade",
                      caption = "estatR") +
        tema_painel
      
      if (requireNamespace("patchwork", quietly = TRUE)) {
        painel <- (p1 | p2) / (p3 | p4)
        suppressMessages(suppressWarnings(print(painel)))
      } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        suppressMessages(suppressWarnings(
          gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
        ))
      } else {
        suppressMessages(suppressWarnings({
          print(p1); print(p2); print(p3); print(p4)
        }))
      }
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o painel de resíduos: ", e$message)
    })
  }
  
  invisible(list(shapiro = st, bp = list(statistic = bp_stat, p.value = p_bp)))
}
