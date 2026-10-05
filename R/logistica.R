# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: REGRESSÃO LOGÍSTICA BINOMIAL
# ─────────────────────────────────────────────────────────────────────────────

# --- Funções Internas de Cálculo ----------------------------------------------

.calc_auc_estatR <- function(prob, y) {
  # AUC via U de Mann-Whitney (muito mais rápido e não depende de pROC)
  y_num <- as.numeric(y)
  if (!all(y_num %in% c(0, 1))) return(NA)
  
  R <- rank(prob)
  n1 <- sum(y_num == 1)
  n0 <- sum(y_num == 0)
  if (n1 == 0 || n0 == 0) return(NA)
  
  u <- sum(R[y_num == 1]) - n1 * (n1 + 1) / 2
  auc <- u / (n1 * n0)
  return(auc)
}

.hosmer_lemeshow_estatR <- function(y, prob, g = 10) {
  y_num <- as.numeric(y)
  q <- unique(quantile(prob, probs = seq(0, 1, by = 1/g)))
  if (length(q) < 3) return(list(statistic = NA, p.value = NA))
  
  q[1] <- -Inf
  q[length(q)] <- Inf
  
  cuty <- cut(prob, breaks = q)
  obs <- xtabs(y_num ~ cuty)
  n_total <- table(cuty)
  exp <- xtabs(prob ~ cuty)
  
  chi <- sum((obs - exp)^2 / (exp * (1 - exp/n_total)), na.rm = TRUE)
  pval <- 1 - pchisq(chi, g - 2)
  return(list(statistic = chi, p.value = pval))
}

.calc_roc_df <- function(prob, y) {
  y_num <- as.numeric(y)
  ord <- order(prob, decreasing = TRUE)
  prob <- prob[ord]
  y_num <- y_num[ord]
  
  tpr <- cumsum(y_num) / sum(y_num)
  fpr <- cumsum(1 - y_num) / sum(1 - y_num)
  
  data.frame(FPR = c(0, fpr), TPR = c(0, tpr))
}

.formata_p_local <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("< 0.001")
  return(sprintf("%.3f", p))
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Regressão Logística Binomial
#' @description Ajusta um modelo de regressão logística binomial, gerando um 
#' painel de diagnóstico completo com Deviance, Pseudo-R², AUC, Odds Ratio (OR),
#' matriz de confusão e teste de Hosmer-Lemeshow.
#' @param formula Fórmula do modelo (ex: y ~ x1 + x2).
#' @param dados Data frame contendo os dados.
#' @param corte_prob Ponto de corte (threshold) para classificação. Padrão 0.5.
#' @param grafico Se TRUE, exibe a Curva ROC do modelo.
#' @return Retorna invisivelmente o modelo ajustado (objeto glm).
#' @export
regressao_logistica <- function(formula, dados, corte_prob = 0.5, grafico = TRUE) {
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  
  y_name <- as.character(formula[[2]])
  if (!(y_name %in% names(dados))) stop("Variável resposta não encontrada.")
  
  # Preparando Y
  Y_raw <- na.omit(dados[[y_name]])
  niveis <- unique(Y_raw)
  
  if (length(niveis) > 2) {
    stop(sprintf("A variável resposta '%s' possui %d níveis. A regressão multinomial será adicionada em breve! Para a regressão logística binomial, a resposta deve ter exatamente 2 categorias.", y_name, length(niveis)))
  }
  if (length(niveis) < 2) {
    stop("A variável resposta deve possuir 2 categorias (variabilidade zero detectada).")
  }
  
  # Ajuste do modelo
  mod <- tryCatch(
    glm(formula, data = dados, family = binomial(link = "logit")),
    error = function(e) stop("Erro ao ajustar o modelo logístico: ", e$message)
  )
  
  # Dados reais extraídos do modelo (para garantir que NAs removidos batam com o mod)
  df_mod <- mod$model
  Y_fit <- df_mod[[1]]
  
  # Definindo Sucesso (1) e Fracasso (0)
  if (is.factor(Y_fit)) {
    sucesso_label <- levels(Y_fit)[2]
    Y_num <- ifelse(Y_fit == sucesso_label, 1, 0)
  } else {
    Y_num <- as.numeric(Y_fit)
    sucesso_label <- as.character(max(Y_num))
    Y_num <- ifelse(Y_num == max(Y_num), 1, 0)
  }
  
  probs <- mod$fitted.values
  
  # ── MÉTRICAS DE AJUSTE GLOBAL ──────────────────────────────────────────────
  dev_nula <- mod$null.deviance
  dev_res <- mod$deviance
  df_nula <- mod$df.null
  df_res <- mod$df.residual
  
  # Likelihood Ratio Test (LRT)
  lrt_chi <- dev_nula - dev_res
  lrt_df <- df_nula - df_res
  lrt_p <- 1 - pchisq(lrt_chi, lrt_df)
  lrt_sig <- if (lrt_p < 0.05) "Modelo globalmente significante" else "Modelo NÃO significante"
  
  aic_val <- AIC(mod)
  bic_val <- BIC(mod)
  
  # ── MÉTRICAS DE PERFORMANCE ────────────────────────────────────────────────
  pseudo_r2 <- 1 - (dev_res / dev_nula) # McFadden
  auc_val <- .calc_auc_estatR(probs, Y_num)
  
  # Matriz de Confusão
  pred_class <- ifelse(probs >= corte_prob, 1, 0)
  vp <- sum(pred_class == 1 & Y_num == 1)
  vn <- sum(pred_class == 0 & Y_num == 0)
  fp <- sum(pred_class == 1 & Y_num == 0)
  fn <- sum(pred_class == 0 & Y_num == 1)
  acc <- (vp + vn) / length(Y_num)
  
  # ── TESTE DE ADEQUAÇÃO ─────────────────────────────────────────────────────
  hl <- .hosmer_lemeshow_estatR(Y_num, probs, g = 10)
  hl_msg <- if (is.na(hl$p.value)) {
    "Não calculável (pouca variação de probabilidade)."
  } else if (hl$p.value > 0.05) {
    "Bom ajuste (As previsões calibram bem com o real)."
  } else {
    "Ajuste Pobre (As previsões desviam do real)."
  }
  
  # ── COEFICIENTES E ODDS RATIO ──────────────────────────────────────────────
  sm <- summary(mod)
  coefs <- sm$coefficients
  
  # Confint default é mais rápido e não trava com perfilagem
  ci <- suppressMessages(confint.default(mod))
  
  df_coef <- data.frame(
    "Variável" = rownames(coefs),
    "Estimativa (log-odds)" = sapply(coefs[,1], .fmt_num, decimais = 4),
    "Erro Padrão" = sapply(coefs[,2], .fmt_num, decimais = 4),
    "Z" = sapply(coefs[,3], .fmt_num, decimais = 2),
    "p-valor" = sapply(coefs[,4], .formata_p_local),
    "OR" = sapply(exp(coefs[,1]), .fmt_num, decimais = 3),
    "IC(OR) 95% Inf" = sapply(exp(ci[,1]), .fmt_num, decimais = 3),
    "IC(OR) 95% Sup" = sapply(exp(ci[,2]), .fmt_num, decimais = 3),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  rownames(df_coef) <- NULL
  
  # AICc
  K_params <- length(coef(mod)) + 1
  n_obs <- length(Y_num)
  aicc_val <- if (n_obs - K_params - 1 > 0) aic_val + (2 * K_params * (K_params + 1)) / (n_obs - K_params - 1) else Inf
  
  df_resumo <- data.frame(
    "Métrica" = c("Deviance Nula", "Deviance Residual", "Teste LRT (\u03c7\u00b2)", "AIC", "AICc", "BIC"),
    "Valor" = c(sprintf("%.2f (df = %d)", dev_nula, df_nula),
                sprintf("%.2f (df = %d)", dev_res, df_res),
                sprintf("%.2f (p = %s)", lrt_chi, .formata_p_local(lrt_p)),
                sprintf("%.2f", aic_val),
                sprintf("%.2f", aicc_val),
                sprintf("%.2f", bic_val)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  df_perf <- data.frame(
    "Métrica" = c("Pseudo-R\u00b2 (McFadden)", "Acurácia Global", "AUC (Curva ROC)"),
    "Valor" = c(.fmt_pct(pseudo_r2, 1),
                sprintf("%s (Corte: %.2f)", .fmt_pct(acc, 1), corte_prob),
                sprintf("%.3f", auc_val)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  if (!is.na(hl$statistic)) {
    df_hl <- data.frame(
      "Teste" = "Hosmer-Lemeshow (g=10)",
      "Estatística (\u03c7\u00b2)" = sprintf("%.2f", hl$statistic),
      "p-valor" = .formata_p_local(hl$p.value),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  
  # ── IMPRESSÃO DO PAINEL ────────────────────────────────────────────────────
  .print_titulo("REGRESSÃO LOGÍSTICA BINOMIAL")
  
  cat(sprintf("  Variável Resposta:   %s (Sucesso = '%s')\n", y_name, sucesso_label))
  cat(sprintf("  Observações Válidas: %d\n\n", length(Y_num)))
  
  .print_topico("RESUMO GERAL (AJUSTE)")
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_resumo, align = c("left", "right"))
  } else {
    print(df_resumo, row.names = FALSE)
  }
  
  .print_topico("MÉTRICAS DE PERFORMANCE")
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_perf, align = c("left", "right"))
  } else {
    print(df_perf, row.names = FALSE)
  }
  
  .print_topico("TESTE DE ADEQUAÇÃO (Pressuposto)")
  if (!is.na(hl$statistic)) {
    if (exists(".print_tabela_estatR", mode = "function")) {
      .print_tabela_estatR(df_hl, align = c("left", "center", "center"))
    } else {
      print(df_hl, row.names = FALSE)
    }
  } else {
    cat(sprintf("  Hosmer-Lemeshow: %s\n\n", hl_msg))
  }
  
  .print_topico("MATRIZ DE CONFUSÃO")
  cat(sprintf("                Predito Não (0)   Predito Sim (1)\n"))
  cat(sprintf("  Real Não (0)  %15d   %15d\n", vn, fp))
  cat(sprintf("  Real Sim (1)  %15d   %15d\n\n", fn, vp))
  
  .print_topico("COEFICIENTES E ODDS RATIO (OR)")
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_coef, align = c("left", "right", "right", "right", "center", "right", "right", "right"))
  } else {
    print(df_coef, row.names = FALSE)
  }
  
  .print_topico("INTERPRETAÇÕES E AVISOS")
  if (lrt_p < 0.05) {
    cat("  \u2022 Ajuste Global: O modelo é estatisticamente significante (LRT p < 0.05).\n")
  } else {
    cat("  \u2022 Ajuste Global: O modelo NÃO é estatisticamente significante (LRT p >= 0.05).\n")
  }
  if (!is.na(hl$p.value)) {
    if (hl$p.value > 0.05) {
      cat("  \u2022 Hosmer-Lemeshow: Bom ajuste. As probabilidades preditas calibram bem com a realidade.\n")
    } else {
      cat("  \u2022 Hosmer-Lemeshow: Ajuste pobre. As previsões desviam significativamente do real.\n")
    }
  }
  cat("\n")
  
  .print_rodape()
  
  # ── GRÁFICO ROC ────────────────────────────────────────────────────────────
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_roc <- .calc_roc_df(probs, Y_num)
    
    p_plot <- ggplot2::ggplot(df_roc, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 1.2) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::labs(
        title = "Curva ROC (Receiver Operating Characteristic)",
        subtitle = sprintf("Área Sob a Curva (AUC): %.3f", auc_val),
        x = "Taxa de Falsos Positivos",
        y = "Taxa de Verdadeiros Positivos"
      )
    
    if (exists("tema_estatR", mode = "function")) {
      p_plot <- p_plot + tema_estatR()
    } else {
      p_plot <- p_plot + ggplot2::theme_minimal()
    }
    
    suppressWarnings(print(p_plot))
  }
  
  invisible(mod)
}
