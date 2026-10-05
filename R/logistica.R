# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: REGRESSÃO LOGÍSTICA BINOMIAL
# ─────────────────────────────────────────────────────────────────────────────

# --- Funções Internas de Cálculo ----------------------------------------------

.calc_auc_estatR <- function(prob, y) {
  y_num <- as.numeric(y)
  if (!all(y_num %in% c(0, 1))) return(NA)
  R   <- rank(prob)
  n1  <- sum(y_num == 1)
  n0  <- sum(y_num == 0)
  if (n1 == 0 || n0 == 0) return(NA)
  u   <- sum(R[y_num == 1]) - n1 * (n1 + 1) / 2
  u / (n1 * n0)
}

.hosmer_lemeshow_estatR <- function(y, prob, g = 10) {
  y_num <- as.numeric(y)
  q     <- unique(quantile(prob, probs = seq(0, 1, by = 1/g)))
  if (length(q) < 3) return(list(statistic = NA, p.value = NA))
  q[1] <- -Inf; q[length(q)] <- Inf
  cuty  <- cut(prob, breaks = q)
  obs   <- as.numeric(xtabs(y_num ~ cuty))
  n_grp <- as.numeric(table(cuty))
  exp   <- as.numeric(xtabs(prob ~ cuty))
  denom <- exp * (1 - exp / n_grp)
  denom[denom <= 0] <- NA
  chi   <- sum((obs - exp)^2 / denom, na.rm = TRUE)
  list(statistic = chi, p.value = 1 - pchisq(chi, g - 2))
}

.calc_roc_df <- function(prob, y) {
  y_num <- as.numeric(y)
  ord   <- order(prob, decreasing = TRUE)
  y_ord <- y_num[ord]
  tpr   <- cumsum(y_ord)     / sum(y_ord)
  fpr   <- cumsum(1 - y_ord) / sum(1 - y_ord)
  data.frame(FPR = c(0, fpr), TPR = c(0, tpr))
}

.formata_p_local <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

.fmt_or <- function(v) {
  if (is.na(v) || is.infinite(v)) return("-")
  if (abs(v) >= 1e4 || (abs(v) < 0.001 && v != 0)) return(formatC(v, format = "e", digits = 2))
  sprintf("%.3f", v)
}

.asterisk_log <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.10) return("*")
  ""
}

.pad_log <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
  paste0(s, strrep(" ", sp))
}

.sep_log <- function(w) paste0("  ", strrep("\u2500", w))

# ─────────────────────────────────────────────────────────────────────────────
#' @title Regressão Logística Binomial
#' @description Ajusta um modelo de regressão logística binomial, gerando um
#' painel completo com Deviance, LRT, Pseudo-R², AUC, Odds Ratio (OR),
#' matriz de confusão, teste de Hosmer-Lemeshow e interpretação automática.
#' @param formula Fórmula do modelo (ex: y ~ x1 + x2).
#' @param dados Data frame contendo os dados.
#' @param corte_prob Ponto de corte para classificação binária. Padrão 0.5.
#' @param grafico Se TRUE, exibe painel com Curva ROC e Gráfico de Calibração.
#' @return Retorna invisivelmente o modelo ajustado (objeto glm).
#' @export
regressao_logistica <- function(formula, dados, corte_prob = 0.5, grafico = TRUE) {
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")

  y_name  <- as.character(formula[[2]])
  x_names <- attr(terms(formula, data = dados), "term.labels")

  if (!(y_name %in% names(dados))) stop("Variável resposta não encontrada.")

  Y_raw  <- na.omit(dados[[y_name]])
  niveis <- sort(unique(as.character(Y_raw)))

  if (length(niveis) > 2)
    stop(sprintf(
      "A variável resposta '%s' possui %d níveis. A regressão multinomial será adicionada em breve! Para a regressão logística binomial, a resposta deve ter exatamente 2 categorias.",
      y_name, length(niveis)
    ))
  if (length(niveis) < 2)
    stop("A variável resposta deve possuir 2 categorias (variabilidade zero detectada).")

  mod <- tryCatch(
    glm(formula, data = dados, family = binomial(link = "logit")),
    error = function(e) stop("Erro ao ajustar o modelo logístico: ", e$message)
  )

  df_mod <- mod$model
  Y_fit  <- df_mod[[1]]

  if (is.factor(Y_fit)) {
    sucesso_label <- levels(Y_fit)[2]
    Y_num <- as.integer(Y_fit == sucesso_label)
  } else {
    vals  <- sort(unique(as.numeric(Y_fit)))
    sucesso_label <- as.character(vals[2])
    Y_num <- as.integer(as.numeric(Y_fit) == vals[2])
  }

  probs     <- mod$fitted.values
  n_obs     <- length(Y_num)
  pred_cls  <- as.integer(probs >= corte_prob)

  # ── Métricas globais ───────────────────────────────────────────────────────
  dev_nula  <- mod$null.deviance
  dev_res   <- mod$deviance
  df_nula   <- mod$df.null
  df_res_df <- mod$df.residual
  lrt_chi   <- dev_nula - dev_res
  lrt_df    <- df_nula  - df_res_df
  lrt_p     <- 1 - pchisq(lrt_chi, lrt_df)
  pseudo_r2 <- 1 - dev_res / dev_nula   # McFadden
  aic_val   <- AIC(mod)
  bic_val   <- BIC(mod)
  K_params  <- length(coef(mod)) + 1
  aicc_val  <- if (n_obs - K_params - 1 > 0)
    aic_val + (2 * K_params * (K_params + 1)) / (n_obs - K_params - 1)
  else Inf

  auc_val   <- .calc_auc_estatR(probs, Y_num)
  vp <- sum(pred_cls == 1 & Y_num == 1)
  vn <- sum(pred_cls == 0 & Y_num == 0)
  fp <- sum(pred_cls == 1 & Y_num == 0)
  fn <- sum(pred_cls == 0 & Y_num == 1)
  acc <- (vp + vn) / n_obs
  sensib <- if ((vp + fn) > 0) vp / (vp + fn) else NA
  especif <- if ((vn + fp) > 0) vn / (vn + fp) else NA

  hl <- .hosmer_lemeshow_estatR(Y_num, probs, g = 10)

  # ── Coeficientes ──────────────────────────────────────────────────────────
  sm    <- summary(mod)
  coefs <- sm$coefficients
  ci    <- suppressMessages(confint.default(mod))
  rn    <- rownames(coefs)
  rn[rn == "(Intercept)"] <- "(Intercepto)"

  # ── Larguras da tabela de coeficientes ────────────────────────────────────
  tem_multiplos <- length(x_names) > 1
  w_c <- c(var = max(12, max(nchar(rn)) + 2),
            est = 14, ep = 12, z = 10, p = 16,
            or = 10, ic_inf = 14, ic_sup = 14)

  hdr_c <- paste0(
    .pad_log("Variável",         w_c["var"],    "left"),
    .pad_log("Estimativa",       w_c["est"],    "center"),
    .pad_log("Erro Pad.",        w_c["ep"],     "center"),
    .pad_log("Z",                w_c["z"],      "center"),
    .pad_log("p-valor",          w_c["p"],      "center"),
    .pad_log("OR",               w_c["or"],     "center"),
    .pad_log("IC(OR) 95% Inf",   w_c["ic_inf"], "center"),
    .pad_log("IC(OR) 95% Sup",   w_c["ic_sup"], "center")
  )
  w_coef <- nchar(hdr_c)

  # ── Impressão ─────────────────────────────────────────────────────────────
  .print_titulo("REGRESSÃO LOGÍSTICA BINOMIAL")
  cat(sprintf("  Variável Resposta:   %s (Sucesso = '%s')\n", y_name, sucesso_label))
  cat(sprintf("  Observações Válidas: %d\n\n", n_obs))

  # ── RESUMO GERAL (AJUSTE) ─────────────────────────────────────────────────
  .print_topico("RESUMO GERAL (AJUSTE)")
  w_aj <- c(med = 26, val = 30)
  aj_sep <- .sep_log(w_aj["med"] + w_aj["val"])
  aj_row <- function(m, v) {
    cat("  ", .pad_log(m, w_aj["med"], "left"), .pad_log(v, w_aj["val"], "right"), "\n", sep = "")
  }
  cat(aj_sep, "\n")
  cat("  ", .pad_log("Métrica", w_aj["med"], "left"), .pad_log("Valor", w_aj["val"], "right"), "\n", sep = "")
  cat(aj_sep, "\n")
  aj_row("Deviance Nula",    sprintf("%.2f  (df = %d)", dev_nula,  df_nula))
  aj_row("Deviance Residual",sprintf("%.2f  (df = %d)", dev_res,   df_res_df))
  aj_row("Teste LRT (\u03c7\u00b2)",   sprintf("%.2f  (p %s)  %s", lrt_chi,
                                               if (lrt_p < 0.001) "< 0,001" else paste0("= ", formatC(lrt_p, format="f", digits=3, decimal.mark=",")),
                                               .asterisk_log(lrt_p)))
  aj_row("AIC",   sprintf("%.2f", aic_val))
  aj_row("AICc",  sprintf("%.2f", aicc_val))
  aj_row("BIC",   sprintf("%.2f", bic_val))
  cat(aj_sep, "\n")
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n\n")

  # ── MÉTRICAS DE PERFORMANCE ───────────────────────────────────────────────
  .print_topico("MÉTRICAS DE PERFORMANCE")
  cat(aj_sep, "\n")
  cat("  ", .pad_log("Métrica", w_aj["med"], "left"), .pad_log("Valor", w_aj["val"], "right"), "\n", sep = "")
  cat(aj_sep, "\n")
  aj_row("Pseudo-R\u00b2 (McFadden)", .fmt_pct(pseudo_r2, 1))
  aj_row("Acurácia Global",    sprintf("%s  (Corte: %.2f)", .fmt_pct(acc, 1), corte_prob))
  aj_row("Sensibilidade",      if (!is.na(sensib))  .fmt_pct(sensib, 1)  else "-")
  aj_row("Especificidade",     if (!is.na(especif)) .fmt_pct(especif, 1) else "-")
  aj_row("AUC (Curva ROC)",    sprintf("%.3f", auc_val))
  cat(aj_sep, "\n\n")

  # ── TESTE DE ADEQUAÇÃO ────────────────────────────────────────────────────
  .print_topico("TESTE DE ADEQUAÇÃO (Pressuposto)")
  w_hl <- c(teste = 30, est = 20, p = 16)
  hl_sep <- .sep_log(w_hl["teste"] + w_hl["est"] + w_hl["p"])
  cat(hl_sep, "\n")
  cat("  ",
      .pad_log("Teste",            w_hl["teste"], "left"),
      .pad_log("Estatística (\u03c7\u00b2)", w_hl["est"],   "center"),
      .pad_log("p-valor",          w_hl["p"],     "center"),
      "\n", sep = "")
  cat(hl_sep, "\n")
  if (!is.na(hl$statistic)) {
    cat("  ",
        .pad_log("Hosmer-Lemeshow (g=10)", w_hl["teste"], "left"),
        .pad_log(sprintf("%.2f", hl$statistic), w_hl["est"], "center"),
        .pad_log(.formata_p_local(hl$p.value),  w_hl["p"],   "center"),
        "\n", sep = "")
  } else {
    cat("  ", .pad_log("Hosmer-Lemeshow", w_hl["teste"], "left"),
        .pad_log("-", w_hl["est"] + w_hl["p"], "center"), "\n", sep = "")
  }
  cat(hl_sep, "\n\n")

  # ── MATRIZ DE CONFUSÃO ────────────────────────────────────────────────────
  .print_topico("MATRIZ DE CONFUSÃO")
  lbl0 <- niveis[1]; lbl1 <- niveis[2]
  w_mx <- c(rot = 18, c0 = 14, c1 = 14)
  mx_sep <- .sep_log(w_mx["rot"] + w_mx["c0"] + w_mx["c1"])
  cat(mx_sep, "\n")
  cat("  ",
      .pad_log("", w_mx["rot"], "left"),
      .pad_log(paste0("Previsto: ", lbl0), w_mx["c0"], "center"),
      .pad_log(paste0("Previsto: ", lbl1), w_mx["c1"], "center"),
      "\n", sep = "")
  cat(mx_sep, "\n")
  cat("  ",
      .pad_log(paste0("Real: ", lbl0), w_mx["rot"], "left"),
      .pad_log(as.character(vn), w_mx["c0"], "center"),
      .pad_log(as.character(fp), w_mx["c1"], "center"),
      "\n", sep = "")
  cat("  ",
      .pad_log(paste0("Real: ", lbl1), w_mx["rot"], "left"),
      .pad_log(as.character(fn), w_mx["c0"], "center"),
      .pad_log(as.character(vp), w_mx["c1"], "center"),
      "\n", sep = "")
  cat(mx_sep, "\n\n")

  # ── COEFICIENTES E ODDS RATIO ─────────────────────────────────────────────
  .print_topico("COEFICIENTES E ODDS RATIO (OR)")
  cat(.sep_log(w_coef), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.sep_log(w_coef), "\n")

  for (i in seq_len(nrow(coefs))) {
    p_ast  <- paste(.formata_p_local(coefs[i, 4]), .asterisk_log(coefs[i, 4]))
    or_val <- exp(coefs[i, 1])
    cat("  ",
        .pad_log(rn[i],                                     w_c["var"],    "left"),
        .pad_log(sprintf("%.4f", coefs[i, 1]),              w_c["est"],    "center"),
        .pad_log(sprintf("%.4f", coefs[i, 2]),              w_c["ep"],     "center"),
        .pad_log(sprintf("%.2f",  coefs[i, 3]),             w_c["z"],      "center"),
        .pad_log(p_ast,                                     w_c["p"],      "center"),
        .pad_log(.fmt_or(or_val),                       w_c["or"],     "center"),
        .pad_log(.fmt_or(exp(ci[i, 1])),                w_c["ic_inf"], "center"),
        .pad_log(.fmt_or(exp(ci[i, 2])),                w_c["ic_sup"], "center"),
        "\n", sep = "")
  }
  cat(.sep_log(w_coef), "\n")
  cat("  Notas:\n")
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n")
  cat("  OR > 1: aumenta as chances de sucesso. OR < 1: reduz as chances.\n\n")

  # ── INTERPRETAÇÃO ─────────────────────────────────────────────────────────
  .print_topico("INTERPRETAÇÃO")

  preds_str <- paste(x_names, collapse = " e ")
  cat(sprintf("  O modelo de regressão logística binomial foi ajustado para\n  classificar %s a partir de %s.\n",
              y_name, preds_str))

  lrt_valido <- if (lrt_p < 0.05) "é globalmente válido" else "não é globalmente válido"
  lrt_p_txt  <- if (lrt_p < 0.001) "< 0,001" else paste0("= ", formatC(lrt_p, format="f", digits=3, decimal.mark=","))
  cat(sprintf("\n  O modelo %s (LRT p %s) e apresenta um\n  Pseudo-R² de McFadden de %s, indicando %s\n  da variabilidade do desfecho explicada pelo modelo.\n",
              lrt_valido, lrt_p_txt, .fmt_pct(pseudo_r2, 1),
              if (pseudo_r2 >= 0.20) "um ajuste satisfatório" else "um ajuste modesto"))

  cat(sprintf("\n  Em termos de discriminação, a AUC foi de %.3f, indicando %s\n  capacidade do modelo em separar os grupos.\n",
              auc_val,
              if      (auc_val >= 0.90) "excelente"
              else if (auc_val >= 0.80) "boa"
              else if (auc_val >= 0.70) "aceitável"
              else                      "baixa"))

  if (!is.na(hl$p.value)) {
    hl_interp <- if (hl$p.value > 0.05)
      "As probabilidades preditas estão bem calibradas com os eventos observados."
    else
      "As probabilidades preditas desviam significativamente dos eventos observados, sugerindo ajuste inadequado."
    cat(sprintf("\n  Quanto à calibração (Hosmer-Lemeshow, g=10): p = %s.\n  %s\n",
                .formata_p_local(hl$p.value), hl_interp))
  }

  # Coeficientes significativos
  for (i in seq_along(x_names)) {
    idx   <- i + 1  # pula intercepto
    if (idx > nrow(coefs)) next
    p_val <- coefs[idx, 4]
    if (p_val < 0.10) {
      or_val    <- exp(coefs[idx, 1])
      direcao   <- if (or_val > 1) "aumenta" else "reduz"
      pct_or    <- abs(or_val - 1) * 100
      mantendo  <- if (tem_multiplos) ", mantendo os demais preditores constantes" else ""
      p_texto   <- if (p_val < 0.001) "p < 0,001" else
        sprintf("p = %s", formatC(p_val, format="f", digits=3, decimal.mark=","))
      cat(sprintf("\n  O preditor '%s' apresentou efeito estatisticamente\n  significativo (%s)%s.\n  Um aumento de uma unidade em '%s' %s as chances de '%s'\n  em %.1f%% (OR = %.3f; IC 95%%: [%.3f; %.3f]).\n",
                  x_names[i], p_texto, mantendo,
                  x_names[i], direcao, sucesso_label,
                  pct_or, or_val,
                  exp(ci[idx, 1]), exp(ci[idx, 2])))
    }
  }

  cat("\n")

  .print_rodape()

  # ── PAINEL GRÁFICO (ROC + Calibração) ─────────────────────────────────────
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tema_painel <- if (exists("meu_tema")) meu_tema(grade = "dupla") else ggplot2::theme_minimal(base_size = 11)

    # Plot 1: Curva ROC
    df_roc <- .calc_roc_df(probs, Y_num)
    p1 <- ggplot2::ggplot(df_roc, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 1.2) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::labs(
        title    = "Curva ROC",
        subtitle = sprintf("AUC = %.3f", auc_val),
        x        = "Taxa de Falsos Positivos",
        y        = "Taxa de Verdadeiros Positivos"
      ) +
      tema_painel

    # Plot 2: Calibração (Observado vs Predito por decil)
    prob_breaks <- quantile(probs, probs = seq(0, 1, by = 0.1))
    prob_breaks[1] <- -Inf; prob_breaks[11] <- Inf
    grp     <- cut(probs, breaks = prob_breaks, labels = FALSE, include.lowest = TRUE)
    df_cal  <- data.frame(grupo = grp, Y = Y_num, prob = probs)
    df_cal_agg <- do.call(rbind, lapply(1:10, function(g) {
      sub <- df_cal[df_cal$grupo == g, ]
      if (nrow(sub) == 0) return(NULL)
      data.frame(obs = mean(sub$Y), pred = mean(sub$prob))
    }))
    df_cal_agg <- df_cal_agg[!is.na(df_cal_agg$obs), ]

    p2 <- ggplot2::ggplot(df_cal_agg, ggplot2::aes(x = pred, y = obs)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      ggplot2::geom_point(color = "#555555", size = 3, alpha = 0.85) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429",
                           linewidth = 0.8, formula = y ~ x) +
      ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::labs(
        title    = "Gráfico de Calibração",
        subtitle = sprintf("Hosmer-Lemeshow: p = %s", .formata_p_local(hl$p.value)),
        x        = "Probabilidade Prevista Média",
        y        = "Proporção Observada",
        caption  = "estatR"
      ) +
      tema_painel

    if (requireNamespace("patchwork", quietly = TRUE)) {
      painel <- p1 | p2
      suppressMessages(suppressWarnings(print(painel)))
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      suppressMessages(suppressWarnings(gridExtra::grid.arrange(p1, p2, ncol = 2)))
    } else {
      suppressMessages(suppressWarnings({ print(p1); print(p2) }))
    }
  }

  invisible(mod)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Métricas do Modelo Logístico
#' @description Exibe as métricas de ajuste e performance de um modelo logístico
#' ajustado por \code{regressao_logistica()}.
#' @param modelo Objeto \code{glm} retornado por \code{regressao_logistica()}.
#' @param corte_prob Ponto de corte para classificação. Padrão 0.5.
#' @return Retorna invisivelmente um data frame com as métricas.
#' @export
metricas_logistica <- function(modelo, corte_prob = 0.5) {
  if (!inherits(modelo, "glm")) stop("'modelo' deve ser um objeto glm (regressao_logistica).")

  Y_num  <- as.integer(modelo$y)
  probs  <- modelo$fitted.values
  n_obs  <- length(Y_num)
  K      <- length(coef(modelo)) + 1

  dev_nula  <- modelo$null.deviance
  dev_res   <- modelo$deviance
  pseudo_r2 <- 1 - dev_res / dev_nula
  aic_val   <- AIC(modelo)
  bic_val   <- BIC(modelo)
  aicc_val  <- if (n_obs - K - 1 > 0) aic_val + (2 * K * (K + 1)) / (n_obs - K - 1) else Inf

  auc_val  <- .calc_auc_estatR(probs, Y_num)
  pred_cls <- as.integer(probs >= corte_prob)
  vp <- sum(pred_cls == 1 & Y_num == 1); vn <- sum(pred_cls == 0 & Y_num == 0)
  fp <- sum(pred_cls == 1 & Y_num == 0); fn <- sum(pred_cls == 0 & Y_num == 1)
  acc      <- (vp + vn) / n_obs
  sensib   <- if ((vp + fn) > 0) vp / (vp + fn) else NA
  especif  <- if ((vn + fp) > 0) vn / (vn + fp) else NA
  ppv      <- if ((vp + fp) > 0) vp / (vp + fp) else NA
  npv      <- if ((vn + fn) > 0) vn / (vn + fn) else NA

  .print_titulo("MÉTRICAS DO MODELO LOGÍSTICO")

  lw <- 26; vw <- 12
  sep <- paste0(.sep_log(lw + vw), "\n")
  f3  <- function(x) if (is.na(x)) "-" else formatC(x, format = "f", digits = 3, decimal.mark = ",")

  cat(sep)
  cat("  ", .pad_log("Métrica", lw, "left"), .pad_log("Valor", vw, "right"), "\n", sep = "")
  cat(sep)
  cat("  ", .pad_log("Pseudo-R\u00b2 (McFadden)", lw, "left"), .pad_log(.fmt_pct(pseudo_r2, 2), vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("AIC",      lw, "left"), .pad_log(sprintf("%.2f", aic_val),  vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("AICc",     lw, "left"), .pad_log(sprintf("%.2f", aicc_val), vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("BIC",      lw, "left"), .pad_log(sprintf("%.2f", bic_val),  vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("AUC",      lw, "left"), .pad_log(f3(auc_val),               vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("Acurácia", lw, "left"), .pad_log(.fmt_pct(acc, 2),          vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("Sensibilidade",   lw, "left"), .pad_log(if (!is.na(sensib))  .fmt_pct(sensib, 2)  else "-", vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("Especificidade",  lw, "left"), .pad_log(if (!is.na(especif)) .fmt_pct(especif, 2) else "-", vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("VPP (Precisão)",  lw, "left"), .pad_log(if (!is.na(ppv))     .fmt_pct(ppv, 2)     else "-", vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("VPN",             lw, "left"), .pad_log(if (!is.na(npv))     .fmt_pct(npv, 2)     else "-", vw, "right"), "\n", sep = "")
  cat(sep)
  cat(sprintf("  * Corte de classificação: %.2f\n\n", corte_prob))

  .print_rodape()
  invisible(data.frame(
    Metrica = c("Pseudo_R2","AIC","AICc","BIC","AUC","Acuracia","Sensibilidade","Especificidade","VPP","VPN"),
    Valor   = c(pseudo_r2, aic_val, aicc_val, bic_val, auc_val, acc, sensib, especif, ppv, npv)
  ))
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Seleção de Modelos Logísticos (Best Subsets)
#' @description Avalia múltiplas combinações de variáveis preditoras para
#' regressão logística binomial, ordenando pelo AIC (critério padrão para GLMs).
#' @param formula Fórmula com todos os candidatos (ex: y ~ x1 + x2 + x3).
#' @param dados Data frame contendo os dados.
#' @param top Número de modelos a exibir no ranking. Padrão 5.
#' @return Retorna invisivelmente um data frame com o ranking completo.
#' @export
selecao_modelos_logistica <- function(formula, dados, top = 5) {
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")

  termos   <- attr(terms(formula, data = dados), "term.labels")
  y_name   <- as.character(formula[[2]])
  cols_req <- c(y_name, termos)

  df_clean <- na.omit(dados[, cols_req, drop = FALSE])
  n        <- nrow(df_clean)
  if (n < 5) stop("Poucas observações válidas após remover NAs (mínimo 5).")

  Y_raw  <- df_clean[[y_name]]
  niveis <- unique(as.character(Y_raw))
  if (length(niveis) != 2) stop("A variável resposta deve ter exatamente 2 categorias.")

  p_total <- length(termos)
  usou_heuristica <- FALSE
  termos_orig     <- termos

  # Screening se p > 10
  if (p_total > 10) {
    usou_heuristica <- TRUE
    Y_bin <- as.integer(factor(df_clean[[y_name]])) - 1
    num_t <- termos[sapply(df_clean[termos], is.numeric)]
    cors  <- sapply(num_t, function(v) abs(cor(Y_bin, df_clean[[v]], use = "complete.obs")))
    top10 <- names(sort(cors, decreasing = TRUE))[1:min(10, length(cors))]
    cat_t <- setdiff(termos, num_t)
    termos <- unique(c(top10, cat_t))[1:min(10, length(unique(c(top10, cat_t))))]
  }

  mods_list <- unlist(lapply(1:length(termos), function(k) combn(termos, k, simplify = FALSE)), recursive = FALSE)
  total_av  <- length(mods_list)

  res <- lapply(mods_list, function(preds) {
    f_tmp <- as.formula(paste(y_name, "~", paste(preds, collapse = " + ")))
    mod   <- tryCatch(glm(f_tmp, data = df_clean, family = binomial()), error = function(e) NULL)
    if (is.null(mod)) return(NULL)
    k_p   <- length(coef(mod)) + 1
    aic_v <- AIC(mod)
    aicc_v <- if (n - k_p - 1 > 0) aic_v + (2 * k_p * (k_p + 1)) / (n - k_p - 1) else Inf
    probs  <- mod$fitted.values
    Y_num  <- as.integer(mod$y)
    auc_v  <- .calc_auc_estatR(probs, Y_num)
    list(modelo = paste(preds, collapse = " + "), k = length(preds),
         AIC = aic_v, AICc = aicc_v, AUC = auc_v)
  })

  res     <- Filter(Negate(is.null), res)
  df_res  <- do.call(rbind, lapply(res, as.data.frame))
  df_res  <- df_res[order(df_res$AIC), ]

  df_top <- if (nrow(df_res) > top) df_res[1:top, ] else df_res

  # Impressão
  .print_titulo("SELEÇÃO DE MODELOS — REGRESSÃO LOGÍSTICA")
  cat(sprintf("  Variável Resposta: %s\n", y_name))
  if (usou_heuristica) {
    cat(sprintf("  Preditores base:   %s\n", paste(termos_orig, collapse = ", ")))
    cat(sprintf("  Método utilizado:  Heurístico (Screening TOP 10 → %d modelos avaliados)\n", total_av))
  } else {
    cat(sprintf("  Preditores base:   %s\n", paste(termos, collapse = ", ")))
    cat(sprintf("  Método utilizado:  Exaustivo (%d modelos avaliados)\n", total_av))
  }
  cat("\n")

  .print_topico(sprintf("TOP %d MODELOS (Ordenados por AIC)", nrow(df_top)))

  df_fmt <- data.frame(
    Ranking = paste0(1:nrow(df_top), "\u00ba"),
    Modelo  = df_top$modelo,
    k       = as.character(df_top$k),
    AIC     = sapply(df_top$AIC,  .fmt_num, decimais = 1),
    AICc    = sapply(df_top$AICc, .fmt_num, decimais = 1),
    AUC     = sapply(df_top$AUC,  function(v) if (is.na(v)) "-" else sprintf("%.3f", v)),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_fmt, align = c("left", "left", "center", "center", "center", "center"))
  } else {
    print(df_fmt, row.names = FALSE)
  }

  cat("  * k = Número de variáveis preditoras.\n")
  cat("  * AIC menor = modelo mais parcimonioso e melhor ajustado.\n\n")
  .print_rodape()

  invisible(df_res)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Análise de Resíduos do Modelo Logístico
#' @description Realiza diagnóstico dos resíduos de um modelo logístico:
#' testa a qualidade do ajuste por Hosmer-Lemeshow, calcula resíduos de Pearson
#' e de desvio, identifica pontos influentes e plota um painel gráfico de diagnóstico.
#' @param modelo Objeto \code{glm} retornado por \code{regressao_logistica()}.
#' @param grafico Se TRUE, exibe o painel 2x2 de diagnóstico. Padrão TRUE.
#' @return Retorna invisivelmente uma lista com os diagnósticos.
#' @export
analise_residual_logistica <- function(modelo, grafico = TRUE) {
  if (!inherits(modelo, "glm")) stop("'modelo' deve ser um objeto glm.")

  Y_num   <- as.integer(modelo$y)
  probs   <- modelo$fitted.values
  n       <- length(Y_num)

  # Resíduos
  res_pearson <- residuals(modelo, type = "pearson")
  res_devian  <- residuals(modelo, type = "deviance")
  alavancagem <- hatvalues(modelo)
  dist_cook   <- cooks.distance(modelo)

  # Hosmer-Lemeshow
  hl <- .hosmer_lemeshow_estatR(Y_num, probs, g = 10)

  # Dispersão de Pearson
  chi_pearson <- sum(res_pearson^2)
  df_pearson  <- n - length(coef(modelo))
  disp_pearson <- chi_pearson / df_pearson

  formata_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }

  # ── Impressão ──────────────────────────────────────────────────────────────
  .print_titulo("ANÁLISE DE RESÍDUOS — REGRESSÃO LOGÍSTICA")

  # Tabela de diagnóstico
  w_t <- c(teste = 34, est = 20, p = 16)
  sep_t <- .sep_log(w_t["teste"] + w_t["est"] + w_t["p"])
  .print_topico("DIAGNÓSTICO DO AJUSTE")
  cat(sep_t, "\n")
  cat("  ",
      .pad_log("Diagnóstico",       w_t["teste"], "left"),
      .pad_log("Valor",             w_t["est"],   "center"),
      .pad_log("p-valor",           w_t["p"],     "center"),
      "\n", sep = "")
  cat(sep_t, "\n")
  if (!is.na(hl$statistic)) {
    cat("  ",
        .pad_log("Hosmer-Lemeshow (g=10)", w_t["teste"], "left"),
        .pad_log(sprintf("\u03c7\u00b2 = %.2f", hl$statistic), w_t["est"], "center"),
        .pad_log(formata_p(hl$p.value), w_t["p"], "center"),
        "\n", sep = "")
  }
  cat("  ",
      .pad_log("Dispersão de Pearson", w_t["teste"], "left"),
      .pad_log(sprintf("%.3f", disp_pearson), w_t["est"], "center"),
      .pad_log("-", w_t["p"], "center"),
      "\n", sep = "")
  cat(sep_t, "\n\n")

  # Pontos influentes (Cook > 4/n)
  infl_idx <- which(dist_cook > 4 / n)
  .print_topico("PONTOS INFLUENTES")
  if (length(infl_idx) > 0) {
    cat(sprintf("  Observações com distância de Cook > %.4f (4/n):\n", 4/n))
    for (idx in infl_idx) {
      cat(sprintf("  Obs. %d: Cook = %.4f | Alavancagem = %.4f | Res. Pearson = %.3f\n",
                  idx, dist_cook[idx], alavancagem[idx], res_pearson[idx]))
    }
    cat("\n")
  } else {
    cat("  Nenhuma observação com distância de Cook > 4/n detectada.\n\n")
  }

  # Interpretação
  .print_topico("INTERPRETAÇÃO DO DIAGNÓSTICO")

  if (!is.na(hl$p.value)) {
    if (hl$p.value > 0.05) {
      cat(sprintf("  [\u2713] : Hosmer-Lemeshow não rejeitou H\u2080 (p = %s).\n  As probabilidades preditas estão bem calibradas com os eventos observados.\n\n",
                  formata_p(hl$p.value)))
    } else {
      cat(sprintf("  [!] : Hosmer-Lemeshow rejeitou H\u2080 (p = %s).\n  As probabilidades preditas desviam dos eventos observados, indicando\n  possível problema de especificação do modelo.\n\n",
                  formata_p(hl$p.value)))
    }
  }

  if (disp_pearson > 1.5) {
    cat(sprintf("  [!] : Dispersão de Pearson = %.3f (> 1,5). Pode indicar\n  sobredispersão ou pontos influentes afetando o ajuste.\n\n",
                disp_pearson))
  } else {
    cat(sprintf("  [\u2713] : Dispersão de Pearson = %.3f. O modelo não apresenta\n  sinais de sobredispersão relevante.\n\n",
                disp_pearson))
  }

  if (length(infl_idx) > 0) {
    cat(sprintf("  [!] : %d observação(ões) com alta influência detectada(s).\n  Considere investigar ou remover essas observações e reajustar o modelo.\n\n",
                length(infl_idx)))
  }

  cat(.sep_log(w_t["teste"] + w_t["est"] + w_t["p"]), "\n\n")
  .print_rodape()

  # ── Painel Gráfico ─────────────────────────────────────────────────────────
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tema_painel <- if (exists("meu_tema")) meu_tema(grade = "dupla") else ggplot2::theme_minimal(base_size = 11)

    df_diag <- data.frame(
      idx        = seq_len(n),
      ajustados  = probs,
      res_pearson = res_pearson,
      res_devian  = res_devian,
      cook        = dist_cook,
      alavancagem = alavancagem
    )

    # P1: Resíduos de Pearson vs Predito
    p1 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = ajustados, y = res_pearson)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "gray70") +
      ggplot2::geom_point(color = "#555555", alpha = 0.8) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429",
                           linewidth = 0.7, formula = y ~ x) +
      ggplot2::labs(title = "Resíduos de Pearson vs Predito",
                    x = "Probabilidade Predita", y = "Resíduo de Pearson") +
      tema_painel

    # P2: Resíduos de Desvio (QQ-Plot)
    p2 <- ggplot2::ggplot(df_diag, ggplot2::aes(sample = res_devian)) +
      ggplot2::stat_qq(color = "#555555", alpha = 0.8) +
      ggplot2::stat_qq_line(color = "#D90429", linewidth = 0.7) +
      ggplot2::labs(title = "QQ-Plot (Resíduos de Desvio)",
                    x = "Quantis Teóricos", y = "Quantis Amostrais") +
      tema_painel

    # P3: Distância de Cook
    p3 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = idx, y = cook)) +
      ggplot2::geom_hline(yintercept = 4/n, linetype = "dashed", color = "#D90429") +
      ggplot2::geom_segment(ggplot2::aes(xend = idx, yend = 0), color = "#555555", alpha = 0.7) +
      ggplot2::geom_point(color = "#555555", alpha = 0.8) +
      ggplot2::labs(title = "Distância de Cook",
                    x = "Observação", y = "Cook") +
      tema_painel

    # P4: Alavancagem vs Resíduo Padronizado
    p4 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = alavancagem, y = res_pearson)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      ggplot2::geom_point(color = "#555555", alpha = 0.8) +
      ggplot2::labs(title = "Alavancagem vs Resíduo de Pearson",
                    x = "Alavancagem (Leverage)", y = "Resíduo de Pearson",
                    caption = "estatR") +
      tema_painel

    if (requireNamespace("patchwork", quietly = TRUE)) {
      painel <- (p1 | p2) / (p3 | p4)
      suppressMessages(suppressWarnings(print(painel)))
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      suppressMessages(suppressWarnings(gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)))
    } else {
      suppressMessages(suppressWarnings({ print(p1); print(p2); print(p3); print(p4) }))
    }
  }

  invisible(list(
    hosmer_lemeshow = hl,
    dispersao_pearson = disp_pearson,
    pontos_influentes = infl_idx
  ))
}
