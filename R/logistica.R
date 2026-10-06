# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: REGRESSÃO LOGÍSTICA (BINOMIAL · MULTINOMIAL · ORDINAL)
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

.tema_log <- function() {
  if (exists("tema_estatR", mode = "function"))
    return(tema_estatR(estilo = 2))
  if (requireNamespace("estatR", quietly = TRUE)) {
    fn <- get("tema_estatR", envir = asNamespace("estatR"))
    return(fn(estilo = 2))
  }
  ggplot2::theme_minimal(base_size = 11)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Regressão Logística (Binomial · Multinomial · Ordinal)
#' @description Ajusta um modelo de regressão logística, detectando automaticamente
#' se a variável resposta é binária (2 categorias) ou policotômica (3+ categorias).
#' Para Y policotômico, use \code{tipo = "nominal"} (multinomial via \code{nnet::multinom})
#' ou \code{tipo = "ordinal"} (logística proporcional via \code{MASS::polr}).
#'
#' @param formula Fórmula do modelo (ex: y ~ x1 + x2).
#' @param dados Data frame contendo os dados.
#' @param tipo Tipo para Y policotômico: \code{"nominal"} (padrão) ou \code{"ordinal"}.
#'   Ignorado quando Y é binário.
#' @param corte_prob Ponto de corte para classificação binária. Padrão 0.5.
#' @param grafico Se TRUE, exibe painel gráfico diagnóstico. Padrão TRUE.
#' @return Retorna invisivelmente o modelo ajustado.
#' @export
regressao_logistica <- function(formula, dados, tipo = c("nominal", "ordinal"),
                                corte_prob = 0.5, grafico = TRUE) {
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  tipo <- match.arg(tipo)

  y_name  <- as.character(formula[[2]])
  x_names <- attr(terms(formula, data = dados), "term.labels")

  if (!(y_name %in% names(dados))) stop("Variável resposta não encontrada.")

  Y_raw  <- na.omit(dados[[y_name]])
  niveis <- sort(unique(as.character(Y_raw)))
  n_niv  <- length(niveis)

  if (n_niv < 2) stop("A variável resposta deve possuir ao menos 2 categorias.")

  # ── Despacho por número de níveis ──────────────────────────────────────────
  if (n_niv == 2) {
    return(.logistica_binomial(formula, dados, corte_prob, grafico,
                               y_name, x_names, niveis))
  } else {
    if (tipo == "ordinal") {
      return(.logistica_ordinal(formula, dados, grafico, y_name, x_names, niveis))
    } else {
      return(.logistica_multinomial(formula, dados, grafico, y_name, x_names, niveis))
    }
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# BINOMIAL (código original preservado integralmente)
# ─────────────────────────────────────────────────────────────────────────────
.logistica_binomial <- function(formula, dados, corte_prob, grafico,
                                y_name, x_names, niveis) {
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

  dev_nula  <- mod$null.deviance
  dev_res   <- mod$deviance
  df_nula   <- mod$df.null
  df_res_df <- mod$df.residual
  lrt_chi   <- dev_nula - dev_res
  lrt_df    <- df_nula  - df_res_df
  lrt_p     <- 1 - pchisq(lrt_chi, lrt_df)
  pseudo_r2 <- 1 - dev_res / dev_nula
  aic_val   <- AIC(mod)
  bic_val   <- BIC(mod)
  K_params  <- length(coef(mod)) + 1
  aicc_val  <- if (n_obs - K_params - 1 > 0)
    aic_val + (2 * K_params * (K_params + 1)) / (n_obs - K_params - 1)
  else Inf

  auc_val   <- .calc_auc_estatR(probs, Y_num)
  vp <- sum(pred_cls == 1 & Y_num == 1); vn <- sum(pred_cls == 0 & Y_num == 0)
  fp <- sum(pred_cls == 1 & Y_num == 0); fn <- sum(pred_cls == 0 & Y_num == 1)
  acc    <- (vp + vn) / n_obs
  sensib <- if ((vp + fn) > 0) vp / (vp + fn) else NA
  especif <- if ((vn + fp) > 0) vn / (vn + fp) else NA

  hl    <- .hosmer_lemeshow_estatR(Y_num, probs, g = 10)
  sm    <- summary(mod)
  coefs <- sm$coefficients
  ci    <- suppressMessages(confint.default(mod))
  rn    <- rownames(coefs)
  rn[rn == "(Intercept)"] <- "(Intercepto)"

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

  .print_titulo("REGRESSÃO LOGÍSTICA BINOMIAL")
  cat(sprintf("  Variável Resposta:   %s (Sucesso = '%s')\n", y_name, sucesso_label))
  cat(sprintf("  Observações Válidas: %d\n\n", n_obs))

  .print_topico("RESUMO GERAL (AJUSTE)")
  w_aj <- c(med = 26, val = 30)
  aj_sep <- .sep_log(w_aj["med"] + w_aj["val"])
  aj_row <- function(m, v) {
    cat("  ", .pad_log(m, w_aj["med"], "left"), .pad_log(v, w_aj["val"], "right"), "\n", sep = "")
  }
  cat(aj_sep, "\n")
  cat("  ", .pad_log("Métrica", w_aj["med"], "left"), .pad_log("Valor", w_aj["val"], "right"), "\n", sep = "")
  cat(aj_sep, "\n")
  aj_row("Deviance Nula",    sprintf("%.2f  (df = %d)", dev_nula, df_nula))
  aj_row("Deviance Residual",sprintf("%.2f  (df = %d)", dev_res,  df_res_df))
  aj_row("Teste LRT (\u03c7\u00b2)",   sprintf("%.2f  (p %s)  %s", lrt_chi,
                                              if (lrt_p < 0.001) "< 0,001" else paste0("= ", formatC(lrt_p, format="f", digits=3, decimal.mark=",")),
                                              .asterisk_log(lrt_p)))
  aj_row("AIC",   sprintf("%.2f", aic_val))
  aj_row("AICc",  sprintf("%.2f", aicc_val))
  aj_row("BIC",   sprintf("%.2f", bic_val))
  cat(aj_sep, "\n")
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n\n")

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

  .print_topico("MATRIZ DE CONFUSÃO")
  lbl0 <- niveis[1]; lbl1 <- niveis[2]
  w_mx <- c(rot = 18, c0 = 14, c1 = 14)
  mx_sep <- .sep_log(w_mx["rot"] + w_mx["c0"] + w_mx["c1"])
  cat(mx_sep, "\n")
  cat("  ", .pad_log("", w_mx["rot"], "left"),
      .pad_log(paste0("Previsto: ", lbl0), w_mx["c0"], "center"),
      .pad_log(paste0("Previsto: ", lbl1), w_mx["c1"], "center"),
      "\n", sep = "")
  cat(mx_sep, "\n")
  cat("  ", .pad_log(paste0("Real: ", lbl0), w_mx["rot"], "left"),
      .pad_log(as.character(vn), w_mx["c0"], "center"),
      .pad_log(as.character(fp), w_mx["c1"], "center"),
      "\n", sep = "")
  cat("  ", .pad_log(paste0("Real: ", lbl1), w_mx["rot"], "left"),
      .pad_log(as.character(fn), w_mx["c0"], "center"),
      .pad_log(as.character(vp), w_mx["c1"], "center"),
      "\n", sep = "")
  cat(mx_sep, "\n\n")

  .print_topico("COEFICIENTES E ODDS RATIO (OR)")
  cat(.sep_log(w_coef), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.sep_log(w_coef), "\n")
  for (i in seq_len(nrow(coefs))) {
    p_ast  <- paste(.formata_p_local(coefs[i, 4]), .asterisk_log(coefs[i, 4]))
    or_val <- exp(coefs[i, 1])
    cat("  ",
        .pad_log(rn[i],                    w_c["var"],    "left"),
        .pad_log(sprintf("%.4f", coefs[i, 1]), w_c["est"],    "center"),
        .pad_log(sprintf("%.4f", coefs[i, 2]), w_c["ep"],     "center"),
        .pad_log(sprintf("%.2f",  coefs[i, 3]), w_c["z"],      "center"),
        .pad_log(p_ast,                    w_c["p"],      "center"),
        .pad_log(.fmt_or(or_val),          w_c["or"],     "center"),
        .pad_log(.fmt_or(exp(ci[i, 1])),   w_c["ic_inf"], "center"),
        .pad_log(.fmt_or(exp(ci[i, 2])),   w_c["ic_sup"], "center"),
        "\n", sep = "")
  }
  cat(.sep_log(w_coef), "\n")
  cat("  Notas:\n")
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n")
  cat("  OR > 1: aumenta as chances de sucesso. OR < 1: reduz as chances.\n\n")

  .print_topico("INTERPRETAÇÃO")
  preds_str <- paste(x_names, collapse = " e ")
  cat(sprintf("  O modelo de regressão logística binomial foi ajustado para\n  classificar %s a partir de %s.\n", y_name, preds_str))
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
  for (i in seq_along(x_names)) {
    idx <- i + 1
    if (idx > nrow(coefs)) next
    p_val <- coefs[idx, 4]
    if (p_val < 0.10) {
      or_val   <- exp(coefs[idx, 1])
      direcao  <- if (or_val > 1) "aumenta" else "reduz"
      pct_or   <- abs(or_val - 1) * 100
      mantendo <- if (tem_multiplos) ", mantendo os demais preditores constantes" else ""
      p_texto  <- if (p_val < 0.001) "p < 0,001" else
        sprintf("p = %s", formatC(p_val, format="f", digits=3, decimal.mark=","))
      cat(sprintf("\n  O preditor '%s' apresentou efeito estatisticamente\n  significativo (%s)%s.\n  Um aumento de uma unidade em '%s' %s as chances de '%s'\n  em %.1f%% (OR = %.3f; IC 95%%: [%.3f; %.3f]).\n",
                  x_names[i], p_texto, mantendo,
                  x_names[i], direcao, sucesso_label,
                  pct_or, or_val, exp(ci[idx, 1]), exp(ci[idx, 2])))
    }
  }
  cat("\n")
  .print_rodape()

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tema_painel <- .tema_log()
    df_roc <- .calc_roc_df(probs, Y_num)
    p1 <- ggplot2::ggplot(df_roc, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 1.2) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::labs(title = "Curva ROC", subtitle = sprintf("AUC = %.3f", auc_val),
                    x = "Taxa de Falsos Positivos", y = "Taxa de Verdadeiros Positivos") +
      tema_painel
    prob_breaks <- quantile(probs, probs = seq(0, 1, by = 0.1))
    prob_breaks[1] <- -Inf; prob_breaks[11] <- Inf
    grp    <- cut(probs, breaks = prob_breaks, labels = FALSE, include.lowest = TRUE)
    df_cal <- data.frame(grupo = grp, Y = Y_num, prob = probs)
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
      ggplot2::labs(title = "Gráfico de Calibração",
                    subtitle = sprintf("Hosmer-Lemeshow: p = %s", .formata_p_local(hl$p.value)),
                    x = "Probabilidade Prevista Média", y = "Proporção Observada",
                    caption = "estatR") +
      tema_painel
    if (requireNamespace("patchwork", quietly = TRUE)) {
      suppressMessages(suppressWarnings(print(p1 | p2)))
    } else {
      suppressMessages(suppressWarnings({ print(p1); print(p2) }))
    }
  }
  invisible(mod)
}

# ─────────────────────────────────────────────────────────────────────────────
# MULTINOMIAL (nnet::multinom)
# ─────────────────────────────────────────────────────────────────────────────
.logistica_multinomial <- function(frm, dados, grafico, y_name, x_names, niveis) {
  if (!requireNamespace("nnet", quietly = TRUE))
    stop("Instale o pacote 'nnet' para regressão multinomial: install.packages('nnet')")

  mod <- tryCatch(
    do.call(nnet::multinom, list(formula = frm, data = dados, trace = FALSE)),
    error = function(e) stop("Erro ao ajustar modelo multinomial: ", e$message)
  )

  df_mod <- na.omit(dados[, c(y_name, x_names), drop = FALSE])
  n_obs  <- nrow(df_mod)
  Y_fac  <- factor(df_mod[[y_name]])
  ref    <- levels(Y_fac)[1]
  cats   <- levels(Y_fac)[-1]

  sm      <- summary(mod)
  coefs   <- sm$coefficients          # cats x (1 + p)
  se_mat  <- sm$standard.errors
  z_mat   <- coefs / se_mat
  p_mat   <- 2 * pnorm(-abs(z_mat))

  # ── Métricas globais ───────────────────────────────────────────────────────
  dev_res   <- mod$deviance
  dev_nula  <- -2 * sum(log(table(Y_fac) / n_obs)[as.character(Y_fac)])
  lrt_chi   <- dev_nula - dev_res
  lrt_df    <- (length(cats)) * length(x_names)
  lrt_p     <- 1 - pchisq(lrt_chi, lrt_df)
  pseudo_r2 <- 1 - dev_res / dev_nula
  aic_val   <- AIC(mod)
  bic_val   <- BIC(mod)
  K         <- length(coef(mod)) + 1
  aicc_val  <- if (n_obs - K - 1 > 0) aic_val + (2*K*(K+1))/(n_obs-K-1) else Inf

  pred_cls  <- as.character(predict(mod, type = "class"))
  acc       <- mean(pred_cls == as.character(df_mod[[y_name]]))

  # ── Impressão ─────────────────────────────────────────────────────────────
  .print_titulo("REGRESSÃO LOGÍSTICA MULTINOMIAL")
  cat(sprintf("  Variável Resposta:   %s  (%d categorias)\n", y_name, length(cats) + 1))
  cat(sprintf("  Categoria Referência: '%s'\n", ref))
  cat(sprintf("  Observações Válidas: %d\n\n", n_obs))

  .print_topico("RESUMO GERAL (AJUSTE)")
  w_aj <- c(med = 26, val = 30)
  aj_sep <- .sep_log(w_aj["med"] + w_aj["val"])
  aj_row <- function(m, v) cat("  ", .pad_log(m, w_aj["med"], "left"), .pad_log(v, w_aj["val"], "right"), "\n", sep = "")
  cat(aj_sep, "\n")
  cat("  ", .pad_log("Métrica", w_aj["med"], "left"), .pad_log("Valor", w_aj["val"], "right"), "\n", sep = "")
  cat(aj_sep, "\n")
  aj_row("Deviance Nula",     sprintf("%.2f", dev_nula))
  aj_row("Deviance Residual", sprintf("%.2f  (df = %d)", dev_res, mod$edf))
  aj_row("Teste LRT (\u03c7\u00b2)",    sprintf("%.2f  (p %s)  %s", lrt_chi,
                                               if (lrt_p < 0.001) "< 0,001" else paste0("= ", formatC(lrt_p, format="f", digits=3, decimal.mark=",")),
                                               .asterisk_log(lrt_p)))
  aj_row("Pseudo-R\u00b2 (McFadden)", .fmt_pct(pseudo_r2, 1))
  aj_row("Acurácia Global",   .fmt_pct(acc, 1))
  aj_row("AIC",  sprintf("%.2f", aic_val))
  aj_row("AICc", sprintf("%.2f", aicc_val))
  aj_row("BIC",  sprintf("%.2f", bic_val))
  cat(aj_sep, "\n")
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n\n")

  # ── Coeficientes por categoria ─────────────────────────────────────────────
  vnames <- c("(Intercepto)", x_names)
  w_c <- c(var = max(14, max(nchar(vnames)) + 2),
            est = 12, ep = 12, z = 10, p = 16, or = 10)
  hdr <- paste0(
    .pad_log("Variável",   w_c["var"], "left"),
    .pad_log("Estimativa", w_c["est"], "center"),
    .pad_log("Erro Pad.",  w_c["ep"],  "center"),
    .pad_log("Z",          w_c["z"],   "center"),
    .pad_log("p-valor",    w_c["p"],   "center"),
    .pad_log("OR",         w_c["or"],  "center")
  )
  w_tot <- nchar(hdr)

  for (cat_i in cats) {
    .print_topico(sprintf("COEFICIENTES — '%s' vs. '%s' (referência)", cat_i, ref))
    cat(.sep_log(w_tot), "\n")
    cat("  ", hdr, "\n", sep = "")
    cat(.sep_log(w_tot), "\n")
    for (j in seq_along(vnames)) {
      vn_j  <- if (j == 1) "(Intercepto)" else x_names[j - 1]
      est_j <- coefs[cat_i, j]
      se_j  <- se_mat[cat_i, j]
      z_j   <- z_mat[cat_i, j]
      p_j   <- p_mat[cat_i, j]
      p_txt <- paste(.formata_p_local(p_j), .asterisk_log(p_j))
      cat("  ",
          .pad_log(vn_j,               w_c["var"], "left"),
          .pad_log(sprintf("%.4f", est_j), w_c["est"], "center"),
          .pad_log(sprintf("%.4f", se_j),  w_c["ep"],  "center"),
          .pad_log(sprintf("%.2f",  z_j),  w_c["z"],   "center"),
          .pad_log(p_txt,            w_c["p"],   "center"),
          .pad_log(.fmt_or(exp(est_j)), w_c["or"], "center"),
          "\n", sep = "")
    }
    cat(.sep_log(w_tot), "\n\n")
  }

  # ── Interpretação ──────────────────────────────────────────────────────────
  .print_topico("INTERPRETAÇÃO")
  cat(sprintf("  O modelo multinomial foi ajustado para classificar '%s' em %d categorias\n  a partir de: %s.\n",
              y_name, length(cats) + 1, paste(x_names, collapse = ", ")))
  cat(sprintf("\n  A categoria de referência é '%s'. Os coeficientes e OR de cada bloco\n  expressam o risco relativo de pertencer àquela categoria vs. '%s'.\n",
              ref, ref))
  cat(sprintf("\n  O modelo %s (LRT: \u03c7\u00b2 = %.2f, p %s),\n  com Pseudo-R\u00b2 de McFadden de %s e acurácia global de %s.\n",
              if (lrt_p < 0.05) "é globalmente significativo" else "não é globalmente significativo",
              lrt_chi,
              if (lrt_p < 0.001) "< 0,001" else paste0("= ", formatC(lrt_p, format="f", digits=3, decimal.mark=",")),
              .fmt_pct(pseudo_r2, 1), .fmt_pct(acc, 1)))
  cat("\n")
  .print_rodape()

  # ── Gráfico: Probabilidades Preditas por Grupo ─────────────────────────────
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tema_g <- .tema_log()
    probs_mat <- predict(mod, type = "probs")
    if (is.vector(probs_mat)) probs_mat <- as.matrix(probs_mat)

    # Só plota gráfico de probabilidades se houver variável contínua
    cont_vars <- x_names[sapply(df_mod[x_names], is.numeric)]

    if (length(cont_vars) > 0) {
      xv <- cont_vars[1]
      xseq <- seq(min(df_mod[[xv]], na.rm = TRUE), max(df_mod[[xv]], na.rm = TRUE), length.out = 200)
      other_vals <- lapply(setdiff(x_names, xv), function(v) {
        col <- df_mod[[v]]
        if (is.numeric(col)) median(col, na.rm = TRUE) else {
          lv <- if (is.factor(col)) levels(col)[1] else sort(unique(col))[1]
          factor(lv, levels = if (is.factor(col)) levels(col) else sort(unique(col)))
        }
      })
      names(other_vals) <- setdiff(x_names, xv)
      newdata <- data.frame(setNames(list(xseq), xv), other_vals, check.names = FALSE)
      pred_p <- predict(mod, newdata = newdata, type = "probs")
      if (is.vector(pred_p)) pred_p <- matrix(pred_p, ncol = 1, dimnames = list(NULL, cats[1]))

      df_plt <- data.frame(x = xseq)
      for (lv in colnames(pred_p)) df_plt[[lv]] <- pred_p[, lv]
      df_long <- stats::reshape(df_plt, varying = colnames(pred_p),
                                v.names = "prob", timevar = "categoria",
                                times = colnames(pred_p), direction = "long")

      p_graf <- ggplot2::ggplot(df_long, ggplot2::aes(x = x, y = prob, color = categoria)) +
        ggplot2::geom_line(linewidth = 1.2) +
        ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        ggplot2::labs(title = "Probabilidades Preditas por Categoria",
                      subtitle = sprintf("Variando '%s', demais fixos na mediana/moda", xv),
                      x = xv, y = "Probabilidade Predita", color = y_name,
                      caption = "estatR") +
        tema_g
      suppressMessages(suppressWarnings(print(p_graf)))
    } else {
      # Sem contínua: barplot das probs médias por nível
      df_bar <- as.data.frame(colMeans(as.matrix(probs_mat)))
      colnames(df_bar) <- "prob"
      df_bar$categoria <- rownames(df_bar)
      p_graf <- ggplot2::ggplot(df_bar, ggplot2::aes(x = categoria, y = prob, fill = categoria)) +
        ggplot2::geom_col(alpha = 0.85, width = 0.6) +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::labs(title = "Probabilidade Predita Média por Categoria",
                      x = y_name, y = "Probabilidade Média", caption = "estatR") +
        ggplot2::guides(fill = "none") +
        tema_g
      suppressMessages(suppressWarnings(print(p_graf)))
    }
  }
  invisible(mod)
}

# ─────────────────────────────────────────────────────────────────────────────
# ORDINAL (MASS::polr)
# ─────────────────────────────────────────────────────────────────────────────
.logistica_ordinal <- function(frm, dados, grafico, y_name, x_names, niveis) {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Instale o pacote 'MASS' para regressão ordinal: install.packages('MASS')")

  df_mod <- na.omit(dados[, c(y_name, x_names), drop = FALSE])
  df_mod[[y_name]] <- factor(df_mod[[y_name]], levels = niveis, ordered = TRUE)
  n_obs  <- nrow(df_mod)

  mod <- tryCatch(
    do.call(MASS::polr, list(formula = frm, data = df_mod, Hess = TRUE)),
    error = function(e) stop("Erro ao ajustar modelo ordinal: ", e$message)
  )

  sm     <- summary(mod)
  coefs  <- sm$coefficients          # (p + J-1) x 3 (coef, se, t)
  ci_mat <- suppressMessages(confint(mod))
  p_vals <- 2 * pnorm(-abs(coefs[, "t value"]))

  # ── Métricas ───────────────────────────────────────────────────────────────
  dev_res  <- deviance(mod)
  # Modelo nulo: só interceptos
  mod_nulo <- tryCatch(
    do.call(MASS::polr, list(formula = stats::as.formula(paste(y_name, "~ 1")), data = df_mod, Hess = FALSE)),
    error = function(e) NULL
  )
  dev_nula  <- if (!is.null(mod_nulo)) deviance(mod_nulo) else NA
  lrt_chi   <- if (!is.na(dev_nula)) dev_nula - dev_res else NA
  lrt_df    <- length(x_names)
  lrt_p     <- if (!is.na(lrt_chi)) 1 - pchisq(lrt_chi, lrt_df) else NA
  pseudo_r2 <- if (!is.na(dev_nula)) 1 - dev_res / dev_nula else NA
  aic_val   <- AIC(mod)
  bic_val   <- BIC(mod)
  K         <- length(coef(mod)) + length(mod$zeta) + 1
  aicc_val  <- if (n_obs - K - 1 > 0) aic_val + (2*K*(K+1))/(n_obs-K-1) else Inf

  pred_cls <- as.character(predict(mod))
  acc      <- mean(pred_cls == as.character(df_mod[[y_name]]))

  # ── Impressão ─────────────────────────────────────────────────────────────
  .print_titulo("REGRESSÃO LOGÍSTICA ORDINAL")
  cat(sprintf("  Variável Resposta:   %s  (%d níveis ordenados)\n", y_name, length(niveis)))
  cat(sprintf("  Ordem das Categorias: %s\n", paste(niveis, collapse = " < ")))
  cat(sprintf("  Observações Válidas: %d\n\n", n_obs))

  .print_topico("RESUMO GERAL (AJUSTE)")
  w_aj <- c(med = 26, val = 30)
  aj_sep <- .sep_log(w_aj["med"] + w_aj["val"])
  aj_row <- function(m, v) cat("  ", .pad_log(m, w_aj["med"], "left"), .pad_log(v, w_aj["val"], "right"), "\n", sep = "")
  cat(aj_sep, "\n")
  cat("  ", .pad_log("Métrica", w_aj["med"], "left"), .pad_log("Valor", w_aj["val"], "right"), "\n", sep = "")
  cat(aj_sep, "\n")
  if (!is.na(dev_nula)) aj_row("Deviance Nula",     sprintf("%.2f", dev_nula))
  aj_row("Deviance Residual", sprintf("%.2f", dev_res))
  if (!is.na(lrt_chi))
    aj_row("Teste LRT (\u03c7\u00b2)", sprintf("%.2f  (p %s)  %s", lrt_chi,
                                              if (lrt_p < 0.001) "< 0,001" else paste0("= ", formatC(lrt_p, format="f", digits=3, decimal.mark=",")),
                                              .asterisk_log(lrt_p)))
  if (!is.na(pseudo_r2)) aj_row("Pseudo-R\u00b2 (McFadden)", .fmt_pct(pseudo_r2, 1))
  aj_row("Acurácia Global",   .fmt_pct(acc, 1))
  aj_row("AIC",  sprintf("%.2f", aic_val))
  aj_row("AICc", sprintf("%.2f", aicc_val))
  aj_row("BIC",  sprintf("%.2f", bic_val))
  cat(aj_sep, "\n")
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n\n")

  # ── Coeficientes dos Preditores ────────────────────────────────────────────
  .print_topico("COEFICIENTES DOS PREDITORES (Odds Ratio Acumulados)")
  pred_idx <- seq_along(x_names)
  w_c <- c(var = max(14, max(nchar(x_names)) + 2),
            est = 12, ep = 12, t = 10, p = 16, or = 10, ic_inf = 12, ic_sup = 12)
  hdr <- paste0(
    .pad_log("Variável",   w_c["var"], "left"),
    .pad_log("Estimativa", w_c["est"], "center"),
    .pad_log("Erro Pad.",  w_c["ep"],  "center"),
    .pad_log("t",          w_c["t"],   "center"),
    .pad_log("p-valor",    w_c["p"],   "center"),
    .pad_log("OR",         w_c["or"],  "center"),
    .pad_log("IC 95% Inf", w_c["ic_inf"], "center"),
    .pad_log("IC 95% Sup", w_c["ic_sup"], "center")
  )
  w_tot <- nchar(hdr)
  cat(.sep_log(w_tot), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.sep_log(w_tot), "\n")
  for (j in pred_idx) {
    vn  <- x_names[j]
    est <- coefs[vn, "Value"]
    se  <- coefs[vn, "Std. Error"]
    tv  <- coefs[vn, "t value"]
    pv  <- p_vals[vn]
    or  <- exp(-est)   # polr usa escala negativa; OR acumulado = exp(-coef)
    ic_inf <- if (vn %in% rownames(ci_mat)) exp(-ci_mat[vn, 2]) else NA
    ic_sup <- if (vn %in% rownames(ci_mat)) exp(-ci_mat[vn, 1]) else NA
    p_txt <- paste(.formata_p_local(pv), .asterisk_log(pv))
    cat("  ",
        .pad_log(vn,               w_c["var"], "left"),
        .pad_log(sprintf("%.4f", est), w_c["est"], "center"),
        .pad_log(sprintf("%.4f", se),  w_c["ep"],  "center"),
        .pad_log(sprintf("%.2f", tv),  w_c["t"],   "center"),
        .pad_log(p_txt,            w_c["p"],   "center"),
        .pad_log(.fmt_or(or),      w_c["or"],  "center"),
        .pad_log(if (!is.na(ic_inf)) .fmt_or(ic_inf) else "-", w_c["ic_inf"], "center"),
        .pad_log(if (!is.na(ic_sup)) .fmt_or(ic_sup) else "-", w_c["ic_sup"], "center"),
        "\n", sep = "")
  }
  cat(.sep_log(w_tot), "\n")
  cat("  OR > 1: associado a categorias mais altas na escala ordinal.\n\n")

  # ── Interceptos (Limiares / Zeta) ─────────────────────────────────────────
  .print_topico("LIMIARES (Interceptos de Corte — \u03b6)")
  zetas <- mod$zeta
  wz <- c(lim = max(16, max(nchar(names(zetas))) + 2), val = 12, or = 12)
  zeta_sep <- .sep_log(wz["lim"] + wz["val"] + wz["or"])
  cat(zeta_sep, "\n")
  cat("  ", .pad_log("Limiar", wz["lim"], "left"),
      .pad_log("Estimativa", wz["val"], "center"),
      .pad_log("OR Acumulado", wz["or"], "center"),
      "\n", sep = "")
  cat(zeta_sep, "\n")
  for (nm in names(zetas)) {
    cat("  ",
        .pad_log(nm,                       wz["lim"], "left"),
        .pad_log(sprintf("%.4f", zetas[nm]),  wz["val"], "center"),
        .pad_log(.fmt_or(exp(zetas[nm])),   wz["or"],  "center"),
        "\n", sep = "")
  }
  cat(zeta_sep, "\n\n")

  # ── Interpretação ──────────────────────────────────────────────────────────
  .print_topico("INTERPRETAÇÃO")
  cat(sprintf("  O modelo ordinal assume proporcionalidade dos odds ao longo das\n  categorias de '%s' (%s).\n",
              y_name, paste(niveis, collapse = " < ")))
  if (!is.na(lrt_p)) {
    cat(sprintf("\n  O modelo %s (LRT: \u03c7\u00b2 = %.2f, p %s),\n  com Pseudo-R\u00b2 de McFadden de %s e acurácia global de %s.\n",
                if (lrt_p < 0.05) "é globalmente significativo" else "não é globalmente significativo",
                lrt_chi,
                if (lrt_p < 0.001) "< 0,001" else paste0("= ", formatC(lrt_p, format="f", digits=3, decimal.mark=",")),
                .fmt_pct(pseudo_r2, 1), .fmt_pct(acc, 1)))
  }
  for (j in pred_idx) {
    vn <- x_names[j]
    pv <- p_vals[vn]
    if (!is.na(pv) && pv < 0.10) {
      or_val  <- exp(-coefs[vn, "Value"])
      direcao <- if (or_val > 1) "aumenta" else "reduz"
      p_txt   <- if (pv < 0.001) "p < 0,001" else sprintf("p = %s", formatC(pv, format="f", digits=3, decimal.mark=","))
      cat(sprintf("\n  '%s' (%s): um aumento de uma unidade %s as chances acumuladas\n  de estar em uma categoria superior em %.1f%% (OR = %.3f).\n",
                  vn, p_txt, direcao, abs(or_val - 1)*100, or_val))
    }
  }
  cat("\n")
  .print_rodape()

  # ── Gráfico: Curvas de Probabilidade Acumulada ─────────────────────────────
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tema_g <- .tema_log()
    cont_vars <- x_names[sapply(df_mod[x_names], is.numeric)]

    if (length(cont_vars) > 0) {
      xv <- cont_vars[1]
      xseq <- seq(min(df_mod[[xv]], na.rm = TRUE), max(df_mod[[xv]], na.rm = TRUE), length.out = 200)
      other_vals <- lapply(setdiff(x_names, xv), function(v) {
        col <- df_mod[[v]]
        if (is.numeric(col)) median(col, na.rm = TRUE) else {
          lv <- if (is.factor(col)) levels(col)[1] else sort(unique(col))[1]
          factor(lv, levels = if (is.factor(col)) levels(col) else sort(unique(col)))
        }
      })
      names(other_vals) <- setdiff(x_names, xv)
      newdata <- data.frame(setNames(list(xseq), xv), other_vals, check.names = FALSE)
      newdata[[y_name]] <- factor(niveis[1], levels = niveis, ordered = TRUE)
      pred_p <- predict(mod, newdata = newdata, type = "probs")

      df_plt <- data.frame(x = xseq)
      for (lv in colnames(pred_p)) df_plt[[lv]] <- pred_p[, lv]
      df_long <- stats::reshape(df_plt, varying = colnames(pred_p),
                                v.names = "prob", timevar = "categoria",
                                times = colnames(pred_p), direction = "long")
      df_long$categoria <- factor(df_long$categoria, levels = niveis)

      p_graf <- ggplot2::ggplot(df_long, ggplot2::aes(x = x, y = prob, color = categoria)) +
        ggplot2::geom_line(linewidth = 1.2) +
        ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        ggplot2::labs(title = "Probabilidades Preditas por Nível Ordinal",
                      subtitle = sprintf("Variando '%s', demais fixos na mediana/moda", xv),
                      x = xv, y = "Probabilidade Predita", color = y_name,
                      caption = "estatR") +
        tema_g
      suppressMessages(suppressWarnings(print(p_graf)))
    }
  }
  invisible(mod)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Métricas do Modelo Logístico
#' @description Exibe as métricas de ajuste e performance de um modelo logístico
#' ajustado por \code{regressao_logistica()}.
#' @param modelo Objeto \code{glm} (binomial) ou \code{multinom}/\code{polr} retornado por \code{regressao_logistica()}.
#' @param corte_prob Ponto de corte para classificação binária. Padrão 0.5.
#' @return Retorna invisivelmente um data frame com as métricas.
#' @export
metricas_logistica <- function(modelo, corte_prob = 0.5) {
  if (inherits(modelo, c("multinom", "polr"))) stop("Esta função de métricas detalhadas (Curva ROC, Sensibilidade, etc.) é exclusiva para regressão logística binomial.\nPara modelos multinomiais e ordinais, as métricas globais (AIC, Acurácia, Pseudo-R2) já são exibidas no resumo principal.")
  if (!inherits(modelo, "glm")) stop("'modelo' deve ser um objeto glm retornado por regressao_logistica().")

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
  acc     <- (vp + vn) / n_obs
  sensib  <- if ((vp + fn) > 0) vp / (vp + fn) else NA
  especif <- if ((vn + fp) > 0) vn / (vn + fp) else NA
  ppv     <- if ((vp + fp) > 0) vp / (vp + fp) else NA
  npv     <- if ((vn + fn) > 0) vn / (vn + fn) else NA

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
  cat("  ", .pad_log("Sensibilidade",  lw, "left"), .pad_log(if (!is.na(sensib))  .fmt_pct(sensib, 2)  else "-", vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("Especificidade", lw, "left"), .pad_log(if (!is.na(especif)) .fmt_pct(especif, 2) else "-", vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("VPP (Precisão)", lw, "left"), .pad_log(if (!is.na(ppv))     .fmt_pct(ppv, 2)     else "-", vw, "right"), "\n", sep = "")
  cat("  ", .pad_log("VPN",            lw, "left"), .pad_log(if (!is.na(npv))     .fmt_pct(npv, 2)     else "-", vw, "right"), "\n", sep = "")
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
#' @description Avalia múltiplas combinações de variáveis para regressão logística,
#' ordenando pelo AIC. Suporta modelos binomiais, multinomiais e ordinais.
#' @param formula Fórmula com todos os candidatos (ex: y ~ x1 + x2 + x3).
#' @param dados Data frame contendo os dados.
#' @param tipo Tipo para Y policotômico: \code{"nominal"} (padrão) ou \code{"ordinal"}.
#' @param top Número de modelos a exibir. Padrão 5.
#' @return Retorna invisivelmente um data frame com o ranking.
#' @export
selecao_modelos_logistica <- function(formula, dados, tipo = c("nominal", "ordinal"), top = 5) {
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  tipo <- match.arg(tipo)

  termos   <- attr(terms(formula, data = dados), "term.labels")
  y_name   <- as.character(formula[[2]])
  cols_req <- c(y_name, termos)
  df_clean <- na.omit(dados[, cols_req, drop = FALSE])
  n        <- nrow(df_clean)
  if (n < 5) stop("Poucas observações válidas (mínimo 5).")

  Y_raw  <- df_clean[[y_name]]
  niveis <- unique(as.character(Y_raw))
  n_niv  <- length(niveis)
  is_bin <- n_niv == 2

  p_total <- length(termos)
  usou_heuristica <- FALSE

  if (p_total > 10) {
    usou_heuristica <- TRUE
    Y_num <- as.integer(factor(df_clean[[y_name]])) - 1
    num_t <- termos[sapply(df_clean[termos], is.numeric)]
    cors  <- sapply(num_t, function(v) abs(cor(Y_num, df_clean[[v]], use = "complete.obs")))
    top10 <- names(sort(cors, decreasing = TRUE))[1:min(10, length(cors))]
    termos <- unique(c(top10, setdiff(termos, num_t)))[1:min(10, p_total)]
  }

  mods_list <- unlist(lapply(seq_along(termos), function(k) combn(termos, k, simplify = FALSE)), recursive = FALSE)
  total_av  <- length(mods_list)

  .fit_aic <- function(preds) {
    f_tmp <- stats::as.formula(paste(y_name, "~", paste(preds, collapse = " + ")))
    mod <- tryCatch({
      if (is_bin) {
        glm(f_tmp, data = df_clean, family = binomial())
      } else if (tipo == "ordinal") {
        df_ord <- df_clean
        df_ord[[y_name]] <- factor(df_ord[[y_name]], levels = sort(niveis), ordered = TRUE)
        MASS::polr(f_tmp, data = df_ord, Hess = FALSE)
      } else {
        nnet::multinom(f_tmp, data = df_clean, trace = FALSE)
      }
    }, error = function(e) NULL)
    if (is.null(mod)) return(NULL)
    k_p    <- length(coef(mod)) + 1
    aic_v  <- AIC(mod)
    aicc_v <- if (n - k_p - 1 > 0) aic_v + (2*k_p*(k_p+1))/(n-k_p-1) else Inf
    list(modelo = paste(preds, collapse = " + "), k = length(preds), AIC = aic_v, AICc = aicc_v)
  }

  res    <- Filter(Negate(is.null), lapply(mods_list, function(p) .fit_aic(p)))
  df_res <- do.call(rbind, lapply(res, as.data.frame))
  df_res <- df_res[order(df_res$AIC), ]
  df_top <- if (nrow(df_res) > top) df_res[1:top, ] else df_res

  tipo_txt <- if (is_bin) "BINOMIAL" else if (tipo == "ordinal") "ORDINAL" else "MULTINOMIAL"
  .print_titulo(sprintf("SELEÇÃO DE MODELOS — REGRESSÃO LOGÍSTICA %s", tipo_txt))
  cat(sprintf("  Variável Resposta: %s\n", y_name))
  cat(sprintf("  Método: %s (%d modelos avaliados)\n\n",
              if (usou_heuristica) "Heurístico (TOP 10)" else "Exaustivo", total_av))

  .print_topico(sprintf("TOP %d MODELOS (Ordenados por AIC)", nrow(df_top)))

  df_fmt <- data.frame(
    Ranking = paste0(1:nrow(df_top), "\u00ba"),
    Modelo  = df_top$modelo,
    k       = as.character(df_top$k),
    AIC     = sapply(df_top$AIC,  .fmt_num, decimais = 1),
    AICc    = sapply(df_top$AICc, .fmt_num, decimais = 1),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_fmt, align = c("left","left","center","center","center"))
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
#' @description Diagnóstico de resíduos de um modelo logístico binomial
#' (Hosmer-Lemeshow, dispersão, influência, painel 2x2).
#' @param modelo Objeto \code{glm} retornado por \code{regressao_logistica()}.
#' @param grafico Se TRUE, exibe o painel de diagnóstico. Padrão TRUE.
#' @return Retorna invisivelmente uma lista com os diagnósticos.
#' @export
analise_residual_logistica <- function(modelo, grafico = TRUE) {
  if (inherits(modelo, c("multinom", "polr"))) stop("A análise de resíduos avançada (Hosmer-Lemeshow, Pearson) atualmente está implementada apenas para a regressão logística binomial no pacote estatR.")
  if (!inherits(modelo, "glm")) stop("'modelo' deve ser um objeto glm retornado por regressao_logistica().")

  Y_num   <- as.integer(modelo$y)
  probs   <- modelo$fitted.values
  n       <- length(Y_num)

  res_pearson <- residuals(modelo, type = "pearson")
  res_devian  <- residuals(modelo, type = "deviance")
  alavancagem <- hatvalues(modelo)
  dist_cook   <- cooks.distance(modelo)
  hl          <- .hosmer_lemeshow_estatR(Y_num, probs, g = 10)
  chi_pearson <- sum(res_pearson^2)
  df_pearson  <- n - length(coef(modelo))
  disp_pearson <- chi_pearson / df_pearson

  formata_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }

  .print_titulo("ANÁLISE DE RESÍDUOS — REGRESSÃO LOGÍSTICA")

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

  .print_topico("INTERPRETAÇÃO DO DIAGNÓSTICO")
  if (!is.na(hl$p.value)) {
    if (hl$p.value > 0.05) {
      cat(sprintf("  [\u2713] : Hosmer-Lemeshow não rejeitou H\u2080 (p = %s).\n  As probabilidades preditas estão bem calibradas.\n\n", formata_p(hl$p.value)))
    } else {
      cat(sprintf("  [!] : Hosmer-Lemeshow rejeitou H\u2080 (p = %s).\n  Possível problema de especificação do modelo.\n\n", formata_p(hl$p.value)))
    }
  }
  if (disp_pearson > 1.5) {
    cat(sprintf("  [!] : Dispersão de Pearson = %.3f (> 1,5). Pode indicar\n  sobredispersão ou pontos influentes.\n\n", disp_pearson))
  } else {
    cat(sprintf("  [\u2713] : Dispersão de Pearson = %.3f. Sem sinais de sobredispersão relevante.\n\n", disp_pearson))
  }
  if (length(infl_idx) > 0) {
    cat(sprintf("  [!] : %d observação(ões) com alta influência detectada(s).\n  Considere investigar ou remover essas observações.\n\n", length(infl_idx)))
  }
  cat(sep_t, "\n\n")
  .print_rodape()

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tema_painel <- .tema_log()
    df_diag <- data.frame(idx = seq_len(n), ajustados = probs,
                          res_pearson = res_pearson, res_devian = res_devian,
                          cook = dist_cook, alavancagem = alavancagem)
    p1 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = ajustados, y = res_pearson)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "gray70") +
      ggplot2::geom_point(color = "#555555", alpha = 0.8) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429", linewidth = 0.7, formula = y ~ x) +
      ggplot2::labs(title = "Resíduos de Pearson vs Predito",
                    x = "Probabilidade Predita", y = "Resíduo de Pearson") + tema_painel
    p2 <- ggplot2::ggplot(df_diag, ggplot2::aes(sample = res_devian)) +
      ggplot2::stat_qq(color = "#555555", alpha = 0.8) +
      ggplot2::stat_qq_line(color = "#D90429", linewidth = 0.7) +
      ggplot2::labs(title = "QQ-Plot (Resíduos de Desvio)",
                    x = "Quantis Teóricos", y = "Quantis Amostrais") + tema_painel
    p3 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = idx, y = cook)) +
      ggplot2::geom_hline(yintercept = 4/n, linetype = "dashed", color = "#D90429") +
      ggplot2::geom_segment(ggplot2::aes(xend = idx, yend = 0), color = "#555555", alpha = 0.7) +
      ggplot2::geom_point(color = "#555555", alpha = 0.8) +
      ggplot2::labs(title = "Distância de Cook", x = "Observação", y = "Cook") + tema_painel
    p4 <- ggplot2::ggplot(df_diag, ggplot2::aes(x = alavancagem, y = res_pearson)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      ggplot2::geom_point(color = "#555555", alpha = 0.8) +
      ggplot2::labs(title = "Alavancagem vs Resíduo de Pearson",
                    x = "Alavancagem (Leverage)", y = "Resíduo de Pearson",
                    caption = "estatR") + tema_painel
    if (requireNamespace("patchwork", quietly = TRUE)) {
      suppressMessages(suppressWarnings(print((p1 | p2) / (p3 | p4))))
    } else {
      suppressMessages(suppressWarnings({ print(p1); print(p2); print(p3); print(p4) }))
    }
  }
  invisible(list(hosmer_lemeshow = hl, dispersao_pearson = disp_pearson, pontos_influentes = infl_idx))
}
