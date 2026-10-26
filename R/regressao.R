#' @title Regressão Linear Simples e Múltipla
#' @description Ajusta um modelo de regressão linear e retorna uma saída formatada com
#'   ANOVA, coeficientes (com IC 95% e VIF), ajuste do modelo e interpretação automática.
#' @param formula Fórmula do modelo (ex: \code{y ~ x} ou \code{y ~ x1 + x2})
#' @param dados Data frame com os dados
#' @param grafico Lógico. Se TRUE, gera gráfico de dispersão (apenas regressão simples).
#' @param ... Outros argumentos (reservado para uso futuro).
#' @export
regressao_linear <- function(formula, dados, grafico = TRUE, ...) {

  # ── 1. Ajuste ─────────────────────────────────────────────────────────────
  modelo  <- lm(formula, data = dados)
  resumo  <- summary(modelo)
  coefs   <- resumo$coefficients
  r2      <- resumo$r.squared
  r2_adj  <- resumo$adj.r.squared
  rmse    <- sqrt(mean(modelo$residuals^2))
  f_stat  <- resumo$fstatistic

  p_global <- if (!is.null(f_stat))
    pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  else NA

  vars   <- all.vars(formula)
  y_name <- vars[1]
  x_names <- vars[-1]

  # ── 2. IC 95% ─────────────────────────────────────────────────────────────
  ic_vals <- suppressMessages(confint(modelo, level = 0.95))

  # ── 3. VIF (somente regressão múltipla, base R) ───────────────────────────
  tem_vif  <- FALSE
  vif_vals <- NULL
  if (length(x_names) > 1) {
    X <- model.matrix(modelo)[, -1, drop = FALSE]
    vif_vals <- numeric(ncol(X))
    for (i in seq_len(ncol(X))) {
      mod_vif    <- lm(X[, i] ~ X[, -i, drop = FALSE])
      r2_vif     <- summary(mod_vif)$r.squared
      vif_vals[i] <- if (r2_vif < 1) 1 / (1 - r2_vif) else Inf
    }
    tem_vif <- TRUE
  }

  # ── 4. Helpers ────────────────────────────────────────────────────────────
  fmt_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("< 0.001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }

  asterisk <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    if (p < 0.10) return("*")
    return("")
  }

  pad <- function(s, w, align = "left") {
    s  <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right")  return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }

  sep_line <- function(w) paste0("  ", strrep("\u2500", w))

  # ── 5. Cabeçalho ──────────────────────────────────────────────────────────
  tipo_reg_txt <- if (length(x_names) > 1) "MÚLTIPLA" else "SIMPLES"
  .print_titulo(paste("RESULTADOS DA REGRESSÃO", tipo_reg_txt))
  cat(sprintf("  Equa\u00e7\u00e3o do modelo: %s\n", deparse(formula)))
  cat(sprintf("  Observa\u00e7\u00f5es: %d\n\n", nrow(modelo$model)))

  # ── 6. ANOVA ──────────────────────────────────────────────────────────────
  .print_topico("QUADRO ANOVA")

  tab_aov <- anova(modelo)
  idx_res <- which(rownames(tab_aov) == "Residuals")
  df_res  <- tab_aov[idx_res, "Df"]
  sq_res  <- tab_aov[idx_res, "Sum Sq"]
  qm_res  <- sq_res / df_res

  df_mod  <- sum(tab_aov[-idx_res, "Df"])
  sq_mod  <- sum(tab_aov[-idx_res, "Sum Sq"])
  qm_mod  <- sq_mod / df_mod
  f_mod   <- if (!is.na(p_global)) f_stat[1] else NA
  df_tot  <- df_mod + df_res
  sq_tot  <- sq_mod + sq_res

  w_a <- c(fv = 22, gl = 11, sq = 12, qm = 13, f = 9, p = 14)
  hdr_a <- paste0(
    pad("Fonte de Varia\u00e7\u00e3o", w_a["fv"], "left"),
    pad("Graus Lib.", w_a["gl"], "center"),
    pad("Soma Quad.", w_a["sq"], "center"),
    pad("Quad. M\u00e9dio", w_a["qm"], "center"),
    pad("F", w_a["f"], "center"),
    pad("p-valor", w_a["p"], "center")
  )
  w_total <- nchar(hdr_a)

  cat(sep_line(w_total), "\n")
  cat("  ", hdr_a, "\n", sep = "")
  cat(sep_line(w_total), "\n")

  fmt_aov_row <- function(fv, gl, sq, qm, f_val, p_val) {
    p_txt <- if (!is.na(p_val)) paste(fmt_p(p_val), asterisk(p_val)) else ""
    f_txt <- if (!is.na(f_val)) sprintf("%.2f", f_val) else ""
    cat("  ",
        pad(fv,  w_a["fv"],  "left"),
        pad(gl,  w_a["gl"],  "center"),
        pad(sprintf("%.1f", sq), w_a["sq"], "center"),
        pad(sprintf("%.1f", qm), w_a["qm"], "center"),
        pad(f_txt, w_a["f"], "center"),
        pad(p_txt, w_a["p"], "center"),
        "\n", sep = "")
  }

  fmt_aov_row("Modelo (Regress\u00e3o)", df_mod, sq_mod, qm_mod, f_mod, p_global)
  cat("  ",
      pad("Res\u00edduos (Erro)", w_a["fv"], "left"),
      pad(df_res, w_a["gl"], "center"),
      pad(sprintf("%.1f", sq_res), w_a["sq"], "center"),
      pad(sprintf("%.1f", qm_res), w_a["qm"], "center"),
      pad("", w_a["f"], "center"),
      pad("", w_a["p"], "center"),
      "\n", sep = "")
  cat("  ",
      pad("Total", w_a["fv"], "left"),
      pad(df_tot, w_a["gl"], "center"),
      pad(sprintf("%.1f", sq_tot), w_a["sq"], "center"),
      pad("", w_a["qm"], "center"),
      pad("", w_a["f"], "center"),
      pad("", w_a["p"], "center"),
      "\n", sep = "")

  cat(sep_line(w_total), "\n")
  cat("  Notas:\n")
  cat("  Signific\u00e2ncia: *** p < 0.01   ** p < 0.05   * p < 0.10\n\n")

  # ── 7. Coeficientes ───────────────────────────────────────────────────────
  .print_topico("COEFICIENTES DO MODELO")

  rn <- rownames(coefs)
  rn[rn == "(Intercept)"] <- "(Intercepto)"

  w_c <- c(var = 16, est = 13, err = 12, ic = 18, p = 16, vif = if (tem_vif) 8 else 0)

  hdr_c <- paste0(
    pad("Vari\u00e1vel",    w_c["var"], "left"),
    pad("Estimativa", w_c["est"], "center"),
    pad("Erro Pad.",  w_c["err"], "center"),
    pad("IC (95%)",   w_c["ic"],  "center"),
    pad("p-valor",    w_c["p"],   "center"),
    if (tem_vif) pad("VIF", w_c["vif"], "center") else ""
  )
  w_coef <- nchar(hdr_c)

  cat(sep_line(w_coef), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(sep_line(w_coef), "\n")

  for (i in seq_len(nrow(coefs))) {
    ic_str  <- sprintf("[%.2f; %.2f]", ic_vals[i, 1], ic_vals[i, 2])
    p_ast   <- paste(fmt_p(coefs[i, 4]), asterisk(coefs[i, 4]))
    vif_str <- ""
    if (tem_vif) {
      vif_str <- if (i == 1) "-" else sprintf("%.2f", vif_vals[i - 1])
    }
    cat("  ",
        pad(rn[i],                           w_c["var"], "left"),
        pad(sprintf("%.3f", coefs[i, 1]),    w_c["est"], "center"),
        pad(sprintf("%.3f", coefs[i, 2]),    w_c["err"], "center"),
        pad(ic_str,                          w_c["ic"],  "center"),
        pad(p_ast,                           w_c["p"],   "center"),
        if (tem_vif) pad(vif_str, w_c["vif"], "center") else "",
        "\n", sep = "")
  }

  cat(sep_line(w_coef), "\n")
  cat("  Notas:\n")
  cat("  Signific\u00e2ncia: *** p < 0.01   ** p < 0.05   * p < 0.10\n")
  if (tem_vif)
    cat("  VIF: Valores > 5 indicam multicolinearidade moderada; > 10 grave.\n")
  cat("\n")

  # ── 8. Ajuste do modelo ───────────────────────────────────────────────────
  .print_topico("AJUSTE DO MODELO")

  w_aj <- c(med = 22, val = 30)
  hdr_aj <- paste0(pad("Medidas de Ajuste", w_aj["med"], "left"),
                   pad("Valor",             w_aj["val"], "left"))
  w_ajuste <- nchar(hdr_aj)

  cat(sep_line(w_ajuste), "\n")
  cat("  ", hdr_aj, "\n", sep = "")
  cat(sep_line(w_ajuste), "\n")

  aj_row <- function(med, val) {
    cat("  ", pad(med, w_aj["med"], "left"), pad(val, w_aj["val"], "left"), "\n", sep = "")
  }

  f_txt <- if (!is.na(p_global))
    sprintf("%.2f (p %s %s)", f_stat[1],
            if (p_global < 0.001) "< 0.001" else paste0("= ", formatC(p_global, format="f", digits=3, decimal.mark=",")),
            asterisk(p_global))
  else "NA"

  aj_row("R\u00b2 (%):", formatC(r2     * 100, format="f", digits=2, decimal.mark=","))
  aj_row("R\u00b2 ajustado (%):", formatC(r2_adj * 100, format="f", digits=2, decimal.mark=","))
  aj_row("RMSE:",           sprintf("%.2f", rmse))
  aj_row("F-estat\u00edstica:",   f_txt)
  cat(sep_line(w_ajuste), "\n\n")

  # ── 9. Interpretação ──────────────────────────────────────────────────────
  .print_topico("INTERPRETA\u00c7\u00c3O")
  cat(sprintf("  O modelo de regress\u00e3o linear foi ajustado para explicar %s a partir de %s.\n",
              y_name, paste(x_names, collapse = " e ")))

  valido_txt <- if (!is.na(p_global) && p_global < 0.05) "\u00e9 globalmente v\u00e1lido" else "n\u00e3o \u00e9 globalmente v\u00e1lido"
  p_global_txt <- if (!is.na(p_global))
    if (p_global < 0.001) "< 0.001" else paste0("= ", formatC(p_global, format="f", digits=3, decimal.mark=","))
  else "NA"

  cat(sprintf("\n  O modelo %s (p %s) e explica aproximadamente %s\n  da variabilidade observada em %s.\n",
              valido_txt, p_global_txt,
              paste0(formatC(r2 * 100, format="f", digits=1, decimal.mark=","), "%"),
              y_name))

  for (i in 2:nrow(coefs)) {
    p_val <- coefs[i, 4]
    if (p_val < 0.10) {
      direcao <- if (coefs[i, 1] > 0) "positiva" else "negativa"
      mantendo <- if (length(x_names) > 1) ", mantendo os demais constantes" else ""
      p_texto <- if (p_val < 0.001) "p < 0.001"
                 else sprintf("p = %s", formatC(p_val, format="f", digits=3, decimal.mark=","))
      cat(sprintf("\n  O preditor '%s' apresentou associa\u00e7\u00e3o %s estatisticamente significativa\n  (%s) com %s%s.\n",
                  rn[i], direcao, p_texto, y_name, mantendo))
    }
  }

  for (i in 2:nrow(coefs)) {
    beta <- coefs[i, 1]
    nome_pred <- rn[i]

    ceteris <- if (length(x_names) > 1) {
      outros <- rn[2:nrow(coefs)]
      outros <- outros[outros != nome_pred]
      sprintf("Mantendo %s constante(s), um", paste(outros, collapse = " e "))
    } else "Um"

    verbo <- if (beta > 0) "um aumento estimado de" else "uma redu\u00e7\u00e3o estimada de"
    cat(sprintf("\n  %s aumento de uma unidade em %s est\u00e1 associado\n  a %s %.2f unidades em %s.\n",
                ceteris, nome_pred, verbo, abs(beta), y_name))
  }

  .print_rodape()

  # ── 10. Gráfico (regressão simples) ───────────────────────────────────────
  if (grafico && length(x_names) == 1) {
    tryCatch({
      intercepto <- coefs[1, 1]
      inclinacao <- coefs[2, 1]
      sinal      <- if (inclinacao >= 0) "+" else "-"
      equacao_texto <- sprintf("%s = %.2f %s %.2f × %s   |   R² = %.1f%%",
                               y_name, intercepto, sinal, abs(inclinacao), x_names[1], r2 * 100)
      df_p <- dados[, c(x_names[1], y_name)]
      names(df_p) <- c("x", "y")
      df_p <- df_p[!is.na(df_p$x) & !is.na(df_p$y), ]

      # Posicionamento inteligente para evitar sobreposição com pontos:
      x_rng <- range(df_p$x, na.rm = TRUE)
      y_rng <- range(df_p$y, na.rm = TRUE)
      if (inclinacao >= 0) {
        ann_x <- x_rng[1] + 0.05 * diff(x_rng)
        ann_hjust <- 0
      } else {
        ann_x <- x_rng[2] - 0.05 * diff(x_rng)
        ann_hjust <- 1
      }
      ann_y <- y_rng[2] - 0.05 * diff(y_rng)

      old_w <- getOption("warn"); options(warn = -1)
      p_plot <- ggplot2::ggplot(df_p, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(color = "#3470a3", size = 2.8, alpha = 0.8) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#D90429", fill = "#85b4d4", alpha = 0.25, linewidth = 1) +
        ggplot2::annotate("text", x = ann_x, y = ann_y, label = equacao_texto,
                          hjust = ann_hjust, vjust = 1, size = 3.8,
                          color = "#222222", fontface = "italic") +
        ggplot2::labs(
          title    = sprintf("Regressão Linear: %s vs %s", y_name, x_names[1]),
          subtitle = sprintf("Modelo estimado com %d observações", nrow(df_p)),
          x        = x_names[1],
          y        = y_name,
          caption  = "estatR"
        ) +
        meu_tema(grade = "dupla")
      suppressMessages(print(p_plot))
      options(warn = old_w)
    }, error = function(e) {
      message("\n[Aviso] Não foi possível gerar o gráfico: ", e$message)
    })
  }

  invisible(modelo)
}


#' @title Predição com o Modelo de Regressão
#' @description Gera previsões a partir de um modelo linear ajustado, exibindo a estimativa
#'   pontual e, opcionalmente, o intervalo de confiança ou de predição.
#'   \itemize{
#'     \item \strong{Intervalo de Confiança:} margem para a \emph{média} da variável resposta
#'       em novas observações (menor incerteza).
#'     \item \strong{Intervalo de Predição:} margem para um valor \emph{individual} futuro
#'       (maior incerteza, sempre mais largo que o de confiança).
#'   }
#' @param modelo Objeto \code{lm} retornado por \code{regressao_linear()}.
#' @param novos_dados Data frame com as variáveis preditoras para as quais se deseja prever.
#' @param intervalo Tipo de intervalo: \code{"confianca"} (padrão), \code{"predicao"} ou \code{"nenhum"}.
#' @param nivel Nível de confiança (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @param grafico Lógico. Se TRUE e o modelo for simples (1 preditor), gera gráfico com
#'   a reta, a faixa do intervalo e os pontos previstos destacados em vermelho.
#' @return Retorna invisivelmente um data frame com as previsões e os limites do intervalo.
#' @export
predicao <- function(modelo, novos_dados,
                     intervalo = c("confianca", "predicao", "nenhum"),
                     nivel = 0.95, decimais = 2, grafico = TRUE) {

  intervalo <- match.arg(intervalo)
  if (!inherits(modelo, "lm"))      stop("'modelo' deve ser da classe 'lm'.")
  if (!is.data.frame(novos_dados))  stop("'novos_dados' deve ser um data.frame.")

  y_name  <- as.character(formula(modelo)[[2]])
  x_names <- all.vars(formula(modelo))[-1]

  # ── Calcular previsões ────────────────────────────────────────────────────
  if (intervalo == "nenhum") {
    preds  <- predict(modelo, newdata = novos_dados)
    df_res <- data.frame(Previsao = preds)
  } else {
    tipo   <- if (intervalo == "confianca") "confidence" else "prediction"
    preds  <- predict(modelo, newdata = novos_dados, interval = tipo, level = nivel)
    df_res <- as.data.frame(preds)
    colnames(df_res) <- c("Previsao", "LI", "LS")
  }

  df_x     <- novos_dados[, x_names, drop = FALSE]
  df_final <- cbind(df_x, df_res)

  # ── Helpers ───────────────────────────────────────────────────────────────
  pad <- function(s, w, align = "left") {
    s  <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right")  return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }

  fmt <- function(x) sprintf(paste0("%.", decimais, "f"), x)

  # ── Cabeçalho ─────────────────────────────────────────────────────────────
  .print_titulo(sprintf("PREVIS\u00d5ES DO MODELO (Vari\u00e1vel: %s)", y_name))

  nome_int <- switch(intervalo,
    confianca = "Confian\u00e7a",
    predicao  = "Predi\u00e7\u00e3o",
    nenhum    = NULL
  )
  if (!is.null(nome_int))
    cat(sprintf("  Tipo de Margem: Intervalo de %s (%.0f%%)\n\n", nome_int, nivel * 100))
  else
    cat("  Margens de erro (intervalos) n\u00e3o solicitadas.\n\n")

  # ── Larguras das colunas ──────────────────────────────────────────────────
  w_x <- pmax(
    sapply(x_names, nchar),
    sapply(df_x, function(col) max(nchar(as.character(col))))
  ) + 3
  w_x <- pmax(w_x, 9)
  w_p <- 14; w_l <- 22

  hdr_x <- paste(mapply(pad, x_names, w_x, "left"), collapse = "")

  hdr <- if (intervalo == "nenhum") {
    paste0(hdr_x, pad("Previs\u00e3o (Y)", w_p, "center"))
  } else {
    paste0(hdr_x,
           pad("Previs\u00e3o (Y)",       w_p, "center"),
           pad("Lim. Inferior (LI)", w_l, "center"),
           pad("Lim. Superior (LS)", w_l, "center"))
  }

  sep_hdr <- paste0("  ", strrep("\u2500", nchar(hdr)))
  cat(sep_hdr, "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(sep_hdr, "\n")

  # ── Linhas de dados ───────────────────────────────────────────────────────
  for (i in seq_len(nrow(df_final))) {
    linha_x <- paste(mapply(function(nm, w) pad(df_final[i, nm], w, "left"),
                            x_names, w_x), collapse = "")
    if (intervalo == "nenhum") {
      cat("  ", linha_x, pad(fmt(df_final[i, "Previsao"]), w_p, "center"), "\n", sep = "")
    } else {
      cat("  ", linha_x,
          pad(fmt(df_final[i, "Previsao"]), w_p, "center"),
          pad(fmt(df_final[i, "LI"]),       w_l, "center"),
          pad(fmt(df_final[i, "LS"]),       w_l, "center"),
          "\n", sep = "")
    }
  }
  cat(sep_hdr, "\n\n")
  .print_rodape()

  # ── Gráfico (regressão simples + intervalo) ───────────────────────────────
  if (grafico && length(x_names) == 1 && intervalo != "nenhum") {
    tryCatch({
      dados_orig <- stats::model.frame(modelo)
      old_w <- getOption("warn"); options(warn = -1)
      p_plot <- ggplot2::ggplot(dados_orig, ggplot2::aes(x = .data[[x_names[1]]], y = .data[[y_name]])) +
        ggplot2::geom_point(color = "#3470a3", size = 2.5, alpha = 0.6) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#2B5C8F", fill = "#85b4d4", alpha = 0.25, linewidth = 0.9) +
        ggplot2::geom_errorbar(data = df_final, ggplot2::aes(x = .data[[x_names[1]]], ymin = LI, ymax = LS),
                              color = "#D90429", width = 0.15, linewidth = 1, inherit.aes = FALSE) +
        ggplot2::geom_point(data = df_final, ggplot2::aes(x = .data[[x_names[1]]], y = Previsao),
                            color = "#D90429", size = 4.5, shape = 19, inherit.aes = FALSE) +
        ggplot2::labs(
          title    = sprintf("Predição: %s vs %s", y_name, x_names[1]),
          subtitle = sprintf("Pontos em vermelho: valores preditos com intervalo de %s (%.0f%%)", nome_int, nivel * 100),
          x        = x_names[1],
          y        = y_name,
          caption  = "estatR"
        ) +
        meu_tema(grade = "dupla")
      suppressMessages(print(p_plot))
      options(warn = old_w)
    }, error = function(e) {
      message("\n[Aviso] Não foi possível gerar o gráfico: ", e$message)
    })
  }

  invisible(df_final)
}


#' @title Métricas do Modelo
#' @description Exibe as principais métricas de avaliação e erro para um modelo de regressão.
#' @param modelo Objeto do tipo lm
#' @export
metricas <- function(modelo, corte_prob = 0.5) {
  if (inherits(modelo, "zeroinfl") && !is.null(modelo$dist)) {
    if (modelo$dist == "poisson" && exists(".metricas_poisson_iz", mode = "function")) {
      return(.metricas_poisson_iz(modelo))
    }
    if (modelo$dist == "negbin" && exists(".metricas_negbin_iz", mode = "function")) {
      return(.metricas_negbin_iz(modelo))
    }
  }
  if (inherits(modelo, c("estatR_cox", "coxph"))) return(metricas_cox(modelo))
  if (inherits(modelo, "betareg"))                return(metricas_beta(modelo))
  if (inherits(modelo, "negbin"))                 return(metricas_binomial_negativa(modelo))
  if (inherits(modelo, c("forecast", "previsao_serie", "Arima", "ets", "HoltWinters", "stlm"))) {
    return(metricas_serie(modelo))
  }
  if (inherits(modelo, "estatR_aprend_maq")) {
    if (exists("metricas_aprend_maq", mode = "function")) return(metricas_aprend_maq(modelo))
  }
  if (inherits(modelo, "glm") && !is.null(modelo$family)) {
    if (modelo$family$family == "binomial") return(metricas_logistica(modelo, corte_prob = corte_prob))
    if (modelo$family$family == "poisson")  return(metricas_poisson(modelo))
    if (modelo$family$family == "Gamma")    return(metricas_gama(modelo))
  }
  resumo <- summary(modelo)
  r2     <- resumo$r.squared
  r2_adj <- resumo$adj.r.squared

  y_obs  <- stats::model.response(stats::model.frame(modelo))
  y_hat  <- stats::fitted(modelo)
  res    <- stats::residuals(modelo)
  n      <- length(res)
  k      <- length(stats::coef(modelo))

  rmse     <- sqrt(mean(res^2))
  mae      <- mean(abs(res))
  mape     <- if (!is.null(y_obs) && all(y_obs != 0)) mean(abs(res / y_obs)) * 100 else NA
  smape    <- if (!is.null(y_obs)) mean(2 * abs(res) / (abs(y_obs) + abs(y_hat))) * 100 else NA
  aic_val  <- stats::AIC(modelo)
  bic_val  <- stats::BIC(modelo)
  aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)

  lw    <- 22
  vw    <- 20
  sep_m <- paste0("  ", strrep("─", lw + vw + 2), "\n")

  f4 <- function(x) formatC(x, format = "f", digits = 4, decimal.mark = ",")
  f2 <- function(x) formatC(x, format = "f", digits = 2, decimal.mark = ",")

  .print_titulo("MÉTRICAS E CRITÉRIOS DE SELEÇÃO")

  cat(sprintf("  %-*s %*s\n", lw, "Métrica", vw, "Valor"))
  cat(sep_m)
  cat(sprintf("  %-*s %*s\n", lw, "R² (%)",            vw, f2(r2 * 100)))
  cat(sprintf("  %-*s %*s\n", lw, "R² ajustado (%)",   vw, f2(r2_adj * 100)))
  cat(sprintf("  %-*s %*s\n", lw, "RMSE",          vw, f4(rmse)))
  cat(sprintf("  %-*s %*s\n", lw, "MAE",           vw, f4(mae)))
  if (!is.na(mape)) {
    cat(sprintf("  %-*s %*s\n", lw, "MAPE (%)",         vw, f2(mape)))
  }
  if (!is.na(smape)) {
    cat(sprintf("  %-*s %*s\n", lw, "sMAPE (%)",        vw, f2(smape)))
  }
  cat(sprintf("  %-*s %*s\n", lw, "AIC",           vw, f2(aic_val)))
  cat(sprintf("  %-*s %*s\n", lw, "AICc",          vw, f2(aicc_val)))
  cat(sprintf("  %-*s %*s\n", lw, "BIC",           vw, f2(bic_val)))
  cat(sep_m, "\n", sep = "")
  .print_rodape()

  invisible(data.frame(
    Metrica = c("R2", "R2_ajustado", "RMSE", "MAE", "MAPE", "sMAPE", "AIC", "AICc", "BIC"),
    Valor   = c(r2, r2_adj, rmse, mae, mape, smape, aic_val, aicc_val, bic_val)
  ))
}



#' @title Análise Residual
#' @description Realiza o diagnóstico completo dos resíduos com testes estatísticos formais,
#' interpretação automatizada e painel gráfico. Para modelos lineares (\code{lm}), avalia
#' normalidade (Shapiro-Wilk ou Kolmogorov-Smirnov para n > 5000), homocedasticidade
#' (Breusch-Pagan) e independência dos resíduos (Durbin-Watson). Para Modelos Lineares
#' Generalizados (\code{glm}), despacha automaticamente para o diagnóstico específico
#' Análise Residual Completa para Modelos Estatísticos
#'
#' @title Diagnóstico de Resíduos Polimórfico
#' @description Realiza o diagnóstico completo de resíduos com testes estatísticos
#' de Normalidade, Homocedasticidade e Independência, além de painel gráfico em ggplot2.
#' Suporta modelos de regressão linear (\code{lm}), delineamentos experimentais
#' (\code{estatR_dic}, \code{estatR_dbc}, \code{estatR_dql}, \code{estatR_fatorial}, \code{aov}),
#' modelos lineares generalizados (\code{glm}) e modelos de sobrevivência (\code{estatR_cox}, \code{coxph}).
#'
#' @param modelo Objeto de modelo ajustado do tipo \code{lm}, \code{glm}, \code{aov},
#' \code{estatR_dic}, \code{estatR_dbc}, \code{estatR_dql}, \code{estatR_fatorial} ou \code{estatR_cox}.
#' @param grafico Lógico. Se TRUE, exibe o painel gráfico de diagnóstico dos resíduos.
#' @return Retorna invisivelmente uma lista com os objetos dos testes estatísticos realizados.
#' @export
analise_residual <- function(modelo, grafico = TRUE) {
  # Suporte a Modelos Inflados de Zero (zeroinfl: ZIP e ZINB)
  if (inherits(modelo, "zeroinfl") && !is.null(modelo$dist)) {
    if (modelo$dist == "poisson" && exists(".analise_residual_poisson_iz", mode = "function")) {
      return(.analise_residual_poisson_iz(modelo, grafico = grafico))
    }
    if (modelo$dist == "negbin" && exists(".analise_residual_negbin_iz", mode = "function")) {
      return(.analise_residual_negbin_iz(modelo, grafico = grafico))
    }
  }

  # 0. Suporte a Modelos de Séries Temporais (previsao_serie, forecast, Arima, ets, stlm)
  if (inherits(modelo, c("previsao_serie", "forecast", "Arima", "ets", "HoltWinters", "stlm"))) {
    if (exists("analise_residual_serie", mode = "function")) {
      return(analise_residual_serie(modelo, grafico = grafico))
    }
  }

  # 1. Suporte a Aprendizado de Máquina (estatR_aprend_maq)
  if (inherits(modelo, "estatR_aprend_maq")) {
    if (exists("analise_residual_aprend_maq", mode = "function")) {
      return(analise_residual_aprend_maq(modelo, grafico = grafico))
    }
  }

  # 2. Suporte a Sobrevivência (Cox)
  if (inherits(modelo, c("estatR_cox", "coxph"))) {
    return(analise_residual_cox(modelo, grafico = grafico))
  }

  # 2. Suporte a Modelos de Experimentos do estatR (DIC, DBC, DQL, Fatorial, Subdividida, Fatorial Adicional)
  if (inherits(modelo, c("estatR_dic", "estatR_dbc", "estatR_dql", "estatR_fatorial", "estatR_subdividida", "estatR_fatorial_adicional"))) {
    if (!is.null(modelo$modelo_aov)) {
      modelo <- modelo$modelo_aov
    } else {
      stop("O objeto de experimento n\u00e3o cont\u00e9m o modelo ANOVA ajustado ($modelo_aov).")
    }
  }

  # 3. Suporte a Modelos Lineares Generalizados (GLM) e correlatos
  if (inherits(modelo, "betareg") && exists("analise_residual_beta", mode = "function")) {
    return(analise_residual_beta(modelo, grafico = grafico))
  }
  if (inherits(modelo, "negbin") && exists("analise_residual_binomial_negativa", mode = "function")) {
    return(analise_residual_binomial_negativa(modelo, grafico = grafico))
  }
  if (inherits(modelo, "glm") && !is.null(modelo$family)) {
    fam <- modelo$family$family
    if (fam == "binomial" && exists("analise_residual_logistica", mode = "function")) {
      return(analise_residual_logistica(modelo, grafico = grafico))
    }
    if (fam == "poisson" && exists("analise_residual_poisson", mode = "function")) {
      return(analise_residual_poisson(modelo, grafico = grafico))
    }
    if (fam == "Gamma" && exists("analise_residual_gama", mode = "function")) {
      return(analise_residual_gama(modelo, grafico = grafico))
    }
    if (grepl("Negative Binomial", fam, ignore.case = TRUE) && exists("analise_residual_binomial_negativa", mode = "function")) {
      return(analise_residual_binomial_negativa(modelo, grafico = grafico))
    }
  }

  # Suporte a objetos com classe 'aov' ou 'lm'
  res <- stats::residuals(modelo)
  if (is.null(res) || !is.numeric(res)) {
    if (!is.null(modelo$residuals) && is.numeric(modelo$residuals)) {
      res <- modelo$residuals
    } else {
      stop("N\u00e3o foi poss\u00edvel extrair res\u00edduos num\u00e9ricos do modelo fornecido.")
    }
  }
  n   <- length(res)

  # 1. Normalidade (Shapiro-Wilk para n <= 5000; Kolmogorov-Smirnov para n > 5000)
  if (n <= 5000) {
    nome_norm <- "Normalidade (Shapiro-Wilk)"
    st        <- shapiro.test(res)
    stat_norm <- sprintf("W = %-9.3f", st$statistic)
    p_norm    <- st$p.value
  } else {
    nome_norm <- "Normalidade (Kolmogorov-S)"
    res_pad   <- (res - mean(res)) / sd(res)
    st        <- suppressWarnings(stats::ks.test(res_pad, "pnorm"))
    stat_norm <- sprintf("D = %-9.3f", st$statistic)
    p_norm    <- st$p.value
  }

  # 2. Homocedasticidade (Breusch-Pagan, base R)
  res_sq       <- res^2
  dados_modelo <- modelo$model
  bp_mod       <- lm(res_sq ~ ., data = dados_modelo[, -1, drop = FALSE])
  r2_bp        <- summary(bp_mod)$r.squared
  bp_stat      <- n * r2_bp
  df_bp        <- length(modelo$coefficients) - 1
  p_bp         <- pchisq(bp_stat, df_bp, lower.tail = FALSE)

  # 3. Independência (Durbin-Watson, base R)
  diff_res <- diff(res)
  dw_stat  <- sum(diff_res^2) / sum(res^2)
  z_dw     <- (dw_stat - 2) / sqrt(4 / n)
  p_dw     <- 2 * pnorm(-abs(z_dw))

  formata_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("<0,001")
    return(formatC(p, format = "f", digits = 3, decimal.mark = ","))
  }

  linha <- paste0(strrep("\u2500", 56), "\n")

  .print_titulo("DIAGN\u00d3STICO DOS RES\u00cdDUOS")

  cat(sprintf("%-30s %-14s %s\n", "Teste", "Estat\u00edstica", "p-valor"))
  cat(linha)
  cat(sprintf("%-30s %-14s %s\n", nome_norm, stat_norm, formata_p(p_norm)))
  cat(sprintf("%-30s BP = %-9.3f %s\n", "Homocedasticidade (B-P)", bp_stat, formata_p(p_bp)))
  cat(sprintf("%-30s DW = %-9.3f %s\n", "Independ\u00eancia (Durbin-Watson)", dw_stat, formata_p(p_dw)))
  cat(linha)

  .print_topico("INTERPRETA\u00c7\u00c3O DO DIAGN\u00d3STICO")

  cor_ok   <- "\033[1;32m"
  cor_erro <- "\033[1;31m"
  reset    <- "\033[0m"

  if (p_norm < 0.05) {
    cat(sprintf("  %s[!]%s O teste de normalidade rejeitou a hip\u00f3tese nula\n    (p = %s). Os res\u00edduos n\u00e3o seguem uma distribui\u00e7\u00e3o normal,\n    o que pode afetar a confiabilidade dos intervalos de confian\u00e7a.\n\n",
                cor_erro, reset, formata_p(p_norm)))
  } else {
    cat(sprintf("  %s[\u2713]%s O teste n\u00e3o encontrou evid\u00eancias para rejeitar a normalidade\n    dos res\u00edduos (p = %s).\n\n",
                cor_ok, reset, formata_p(p_norm)))
  }

  if (p_bp < 0.05) {
    cat(sprintf("  %s[!]%s O teste detectou heterocedasticidade (p = %s),\n    indicando que a vari\u00e2ncia dos res\u00edduos n\u00e3o \u00e9 constante.\n\n",
                cor_erro, reset, formata_p(p_bp)))
  } else {
    cat(sprintf("  %s[\u2713]%s O teste indicou homocedasticidade (p = %s),\n    ou seja, a vari\u00e2ncia dos res\u00edduos \u00e9 considerada constante\n    ao longo dos valores ajustados.\n\n",
                cor_ok, reset, formata_p(p_bp)))
  }

  if (p_dw < 0.05) {
    direcao_ac <- if (dw_stat < 2) "positiva" else "negativa"
    cat(sprintf("  %s[!]%s O teste detectou autocorrela\u00e7\u00e3o %s (p = %s),\n    indicando depend\u00eancia entre os res\u00edduos (DW = %.3f).\n\n",
                cor_erro, reset, direcao_ac, formata_p(p_dw), dw_stat))
  } else {
    cat(sprintf("  %s[\u2713]%s O teste indicou independ\u00eancia dos res\u00edduos (p = %s,\n    DW = %.3f), ou seja, n\u00e3o h\u00e1 evid\u00eancia de autocorrela\u00e7\u00e3o.\n\n",
                cor_ok, reset, formata_p(p_dw), dw_stat))
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
        ggplot2::labs(title = "Res\u00edduos vs Ajustados",
                      x = "Valores Ajustados", y = "Res\u00edduos") +
        tema_painel

      p2 <- ggplot2::ggplot(df_res, ggplot2::aes(sample = res_pad)) +
        ggplot2::stat_qq(color = "#555555", alpha = 0.8) +
        ggplot2::stat_qq_line(color = "#D90429", linewidth = 0.7) +
        ggplot2::labs(title = "QQ-Plot Normal",
                      x = "Quantis Te\u00f3ricos", y = "Quantis Amostrais") +
        tema_painel

      p3 <- ggplot2::ggplot(df_res, ggplot2::aes(x = ajustados, y = sqrt(abs(res_pad)))) +
        ggplot2::geom_point(color = "#555555", alpha = 0.8) +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429",
                             linewidth = 0.7, formula = y ~ x) +
        ggplot2::labs(title = "Escala-Localiza\u00e7\u00e3o",
                      x = "Valores Ajustados", y = "\u221a|Res. Padronizados|") +
        tema_painel

      p4 <- ggplot2::ggplot(df_res, ggplot2::aes(x = res_pad)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                bins = 10, fill = "#555555", color = "black", alpha = 0.8) +
        ggplot2::geom_density(color = "#D90429", linewidth = 0.8) +
        ggplot2::labs(title = "Histograma dos Res\u00edduos",
                      x = "Res\u00edduos Padronizados", y = "Densidade",
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
    }, error = function(e) {
      message("[Aviso] N\u00e3o foi poss\u00edvel gerar o painel de res\u00edduos: ", e$message)
    })
  }

  .print_rodape()

  invisible(list(normalidade = st, shapiro = st,
                 bp = list(statistic = bp_stat, p.value = p_bp),
                 dw = list(statistic = dw_stat, p.value = p_dw)))
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Seleção Múltipla de Modelos (Best Subsets / Stepwise)
#' @description Avalia múltiplas combinações de variáveis preditoras para encontrar
#' os melhores modelos baseando-se em métricas de ajuste e parcimônia. Utiliza 
#' algoritmo exaustivo para até 10 preditores e heurística inteligente para mais.
#' @param formula Fórmula do modelo inicial (ex: y ~ .).
#' @param dados Data frame contendo os dados.
#' @param comparar_modelos Famílias de modelos a comparar: \code{"estatistico"} (padrão,
#'   avalia subconjuntos de variáveis da regressão estatística), \code{"aprend_maq"}
#'   ou ambos \code{c("estatistico", "aprend_maq")}.
#' @param top Número máximo de modelos a serem exibidos no ranking (padrão: 5).
#' @return Retorna invisivelmente um data frame com o ranking completo.
#' @export
selecao_modelos <- function(formula, dados, comparar_modelos = c("estatistico", "aprend_maq"), top = 5, ...) {
  # Polimorfismo: suporte direto a Séries Temporais (objeto 'ts' ou vetor numérico sem 'dados')
  if (stats::is.ts(formula) || (missing(dados) && is.numeric(formula))) {
    if (exists('selecao_modelos_serie', mode = 'function')) {
      dots <- list(...)
      h_val <- if ('h' %in% names(dots)) dots$h else 12
      fam_val <- if ('familia' %in% names(dots)) dots$familia else 'todas'
      top_val <- if ('top' %in% names(dots)) dots$top else 5
      ord_val <- if ('ordenar_por' %in% names(dots)) dots$ordenar_por else 'sMAPE'
      return(selecao_modelos_serie(serie = formula, h = h_val, familia = fam_val, top = top_val, ordenar_por = ord_val))
    }
  }
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")

  dots <- list(...)
  if ("familia" %in% names(dots)) {
    if (dots$familia == "poisson_iz" && exists(".selecao_modelos_poisson_iz", mode = "function")) {
      return(.selecao_modelos_poisson_iz(formula, dados, top = top))
    }
    if (dots$familia %in% c("negbin_iz", "binomial_negativa_iz") && exists(".selecao_modelos_negbin_iz", mode = "function")) {
      return(.selecao_modelos_negbin_iz(formula, dados, top = top))
    }
  }

  # Se o usuário solicitou comparar com ou exclusivamente Aprendizado de Máquina
  if (identical(comparar_modelos, "aprend_maq") ||
      (is.character(comparar_modelos) && "aprend_maq" %in% comparar_modelos && length(comparar_modelos) > 1)) {
    if (exists(".selecao_comparativa_modelos", mode = "function")) {
      return(.selecao_comparativa_modelos(formula = formula, dados = dados,
                                          comparar_modelos = comparar_modelos, top = top))
    }
  }

  # Preparação de variáveis

  resp_dep <- deparse(formula[[2]])
  if (grepl("Surv\\(", resp_dep)) {
    return(selecao_modelos_cox(formula, dados, top = top))
  }
  
  termos <- attr(terms(formula, data = dados), "term.labels")
  y_name <- as.character(formula[[2]])
  
  if (!(y_name %in% names(dados))) stop("Variável resposta não encontrada nos dados.")
  # Polimorfismo: detecta variável resposta binária -> logística
  Y_chk <- na.omit(dados[[y_name]])
  n_niveis_chk <- length(unique(as.character(Y_chk)))
  is_binario <- (n_niveis_chk == 2) && (is.factor(Y_chk) || is.character(Y_chk) ||
                  (is.numeric(Y_chk) && all(Y_chk %in% c(0, 1))))
  if (is_binario) {
    return(selecao_modelos_logistica(formula, dados, top = top))
  }
  
  # Polimorfismo: detecta contagem (Regressao Poisson)
  # Variavel resposta precisa ser numerica, inteira, nao negativa e com mais de 2 niveis
  is_count <- is.numeric(Y_chk) && all(Y_chk >= 0) && all(Y_chk == floor(Y_chk)) && n_niveis_chk > 2
  if (is_count) {
    return(selecao_modelos_poisson(formula, dados, top = top))
  }
  if (!is.numeric(dados[[y_name]])) stop("A variável resposta deve ser numérica.")
  
  # Limpeza de NAs (apenas colunas envolvidas) para garantir amostras idênticas
  cols_req <- c(y_name, termos)
  missing_cols <- setdiff(cols_req, names(dados))
  if (length(missing_cols) > 0) {
    stop("Variáveis preditoras não encontradas: ", paste(missing_cols, collapse = ", "))
  }
  
  df_clean <- na.omit(dados[, cols_req, drop = FALSE])
  n <- nrow(df_clean)
  if (n < 5) stop("Poucas observações válidas após remover NAs (mínimo 5).")
  
  p_total <- length(termos)
  if (p_total == 0) stop("O modelo deve ter pelo menos 1 variável preditora.")
  
  usou_heuristica <- FALSE
  termos_originais <- termos
  
  # Algoritmo Inteligente (Heurística de Screening se p > 10)
  if (p_total > 10) {
    usou_heuristica <- TRUE
    num_termos <- termos[sapply(df_clean[termos], is.numeric)]
    y_num <- df_clean[[y_name]]
    
    cor_vals <- sapply(num_termos, function(x) {
      abs(cor(y_num, df_clean[[x]], use = "complete.obs"))
    })
    
    top_10_names <- names(sort(cor_vals, decreasing = TRUE))[1:min(10, length(cor_vals))]
    cat_termos <- setdiff(termos, num_termos)
    termos <- unique(c(top_10_names, cat_termos))
    if (length(termos) > 10) termos <- termos[1:10]
  }
  
  p_avaliar <- length(termos)
  
  # Gerador Exaustivo de Modelos
  modelos_list <- list()
  for (k in 1:p_avaliar) {
    comb <- combn(termos, k, simplify = FALSE)
    modelos_list <- c(modelos_list, comb)
  }
  
  total_avaliados <- length(modelos_list)
  
  # Função interna de VIF
  .calc_vif <- function(mod) {
    if (length(coef(mod)) <= 2) return(1.0)
    X <- model.matrix(mod)[, -1, drop = FALSE]
    if (ncol(X) <= 1) return(1.0)
    
    vars <- apply(X, 2, var)
    if (any(vars == 0)) return(Inf)
    
    vif_vals <- tryCatch({
      diag(solve(cor(X)))
    }, error = function(e) Inf)
    return(max(vif_vals))
  }
  
  # Pre-allocate
  res_mod <- character(total_avaliados)
  res_k <- integer(total_avaliados)
  res_r2 <- numeric(total_avaliados)
  res_r2adj <- numeric(total_avaliados)
  res_aic <- numeric(total_avaliados)
  res_aicc <- numeric(total_avaliados)
  res_vif <- numeric(total_avaliados)
  
  for (i in seq_len(total_avaliados)) {
    preds <- modelos_list[[i]]
    f_mod <- as.formula(paste(y_name, "~", paste(preds, collapse = " + ")))
    
    mod <- tryCatch(lm(f_mod, data = df_clean), error = function(e) NULL)
    if (is.null(mod)) next
    
    sm <- summary(mod)
    aic_val <- AIC(mod)
    
    K_params <- length(coef(mod)) + 1
    if (n - K_params - 1 > 0) {
      aicc_val <- aic_val + (2 * K_params * (K_params + 1)) / (n - K_params - 1)
    } else {
      aicc_val <- Inf
    }
    
    res_mod[i] <- paste(preds, collapse = " + ")
    res_k[i] <- length(preds)
    res_r2[i] <- if (is.null(sm$r.squared)) NA else sm$r.squared
    res_r2adj[i] <- if (is.null(sm$adj.r.squared)) NA else sm$adj.r.squared
    res_aic[i] <- aic_val
    res_aicc[i] <- aicc_val
    res_vif[i] <- .calc_vif(mod)
  }
  
  # Montar dataframe de resultados
  df_res <- data.frame(
    Modelo = res_mod,
    k = res_k,
    R2 = res_r2,
    R2_adj = res_r2adj,
    AIC = res_aic,
    AICc = res_aicc,
    VIF_Max = res_vif,
    stringsAsFactors = FALSE
  )
  
  df_res <- df_res[!is.na(df_res$R2_adj), ]
  df_res <- df_res[order(df_res$R2_adj, decreasing = TRUE), ]
  
  df_top <- df_res
  if (nrow(df_top) > top) df_top <- df_top[1:top, ]
  
  # Formatando para exibição
  df_fmt <- data.frame(
    Ranking = paste0(1:nrow(df_top), "\u00BA"),
    Modelo = df_top$Modelo,
    k = as.character(df_top$k),
    "R²" = sapply(df_top$R2, .fmt_pct, decimais = 1),
    "R² Adj" = sapply(df_top$R2_adj, .fmt_pct, decimais = 1),
    "AIC" = sapply(df_top$AIC, .fmt_num, decimais = 1),
    "AICc" = sapply(df_top$AICc, .fmt_num, decimais = 1),
    "VIF Máx" = sapply(df_top$VIF_Max, function(v) if(is.infinite(v) || is.na(v)) "-" else .fmt_num(v, decimais = 1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  # Imprimir
  .print_titulo("SELEÇÃO DE MODELOS DE REGRESSÃO")
  
  cat(sprintf("  Variável Resposta: %s\n", y_name))
  if (usou_heuristica) {
    cat(sprintf("  Preditores base:   %s\n", paste(termos_originais, collapse = ", ")))
    cat(sprintf("  Método utilizado:  Heurístico (Screening TOP 10 -> %d modelos avaliados)\n", total_avaliados))
  } else {
    cat(sprintf("  Preditores base:   %s\n", paste(termos, collapse = ", ")))
    cat(sprintf("  Método utilizado:  Exaustivo (%d modelos avaliados)\n", total_avaliados))
  }
  cat("\n")
  
  .print_topico(sprintf("TOP %d MODELOS (Ordenados por R² Ajustado)", nrow(df_fmt)))
  
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_fmt, align = c("left", "left", "center", "center", "center", "center", "center", "center"))
  } else {
    print(df_fmt, row.names = FALSE)
  }
  
  cat("  * k = Número de variáveis preditoras no modelo.\n")
  cat("  * VIF Máx = Maior fator de inflação da variância no modelo (ideal < 5).\n\n")
  
  .print_rodape()
  
  invisible(df_res)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Diagnóstico de Alavancagem e Pontos Influentes
#' @description Identifica observações atípicas, pontos de alta alavancagem (leverage) e
#'   observações influentes que afetam desproporcionalmente a estimação dos parâmetros
#'   do modelo de regressão linear (\code{lm}) ou modelo linear generalizado (\code{glm}).
#'   Opcionalmente reestima o modelo omitindo tais observações para comparar o impacto
#'   nos coeficientes e exibe gráficos diagnósticos com grade dupla.
#'
#' @param modelo Objeto de modelo ajustado da classe \code{lm} ou \code{glm}.
#' @param corte Critério para limites críticos: \code{"automatico"} (padrão:
#'   Cook > 4/n, Leverage > 2p/n, |DFFITS| > 2*sqrt(p/n)) ou \code{"conservador"}
#'   (Cook > 1, Leverage > 3p/n, |DFFITS| > 3*sqrt(p/n)).
#' @param comparar Lógico. Se TRUE (padrão) e forem detectados pontos influentes,
#'   reajusta o modelo omitindo tais observações e exibe a comparação dos coeficientes.
#' @param grafico Lógico. Se TRUE (padrão), plota diagnóstico visual com \code{grade = "dupla"}.
#'   Para regressão linear simples (1 preditor), sobrepõe as retas ajustadas com e sem
#'   os pontos influentes. Para regressão múltipla ou GLM, gera gráfico de dispersão
#'   Alavancagem vs. Distância de Cook.
#'
#' @return Retorna invisivelmente uma lista com:
#'   \item{tabela_pontos}{Data frame com as observações destacadas e suas métricas.}
#'   \item{n_influentes}{Número de observações influentes encontradas.}
#'   \item{indices_influentes}{Vetor de índices das observações influentes.}
#'   \item{comparacao_coefs}{Data frame comparativo de coeficientes (quando aplicável).}
#'   \item{modelo_sem_influentes}{Modelo reajustado sem os pontos influentes (quando aplicável).}
#' @export
pontos_influentes <- function(modelo,
                              corte = c("automatico", "conservador"),
                              comparar = TRUE,
                              grafico = TRUE) {
  if (!inherits(modelo, c("lm", "glm"))) {
    stop("O argumento 'modelo' deve ser um objeto da classe 'lm' ou 'glm'.")
  }
  corte <- match.arg(corte)

  n <- stats::nobs(modelo)
  p <- length(stats::coef(modelo))

  if (n <= p + 1) {
    stop("Número insuficiente de graus de liberdade residuais para cálculo de influência.")
  }

  hat_vals  <- stats::hatvalues(modelo)
  cooks_d   <- stats::cooks.distance(modelo)
  res_stud  <- tryCatch(stats::rstudent(modelo), error = function(e) stats::rstandard(modelo))

  dffits_v  <- tryCatch({
    stats::dffits(modelo)
  }, error = function(e) {
    res_stud * sqrt(hat_vals / pmax(1e-10, 1 - hat_vals))
  })

  if (corte == "automatico") {
    lim_hat   <- 2 * p / n
    lim_cook  <- 4 / n
    lim_dff   <- 2 * sqrt(p / n)
    lim_res   <- 2.0
  } else {
    lim_hat   <- 3 * p / n
    lim_cook  <- 1.0
    lim_dff   <- 3 * sqrt(p / n)
    lim_res   <- 2.5
  }

  is_cook  <- !is.na(cooks_d) & (cooks_d > lim_cook)
  is_dff   <- !is.na(dffits_v) & (abs(dffits_v) > lim_dff)
  is_infl  <- is_cook | is_dff
  is_alav  <- !is.na(hat_vals) & (hat_vals > lim_hat)
  is_outl  <- !is.na(res_stud) & (abs(res_stud) > lim_res)

  classificacao <- character(n)
  for (i in seq_len(n)) {
    if (is_infl[i] && is_alav[i]) {
      classificacao[i] <- "Alavancagem e Influ\u00eancia"
    } else if (is_infl[i]) {
      classificacao[i] <- "Ponto Influente"
    } else if (is_alav[i]) {
      classificacao[i] <- "Alta Alavancagem"
    } else if (is_outl[i]) {
      classificacao[i] <- "Res\u00edduo At\u00edpico"
    } else {
      classificacao[i] <- "T\u00edpico"
    }
  }

  idx_destaque <- which(classificacao != "T\u00edpico")
  idx_infl_pur <- which(is_infl)
  n_infl       <- length(idx_infl_pur)

  nomes_obs <- names(hat_vals)
  if (is.null(nomes_obs)) nomes_obs <- as.character(seq_len(n))

  df_todos <- data.frame(
    Indice        = seq_len(n),
    Obs           = nomes_obs,
    Resid_Stud    = round(res_stud, 3),
    Leverage      = round(hat_vals, 3),
    Cook_D        = round(cooks_d, 4),
    DFFITS        = round(dffits_v, 3),
    Classificacao = classificacao,
    stringsAsFactors = FALSE
  )

  df_destaque <- df_todos[idx_destaque, , drop = FALSE]
  if (nrow(df_destaque) > 0) {
    df_destaque <- df_destaque[order(df_destaque$Cook_D, decreasing = TRUE), ]
  }

  tabela_comp  <- NULL
  modelo_sem   <- NULL

  if (comparar && n_infl > 0 && (n - n_infl > p)) {
    dados_orig <- tryCatch(stats::model.frame(modelo), error = function(e) NULL)
    if (!is.null(dados_orig)) {
      dados_sem <- dados_orig[-idx_infl_pur, , drop = FALSE]
      modelo_sem <- tryCatch(stats::update(modelo, data = dados_sem), error = function(e) NULL)

      if (!is.null(modelo_sem)) {
        cf_com <- stats::coef(modelo)
        cf_sem <- stats::coef(modelo_sem)
        sm_sem <- summary(modelo_sem)
        p_sem  <- sm_sem$coefficients[, 4]

        var_pct <- ((cf_sem - cf_com) / pmax(1e-10, abs(cf_com))) * 100

        tabela_comp <- data.frame(
          "Termo"         = names(cf_com),
          "Modelo Completo" = sapply(cf_com, .fmt_num, decimais = 3),
          "Sem Influentes"= sapply(cf_sem, .fmt_num, decimais = 3),
          "Variação (%)"  = sprintf("%+.1f%%", var_pct),
          "p-valor (Sem)" = sapply(p_sem, function(pv) if (pv < 0.001) "< 0,001" else .fmt_p(pv)),
          stringsAsFactors = FALSE, check.names = FALSE
        )
      }
    }
  }

  .print_titulo("DIAGN\u00d3STICO DE ALAVANCAGEM E PONTOS INFLUENTES")
  cat(sprintf("  Modelo avaliado:   %s (%d observa\u00e7\u00f5es, %d par\u00e2metros)\n",
              class(modelo)[1], n, p))
  cat(sprintf("  Crit\u00e9rio de corte:  %s\n", ifelse(corte == "automatico", "Autom\u00e1tico (Literatura Padr\u00e3o)", "Conservador")))
  cat(sprintf("  Limites de Alerta: Cook > %s | Alavancagem > %s | |DFFITS| > %s\n\n",
              .fmt_num(lim_cook, 4), .fmt_num(lim_hat, 3), .fmt_num(lim_dff, 3)))

  .print_topico(sprintf("OBSERVA\u00c7\u00d5ES DESTACADAS (%d detectadas)", nrow(df_destaque)))

  if (nrow(df_destaque) > 0) {
    df_print <- data.frame(
      "Obs."           = df_destaque$Obs,
      "Res\u00edduo Stud." = sprintf("%.3f", df_destaque$Resid_Stud),
      "Alavancagem (h)" = sprintf("%.3f", df_destaque$Leverage),
      "Cook (D)"       = sprintf("%.4f", df_destaque$Cook_D),
      "DFFITS"         = sprintf("%.3f", df_destaque$DFFITS),
      "Diagn\u00f3stico" = df_destaque$Classificacao,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    .print_tabela_estatR(head(df_print, 15), align = c("center", "center", "center", "center", "center", "left"))
    if (nrow(df_print) > 15) {
      cat(sprintf("  * Exibindo as 15 principais de %d observa\u00e7\u00f5es at\u00edpicas ordenadas por Cook (D).\n\n", nrow(df_print)))
    }
  } else {
    cat("  Nenhuma observa\u00e7\u00e3o at\u00edpica ou influente ultrapassou os limites de alerta.\n\n")
  }

  if (!is.null(tabela_comp)) {
    .print_topico(sprintf("IMPACTO NOS COEFICIENTES (Removendo %d ponto(s) influente(s))", n_infl))
    .print_tabela_estatR(tabela_comp, align = c("left", "center", "center", "center", "center"))
  }

  cor_ok   <- "\033[1;32m"
  cor_warn <- "\033[1;33m"
  cor_erro <- "\033[1;31m"
  reset    <- "\033[0m"

  .print_topico("CONCLUS\u00c3O DO DIAGN\u00d3STICO")
  if (n_infl == 0) {
    cat(sprintf("  %s[\u2713]%s Aus\u00eancia de pontos influentes desestabilizadores. As estimativas s\u00e3o robustas.\n\n", cor_ok, reset))
  } else if (n_infl <= max(1, floor(0.02 * n))) {
    cat(sprintf("  %s[!]%s Detectado(s) %d ponto(s) influente(s) com impacto nos coeficientes. Recomenda-se investigar a origem destas observa\u00e7\u00f5es.\n\n", cor_warn, reset, n_infl))
  } else {
    cat(sprintf("  %s[!]%s Detectada presen\u00e7a expressiva de observa\u00e7\u00f5es influentes (%d obs). Considere modelos robustos ou investiga\u00e7\u00e3o aprofundada.\n\n", cor_erro, reset, n_infl))
  }

  .print_rodape()

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch({
      tema_padrao <- if (exists("meu_tema", mode = "function")) {
        meu_tema(grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }

      termos <- attr(stats::terms(modelo), "term.labels")
      is_reg_simples <- inherits(modelo, "lm") && length(termos) == 1 &&
                        is.numeric(stats::model.frame(modelo)[[termos[1]]])

      if (is_reg_simples && !is.null(modelo_sem)) {
        mf <- stats::model.frame(modelo)
        y_col <- names(mf)[1]
        x_col <- termos[1]

        df_scatter <- data.frame(
          x = mf[[x_col]],
          y = mf[[y_col]],
          Tipo = ifelse(seq_len(n) %in% idx_infl_pur, "Ponto Influente", "Observação Típica"),
          Obs = nomes_obs
        )

        p_plot <- ggplot2::ggplot(df_scatter, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_point(ggplot2::aes(color = Tipo, shape = Tipo, size = Tipo), alpha = 0.85) +
          ggplot2::scale_color_manual(values = c("Ponto Influente" = "#D90429", "Observação Típica" = "#2B5C8F")) +
          ggplot2::scale_shape_manual(values = c("Ponto Influente" = 17, "Observação Típica" = 16)) +
          ggplot2::scale_size_manual(values = c("Ponto Influente" = 3.5, "Observação Típica" = 2.2)) +
          ggplot2::geom_abline(intercept = stats::coef(modelo)[1], slope = stats::coef(modelo)[2],
                               color = "#2B5C8F", linewidth = 1.0, linetype = "solid") +
          ggplot2::geom_abline(intercept = stats::coef(modelo_sem)[1], slope = stats::coef(modelo_sem)[2],
                               color = "#D90429", linewidth = 1.1, linetype = "dashed") +
          ggplot2::labs(
            title = "Ajuste Linear: Modelo Completo vs. Sem Pontos Influentes",
            subtitle = sprintf("Linha contínua azul: Completo | Linha tracejada vermelha: Sem %d influente(s)", n_infl),
            x = x_col, y = y_col
          ) +
          tema_padrao

        print(p_plot)

      } else {
        df_diag_plot <- data.frame(
          Leverage = hat_vals,
          Cook     = cooks_d,
          Tipo     = df_todos$Classificacao,
          Obs      = nomes_obs
        )
        df_diag_plot$Destaque <- df_diag_plot$Tipo != "Típico"

        p_plot <- ggplot2::ggplot(df_diag_plot, ggplot2::aes(x = Leverage, y = Cook)) +
          ggplot2::geom_hline(yintercept = lim_cook, linetype = "dashed", color = "#D90429", linewidth = 0.8) +
          ggplot2::geom_vline(xintercept = lim_hat, linetype = "dashed", color = "#E69F00", linewidth = 0.8) +
          ggplot2::geom_point(ggplot2::aes(color = Tipo, size = Cook), alpha = 0.8) +
          ggplot2::scale_color_manual(values = c(
            "Alavancagem e Influência" = "#D90429",
            "Ponto Influente"          = "#E69F00",
            "Alta Alavancagem"         = "#56B4E9",
            "Resíduo Atípico"          = "#CC79A7",
            "Típico"                   = "#999999"
          )) +
          ggplot2::labs(
            title = "Alavancagem (Leverage) vs. Distância de Cook",
            subtitle = sprintf("Linhas pontilhadas: Limite Cook (%.4f) e Limite Alavancagem (%.3f)", lim_cook, lim_hat),
            x = "Alavancagem (hatvalues)", y = "Distância de Cook"
          ) +
          tema_padrao

        print(p_plot)
      }
    }, error = function(e) {
      message("[Aviso] Não foi possível renderizar o gráfico diagnóstico: ", e$message)
    })
  }

  resultado <- list(
    tabela_pontos         = df_destaque,
    n_influentes          = n_infl,
    indices_influentes    = idx_infl_pur,
    comparacao_coefs      = tabela_comp,
    modelo_sem_influentes = modelo_sem
  )
  invisible(resultado)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Diagnóstico de Multicolinearidade (VIF e GVIF)
#' @description Avalia a presença e a gravidade de multicolinearidade entre as variáveis
#'   preditoras de um modelo linear (\code{lm}) ou modelo linear generalizado (\code{glm}).
#'   Calcula o Fator de Inflação da Variância clássico (VIF) para variáveis contínuas
#'   e o VIF Generalizado (GVIF) para fatores com mais de dois níveis (Fox & Monette, 1992),
#'   além da Tolerância e do Número de Condição da matriz de planejamento.
#'
#' @param modelo Objeto ajustado da classe \code{lm} ou \code{glm}.
#' @param grafico Lógico. Se TRUE, plota gráfico de barras horizontais do VIF/GVIF
#'   com linhas de corte em 5 (moderado) e 10 (severo) com \code{grade = "dupla"}.
#'
#' @return Retorna invisivelmente um data frame com as métricas de multicolinearidade.
#' @export
multicolinearidade <- function(modelo, grafico = FALSE) {
  if (!inherits(modelo, c("lm", "glm"))) {
    stop("O argumento 'modelo' deve ser um objeto da classe 'lm' ou 'glm'.")
  }

  termos_obj <- stats::terms(modelo)
  term_labels <- attr(termos_obj, "term.labels")

  if (length(term_labels) < 2) {
    message("[Informação] Modelos com apenas uma variável preditora não apresentam multicolinearidade.")
    return(invisible(NULL))
  }

  X <- stats::model.matrix(modelo)
  tem_intercepto <- attr(termos_obj, "intercept") == 1

  assign_vec <- attr(X, "assign")
  if (tem_intercepto) {
    X <- X[, -1, drop = FALSE]
    assign_vec <- assign_vec[-1]
  }

  p_cols <- ncol(X)
  if (p_cols < 2) {
    message("[Informação] Menos de duas colunas de preditores na matriz de planejamento.")
    return(invisible(NULL))
  }

  sd_cols <- apply(X, 2, stats::sd)
  if (any(sd_cols == 0)) {
    warning("Colunas constantes detectadas entre os preditores.")
  }

  R <- tryCatch(stats::cor(X), error = function(e) NULL)
  det_R <- if (!is.null(R)) tryCatch(det(R), error = function(e) 0) else 0

  u_assign <- sort(unique(assign_vec))
  n_termos <- length(u_assign)

  df_list <- list()

  for (i in seq_along(u_assign)) {
    as_val <- u_assign[i]
    nome_termo <- term_labels[as_val]
    cols_termo <- which(assign_vec == as_val)
    df_termo   <- length(cols_termo)

    if (df_termo == 1 && !is.null(R) && det_R > 1e-12) {
      inv_R <- tryCatch(solve(R), error = function(e) NULL)
      if (!is.null(inv_R)) {
        vif_val <- inv_R[cols_termo, cols_termo]
      } else {
        m_aux <- stats::lm(X[, cols_termo] ~ X[, -cols_termo, drop = FALSE])
        r2_aux <- summary(m_aux)$r.squared
        vif_val <- if (r2_aux < 0.99999) 1 / (1 - r2_aux) else 999
      }
      gvif_val <- vif_val
      gvif_adj <- sqrt(vif_val)
      vif_equiv <- vif_val
    } else if (df_termo > 1 && !is.null(R) && det_R > 1e-12) {
      R11 <- R[cols_termo, cols_termo, drop = FALSE]
      R22 <- R[-cols_termo, -cols_termo, drop = FALSE]
      det_R11 <- det(R11)
      det_R22 <- det(R22)
      gvif_val <- (det_R11 * det_R22) / det_R
      gvif_adj <- gvif_val^(1 / (2 * df_termo))
      vif_equiv <- gvif_adj^2
    } else {
      gvif_val <- Inf
      gvif_adj <- Inf
      vif_equiv <- Inf
    }

    tol_val <- if (is.finite(vif_equiv) && vif_equiv > 0) 1 / vif_equiv else 0

    diag_txt <- if (vif_equiv < 5) {
      "Baixa (Ideal)"
    } else if (vif_equiv < 10) {
      "Moderada (Atenção)"
    } else {
      "Severa (Crítica)"
    }

    df_list[[i]] <- data.frame(
      Termo        = nome_termo,
      gl           = df_termo,
      VIF_GVIF     = gvif_val,
      GVIF_Ajust   = gvif_adj,
      VIF_Equiv    = vif_equiv,
      Tolerancia   = tol_val,
      Diagnostico  = diag_txt,
      stringsAsFactors = FALSE
    )
  }

  tabela_vif <- do.call(rbind, df_list)

  X_scaled <- scale(X, center = FALSE, scale = sqrt(colSums(X^2)))
  sv_vals  <- tryCatch(svd(X_scaled)$d, error = function(e) c(1, 1))
  kappa_val <- if (min(sv_vals) > 1e-12) max(sv_vals) / min(sv_vals) else Inf

  .print_titulo("DIAGN\u00d3STICO DE MULTICOLINEARIDADE")
  cat(sprintf("  Modelo avaliado:     %s (%d preditores na f\u00f3rmula)\n", class(modelo)[1], n_termos))
  cat(sprintf("  N\u00famero de Condi\u00e7\u00e3o: \u03ba = %s (%s)\n\n",
              if (is.infinite(kappa_val)) "Inf" else .fmt_num(kappa_val, 2),
              if (kappa_val < 15) "estabilidade global excelente" else
              if (kappa_val < 30) "colinearidade global moderada" else "colinearidade severa"))

  .print_topico("FATORES DE INFLA\u00c7\u00c3O DA VARI\u00c2NCIA (VIF / GVIF)")

  tem_fator_mult <- any(tabela_vif$gl > 1)

  df_print <- if (tem_fator_mult) {
    data.frame(
      'Preditor'        = tabela_vif$Termo,
      'gl'              = tabela_vif$gl,
      'GVIF'            = sapply(tabela_vif$VIF_GVIF, function(v) if (is.infinite(v)) 'Inf' else .fmt_num(v, 2)),
      'GVIF^(1/(2*gl))' = sapply(tabela_vif$GVIF_Ajust, function(v) if (is.infinite(v)) 'Inf' else .fmt_num(v, 3)),
      'VIF Equivalente' = sapply(tabela_vif$VIF_Equiv, function(v) if (is.infinite(v)) 'Inf' else .fmt_num(v, 2)),
      'Tolerância'      = sapply(tabela_vif$Tolerancia, .fmt_num, decimais = 3),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  } else {
    data.frame(
      'Preditor'   = tabela_vif$Termo,
      'VIF'        = sapply(tabela_vif$VIF_GVIF, function(v) if (is.infinite(v)) 'Inf' else .fmt_num(v, 2)),
      'Tolerância' = sapply(tabela_vif$Tolerancia, .fmt_num, decimais = 3),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }

  align_vec <- if (tem_fator_mult) c('left', 'center', 'center', 'center', 'center', 'center') else c('left', 'center', 'center')
  .print_tabela_estatR(df_print, align = align_vec)

  cat("  * Crit\u00e9rios: VIF < 5 (Ideal) | 5 a 10 (Moderado / Alerta) | >= 10 (Severo / Cr\u00edtico).\n")
  if (tem_fator_mult) {
    cat("  * Para fatores com gl > 1, avalia-se GVIF^(1/(2*gl)) comparado \u00e0 raiz do VIF padr\u00e3o.\n")
  }
  cat("\n")

  cor_ok   <- "\033[1;32m"
  cor_warn <- "\033[1;33m"
  cor_erro <- "\033[1;31m"
  reset    <- "\033[0m"

  max_vif <- max(tabela_vif$VIF_Equiv)
  .print_topico("CONCLUS\u00c3O DO DIAGN\u00d3STICO")
  if (max_vif < 5 && kappa_val < 15) {
    cat(sprintf("  %s[\u2713]%s Aus\u00eancia de multicolinearidade prejudicial. O modelo preserva estabilidade param\u00e9trica.\n\n", cor_ok, reset))
  } else if (max_vif < 10 && kappa_val < 30) {
    cat(sprintf("  %s[!]%s Multicolinearidade moderada detectada. Recomenda-se aten\u00e7\u00e3o aos erros padr\u00e3o dos preditores destacados.\n\n", cor_warn, reset))
  } else {
    cat(sprintf("  %s[!]%s Multicolinearidade severa detectada (VIF >= 10). Recomenda-se eliminar preditores redundantes ou utilizar regulariza\u00e7\u00e3o (Lasso/Ridge).\n\n", cor_erro, reset))
  }

  .print_rodape()

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch({
      tema_padrao <- if (exists("meu_tema", mode = "function")) {
        meu_tema(grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }

      df_g <- tabela_vif
      df_g$Termo <- factor(df_g$Termo, levels = rev(df_g$Termo))

      p_vif <- ggplot2::ggplot(df_g, ggplot2::aes(x = Termo, y = VIF_Equiv)) +
        ggplot2::geom_col(ggplot2::aes(fill = Diagnostico), width = 0.6) +
        ggplot2::geom_hline(yintercept = 5, linetype = "dashed", color = "#E69F00", linewidth = 0.8) +
        ggplot2::geom_hline(yintercept = 10, linetype = "dashed", color = "#D90429", linewidth = 0.8) +
        ggplot2::scale_fill_manual(values = c(
          "Baixa (Ideal)"       = "#2B5C8F",
          "Moderada (Atenção)"  = "#E69F00",
          "Severa (Crítica)"    = "#D90429"
        )) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Fator de Inflação da Variância (VIF / GVIF)",
          subtitle = "Linhas pontilhadas: Limite de Alerta (5) e Limite Crítico (10)",
          x = "", y = "VIF Equivalente"
        ) +
        tema_padrao

      print(p_vif)
    }, error = function(e) {
      message("[Aviso] Não foi possível renderizar o gráfico de VIF: ", e$message)
    })
  }

  invisible(tabela_vif)
}
