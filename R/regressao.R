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
  tipo_reg <- if (length(x_names) > 1) "M\u00daNTIPLA" else "SIMPLES"
  cat(sprintf("\n\u2500\u2500 RESULTADOS DA REGRESS\u00c3O %s ", tipo_reg))
  cat(strrep("\u2500", max(0, 55 - nchar(tipo_reg))), "\n", sep = "")
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

  aj_row("R\u00b2 (Explica\u00e7\u00e3o):", paste0(formatC(r2     * 100, format="f", digits=2, decimal.mark=","), "%"))
  aj_row("R\u00b2 ajustado:",  paste0(formatC(r2_adj * 100, format="f", digits=2, decimal.mark=","), "%"))
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

  cat("\n", strrep("\u2500", w_total + 2), "\n", sep = "")

  # ── 10. Gráfico (regressão simples) ───────────────────────────────────────
  if (grafico && length(x_names) == 1) {
    tryCatch({
      if (exists("grafico_de_dispersao") && exists("meu_tema")) {
        intercepto <- coefs[1, 1]
        inclinacao <- coefs[2, 1]
        sinal      <- if (inclinacao >= 0) "+" else "-"
        equacao_texto <- sprintf("%s = %.2f %s %.2f * %s\nR\u00b2 = %.1f%%",
                                 y_name, intercepto, sinal, abs(inclinacao), x_names[1], r2 * 100)

        if (inclinacao < 0) {
          pos_x <- Inf;  pos_y <- Inf; anc_h <- 1.1;  anc_v <- 1.5
        } else {
          pos_x <- -Inf; pos_y <- Inf; anc_h <- -0.1; anc_v <- 1.5
        }

        chamada_grafico <- bquote(
          grafico_de_dispersao(data = dados, x = .(as.name(x_names[1])), y = .(as.name(y_name)),
                               reta = TRUE, ic = TRUE, outlier = TRUE,
                               deteccao_outliers = "residuos", paleta = 1)
        )

        old_w <- getOption("warn"); options(warn = -1)
        p_plot <- ggplot2::ggplot() + eval(chamada_grafico) +
          ggplot2::annotate("text", x = pos_x, y = pos_y, label = equacao_texto,
                            hjust = anc_h, vjust = anc_v, size = 6.5,
                            fontface = "bold", color = "#333333") +
          ggplot2::labs(title = sprintf("Regress\u00e3o Linear: %s vs %s", y_name, x_names[1]),
                        caption = "estatR") +
          meu_tema(grade = "dupla")
        suppressMessages(print(p_plot))
        options(warn = old_w)
      }
    }, error = function(e) {
      message("\n[Aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico: ", e$message)
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
  cat(sprintf("\n\u2500\u2500 PREVIS\u00d5ES DO MODELO (Vari\u00e1vel: %s) \u2500\u2500\n", y_name))

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

  # ── Gráfico (regressão simples + intervalo) ───────────────────────────────
  if (grafico && length(x_names) == 1 && intervalo != "nenhum") {
    tryCatch({
      if (exists("grafico_de_dispersao") && exists("meu_tema")) {
        dados_orig <- modelo$model

        chamada_grafico <- bquote(
          grafico_de_dispersao(data = dados_orig,
                               x = .(as.name(x_names[1])),
                               y = .(as.name(y_name)),
                               reta = TRUE, ic = TRUE,
                               nivel_ic = .(nivel), paleta = 1)
        )

        old_w <- getOption("warn"); options(warn = -1)
        p_plot <- ggplot2::ggplot() + eval(chamada_grafico) +
          ggplot2::geom_point(
            data = df_final,
            ggplot2::aes(x = .data[[x_names[1]]], y = Previsao),
            color = "#D90429", size = 4.5, shape = 19
          ) +
          ggplot2::geom_errorbar(
            data = df_final,
            ggplot2::aes(x = .data[[x_names[1]]], ymin = LI, ymax = LS),
            color = "#D90429", width = 0.2, linewidth = 1
          ) +
          ggplot2::labs(
            title    = sprintf("Predi\u00e7\u00e3o: %s vs %s", y_name, x_names[1]),
            subtitle = sprintf("Intervalo de %s (%.0f%%)", nome_int, nivel * 100),
            caption  = "estatR"
          ) +
          meu_tema(grade = "dupla")
        suppressMessages(print(p_plot))
        options(warn = old_w)
      }
    }, error = function(e) {
      message("\n[Aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico: ", e$message)
    })
  }

  invisible(df_final)
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

  lw   <- 20
  vw   <- 10
  linha <- paste0(strrep("\u2500", lw + vw + 1), "\n")

  f3 <- function(x) formatC(x, format = "f", digits = 3, decimal.mark = ",")
  f2 <- function(x) formatC(x, format = "f", digits = 2, decimal.mark = ",")

  cat("\nM\u00c9TRICAS E CRIT\u00c9RIOS DE SELE\u00c7\u00c3O\n\n")
  cat(linha)
  cat(sprintf("%-*s %*s\n", lw, "M\u00e9trica", vw, "Valor"))
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

  cat("\nDIAGN\u00d3STICO DOS RES\u00cdDUOS\n\n")
  cat(linha)
  cat(sprintf("%-28s %-14s %s\n", "Teste", "Estat\u00edstica", "p-valor"))
  cat(linha)
  cat(sprintf("%-28s W = %-10.3f %s\n", "Normalidade (Shapiro-Wilk)", w_norm, formata_p(p_norm)))
  cat(sprintf("%-28s BP = %-9.3f %s\n", "Homocedasticidade (B-P)",   bp_stat, formata_p(p_bp)))
  cat(linha)

  cat("\nINTERPRETAÇÃO DO DIAGNÓSTICO\n\n")

  if (p_norm < 0.05) {
    cat(sprintf("[!] : O teste de normalidade rejeitou a hip\u00f3tese nula\n    (p = %s). Os res\u00edduos n\u00e3o seguem uma distribui\u00e7\u00e3o normal,\n    o que pode afetar a confiabilidade dos intervalos de confian\u00e7a.\n\n",
                formata_p(p_norm)))
  } else {
    cat(sprintf("[\u2713] : O teste n\u00e3o encontrou evid\u00eancias para rejeitar a normalidade\n    dos res\u00edduos (p = %s).\n\n",
                formata_p(p_norm)))
  }

  if (p_bp < 0.05) {
    cat(sprintf("[!] : O teste detectou heterocedasticidade (p = %s),\n    indicando que a vari\u00e2ncia dos res\u00edduos n\u00e3o \u00e9 constante.\n\n",
                formata_p(p_bp)))
  } else {
    cat(sprintf("[\u2713] : O teste indicou homocedasticidade (p = %s),\n    ou seja, a vari\u00e2ncia dos res\u00edduos \u00e9 considerada constante\n    ao longo dos valores ajustados.\n\n",
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

  invisible(list(shapiro = st, bp = list(statistic = bp_stat, p.value = p_bp)))
}
