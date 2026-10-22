# ==============================================================================
# MÓDULO: REGRESSÃO BINOMIAL NEGATIVA INFLADA DE ZERO (NEGBIN IZ)
# Funções internas — não exportadas
# ==============================================================================

# Função interna principal: ajusta modelo NEGBIN IZ via pscl::zeroinfl
.regressao_negbin_iz <- function(formula, dados, formula_iz = ~1, grafico = TRUE) {

  if (!requireNamespace("pscl", quietly = TRUE)) {
    stop("O pacote 'pscl' é necessário para modelos NEGBIN IZ. Instale com install.packages('pscl').")
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("O pacote 'MASS' é necessário para o Teste de Vuong. Instale com install.packages('MASS').")
  }
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")

  y_name  <- as.character(formula[[2]])
  x_names <- attr(terms(formula), "term.labels")

  # Formula combinada para pscl::zeroinfl (contagem | inflacao)
  f_str <- deparse(formula)
  if (!grepl("\\|", f_str)) {
    iz_rhs <- gsub("^~\\s*", "", deparse(formula_iz))
    formula_usar <- stats::as.formula(paste(f_str, "|", iz_rhs), env = environment(formula))
  } else {
    formula_usar <- formula
  }

  # Ajuste do modelo NEGBIN IZ
  mod_niz <- tryCatch(
    pscl::zeroinfl(formula_usar, data = dados, dist = "negbin", link = "logit"),
    error = function(e) stop("Erro ao ajustar o modelo NEGBIN IZ: ", e$message)
  )

  # Modelo Binomial Negativo simples para o Teste de Vuong
  mod_nb <- tryCatch(
    MASS::glm.nb(formula, data = dados),
    error = function(e) NULL
  )

  # Extração de métricas
  resumo      <- summary(mod_niz)
  coefs_count <- resumo$coefficients$count
  coefs_zero  <- resumo$coefficients$zero
  loglik_val  <- as.numeric(logLik(mod_niz))
  aic_val     <- AIC(mod_niz)
  bic_val     <- BIC(mod_niz)
  n_obs       <- nrow(dados)
  theta_val   <- mod_niz$theta

  # Zeros observados vs esperados
  y_obs     <- dados[[y_name]]
  zeros_obs <- sum(y_obs == 0)
  pct_obs   <- zeros_obs / n_obs * 100

  # Zeros esperados pelo modelo
  zeros_esp <- sum(predict(mod_niz, type = "prob")[, 1])
  pct_esp   <- zeros_esp / n_obs * 100

  # Pseudo-R² de McFadden
  mod_nulo    <- tryCatch(
    glm(reformulate("1", y_name), data = dados, family = poisson()),
    error = function(e) NULL
  )
  loglik_nulo <- if (!is.null(mod_nulo)) as.numeric(logLik(mod_nulo)) else NA
  pseudo_r2   <- if (!is.na(loglik_nulo) && loglik_nulo != 0) 1 - (loglik_val / loglik_nulo) else NA

  # Teste de Vuong
  vuong_stat <- NA_real_
  vuong_p    <- NA_real_
  vuong_txt  <- "Não calculado"
  if (!is.null(mod_nb)) {
    vuong_res <- tryCatch(pscl::vuong(mod_niz, mod_nb), error = function(e) NULL)
    if (!is.null(vuong_res)) {
      vuong_stat <- tryCatch(vuong_res$vuong["Vuong z-statistic"], error = function(e) NA_real_)
      vuong_p    <- tryCatch(vuong_res$p.value["one.sided (H_A: model1 > model2)"], error = function(e) NA_real_)
      if (!is.na(vuong_stat) && !is.na(vuong_p)) {
        vuong_txt <- if (vuong_p < 0.05) "NEGBIN IZ preferível à BN" else "Equivalente à BN"
      }
    }
  }

  # 3. Helpers de Formatação
  .fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }
  .ast <- function(p) {
    if (is.na(p) || p >= 0.10) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    return("*")
  }
  .fmt_num <- function(v, digits = 3) {
    if (is.na(v) || is.infinite(v)) return("-")
    formatC(v, format = "f", digits = digits, decimal.mark = ",")
  }
  .fmt_irr <- function(v) {
    if (is.na(v) || is.infinite(v)) return("-")
    if (abs(v) >= 1e4 || (abs(v) < 0.001 && v != 0)) return(formatC(v, format = "e", digits = 2))
    sprintf("%.3f", v)
  }
  .pad <- function(s, w, align = "center") {
    s   <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left")  return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }

  # 4. Impressão Console
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("REGRESSÃO BINOMIAL NEGATIVA INFLADA DE ZERO (NEGBIN IZ)")
  } else {
    cat("\n── REGRESSÃO BINOMIAL NEGATIVA INFLADA DE ZERO (NEGBIN IZ) ──\n")
  }

  cat(sprintf("  Fórmula (Contagem):     %s\n", deparse(formula)))
  cat(sprintf("  Fórmula (Inflação):     %s\n", deparse(formula_iz)))
  cat(sprintf("  Observações:            %d\n\n", n_obs))

  # DIAGNÓSTICO DO AJUSTE
  if (exists(".print_topico", mode = "function")) .print_topico("DIAGNÓSTICO DO AJUSTE")

  w_diag   <- c(met = 38, val = 16)
  sep_diag <- paste0("  ", strrep("─", sum(w_diag) + 3))
  cat(sep_diag, "\n")
  cat("  ", .pad("Métrica", w_diag["met"], "left"), .pad("Valor", w_diag["val"], "center"), "\n", sep = "")
  cat(sep_diag, "\n")
  cat("  ", .pad("Log-Verossimilhança",       w_diag["met"], "left"), .pad(.fmt_num(loglik_val, 2), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("AIC",                        w_diag["met"], "left"), .pad(.fmt_num(aic_val, 2),   w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("BIC",                        w_diag["met"], "left"), .pad(.fmt_num(bic_val, 2),   w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Parâmetro \u03b8 (theta)",   w_diag["met"], "left"), .pad(.fmt_num(theta_val, 3), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Zeros Observados",           w_diag["met"], "left"),
      .pad(sprintf("%d (%.1f%%)", zeros_obs, pct_obs), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Zeros Esperados (Modelo)",   w_diag["met"], "left"),
      .pad(sprintf("%.1f (%.1f%%)", zeros_esp, pct_esp), w_diag["val"], "center"), "\n", sep = "")
  cat(sep_diag, "\n\n")

  # COEFICIENTES — COMPONENTE DE CONTAGEM
  if (exists(".print_topico", mode = "function")) .print_topico("COEFICIENTES — COMPONENTE DE CONTAGEM (Binomial Negativa / Log)")

  w_c      <- c(var = 16, est = 13, err = 12, z = 13, pval = 16, irr = 8)
  sep_coef <- paste0("  ", strrep("─", sum(w_c) + length(w_c) - 1))
  cat(sep_coef, "\n")
  cat("  ",
      .pad("Variável",    w_c["var"], "left"),
      .pad("Estimativa",  w_c["est"], "center"),
      .pad("Erro Pad.",   w_c["err"], "center"),
      .pad("Estatística", w_c["z"],   "center"),
      .pad("p-valor",     w_c["pval"],"center"),
      .pad("IRR",         w_c["irr"], "center"),
      "\n", sep = "")
  cat(sep_coef, "\n")

  for (i in seq_len(nrow(coefs_count))) {
    nv   <- rownames(coefs_count)[i]
    if (nv == "(Intercept)") nv <- "(Intercepto)"
    est  <- coefs_count[i, 1]
    err  <- coefs_count[i, 2]
    zval <- coefs_count[i, 3]
    pval <- coefs_count[i, 4]
    cat("  ",
        .pad(nv,                              w_c["var"],  "left"),
        .pad(.fmt_num(est, 3),                w_c["est"],  "center"),
        .pad(.fmt_num(err, 3),                w_c["err"],  "center"),
        .pad(.fmt_num(zval, 2),               w_c["z"],    "center"),
        .pad(paste(.fmt_p(pval), .ast(pval)), w_c["pval"], "center"),
        .pad(.fmt_irr(exp(est)),              w_c["irr"],  "center"),
        "\n", sep = "")
  }
  cat(sep_coef, "\n\n")

  # COEFICIENTES — COMPONENTE DE INFLAÇÃO DE ZERO
  if (exists(".print_topico", mode = "function")) .print_topico("COEFICIENTES — COMPONENTE DE INFLAÇÃO DE ZERO (Binomial / Logit)")

  cat(sep_coef, "\n")
  cat("  ",
      .pad("Variável",    w_c["var"], "left"),
      .pad("Estimativa",  w_c["est"], "center"),
      .pad("Erro Pad.",   w_c["err"], "center"),
      .pad("Estatística", w_c["z"],   "center"),
      .pad("p-valor",     w_c["pval"],"center"),
      .pad("OR",          w_c["irr"], "center"),
      "\n", sep = "")
  cat(sep_coef, "\n")

  for (i in seq_len(nrow(coefs_zero))) {
    nv   <- rownames(coefs_zero)[i]
    if (nv == "(Intercept)") nv <- "(Intercepto)"
    est  <- coefs_zero[i, 1]
    err  <- coefs_zero[i, 2]
    zval <- coefs_zero[i, 3]
    pval <- coefs_zero[i, 4]
    cat("  ",
        .pad(nv,                              w_c["var"],  "left"),
        .pad(.fmt_num(est, 3),                w_c["est"],  "center"),
        .pad(.fmt_num(err, 3),                w_c["err"],  "center"),
        .pad(.fmt_num(zval, 2),               w_c["z"],    "center"),
        .pad(paste(.fmt_p(pval), .ast(pval)), w_c["pval"], "center"),
        .pad(.fmt_irr(exp(est)),              w_c["irr"],  "center"),
        "\n", sep = "")
  }
  cat(sep_coef, "\n")
  cat("  Notas: IRR = exp(\u03b2) \u2014 Razão de Taxa de Incidência (componente contagem)\n")
  cat("         OR  = exp(\u03b2) \u2014 Razão de Chances de ser zero estrutural (inflação)\n")
  cat("         \u03b8  = parâmetro de dispersão da Binomial Negativa\n")
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n\n")

  # TESTE DE VUONG
  if (exists(".print_topico", mode = "function")) .print_topico("TESTE DE VUONG (NEGBIN IZ vs Binomial Negativa)")

  w_v   <- c(met = 38, val = 16)
  sep_v <- paste0("  ", strrep("─", sum(w_v) + 3))
  cat(sep_v, "\n")
  cat("  ", .pad("Métrica",      w_v["met"], "left"), .pad("Valor",    w_v["val"], "center"), "\n", sep = "")
  cat(sep_v, "\n")
  cat("  ", .pad("Estatística V", w_v["met"], "left"), .pad(.fmt_num(vuong_stat, 2), w_v["val"], "center"), "\n", sep = "")
  cat("  ", .pad("p-valor",       w_v["met"], "left"), .pad(.fmt_p(vuong_p),         w_v["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Conclusão",     w_v["met"], "left"), .pad(vuong_txt,               w_v["val"], "center"), "\n", sep = "")
  cat(sep_v, "\n")
  cat("  * H\u2080: os dois modelos têm desempenho equivalente\n\n")

  # INTERPRETAÇÃO
  if (exists(".print_topico", mode = "function")) .print_topico("INTERPRETAÇÃO")

  cat(sprintf("  O modelo NEGBIN IZ foi ajustado para explicar %s a partir de %s.\n",
              y_name, paste(x_names, collapse = ", ")))
  cat(sprintf("  O parâmetro de dispersão \u03b8 = %s controla a sobredispersão da\n", .fmt_num(theta_val, 3)))
  cat("  componente de contagem (variância = \u03bc + \u03bc²/\u03b8).\n")
  cat(sprintf("  O modelo apresenta %.0f zeros estruturais (parte de inflação) e\n", zeros_esp))
  cat(sprintf("  %.0f zeros Binomial Negativa (parte de contagem).\n", max(0, zeros_obs - zeros_esp)))
  if (!is.na(vuong_stat) && !is.na(vuong_p)) {
    pref <- if (vuong_p < 0.05) "preferível" else "equivalente"
    cat(sprintf("  O Teste de Vuong indica que o modelo NEGBIN IZ é %s\n", pref))
    cat(sprintf("  à Binomial Negativa simples (V = %s, p = %s).\n\n",
                .fmt_num(vuong_stat, 2), .fmt_p(vuong_p)))
  } else {
    cat("  O Teste de Vuong não pôde ser calculado.\n\n")
  }

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(strrep("─", 78), "\n\n")
  }

  # 5. Gráfico: frequências observadas vs esperadas
  if (grafico) {
    tryCatch({
      max_y    <- max(y_obs)
      freq_obs <- tabulate(y_obs + 1L, nbins = max_y + 1)
      freq_obs <- freq_obs / n_obs
      prob_mat <- predict(mod_niz, type = "prob")
      freq_esp <- colMeans(prob_mat)[seq_len(max_y + 1)]
      contagens <- 0:max_y

      df_plot <- data.frame(
        Contagem  = contagens,
        Observado = freq_obs,
        Esperado  = freq_esp,
        Destaque  = ifelse(contagens == 0, "zero", "positivo")
      )

      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Contagem)) +
        ggplot2::geom_col(ggplot2::aes(y = Observado, fill = Destaque),
                          alpha = 0.75, width = 0.6) +
        ggplot2::geom_point(ggplot2::aes(y = Esperado), color = "#D90429", size = 3) +
        ggplot2::geom_line(ggplot2::aes(y = Esperado), color = "#D90429",
                           linewidth = 0.8, linetype = "dashed") +
        ggplot2::scale_fill_manual(values = c(zero = "#D90429", positivo = "#555555"),
                                   name = "Categoria") +
        ggplot2::labs(
          title    = "Frequências Observadas vs Esperadas (NEGBIN IZ)",
          subtitle = "Barras = Observado | Ponto/Linha = Esperado pelo Modelo",
          x = "Contagem", y = "Proporção"
        )

      if (exists("meu_tema", mode = "function")) {
        p <- p + meu_tema(estilo = 1, grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        p <- p + tema_estatR(estilo = 1)
      } else {
        p <- p + ggplot2::theme_minimal()
      }

      print(p)
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o gráfico NEGBIN IZ: ", e$message)
    })
  }

  invisible(mod_niz)
}


# Métricas do modelo NEGBIN IZ
.metricas_negbin_iz <- function(modelo) {

  if (!inherits(modelo, "zeroinfl") || modelo$dist != "negbin") {
    stop("O modelo fornecido não é um modelo NEGBIN IZ (zeroinfl com dist='negbin').")
  }

  loglik_val <- as.numeric(logLik(modelo))
  aic_val    <- AIC(modelo)
  bic_val    <- BIC(modelo)
  k          <- length(coef(modelo)) + 1  # +1 para o theta
  n          <- nrow(modelo$model)
  aicc_val   <- aic_val + (2 * k * (k + 1)) / max(1, n - k - 1)
  theta_val  <- modelo$theta

  # Pseudo-R² McFadden (referência: Poisson nulo)
  y_name    <- as.character(formula(modelo)[[2]])
  mod_nulo  <- tryCatch(
    glm(reformulate("1", y_name), data = modelo$model, family = poisson()),
    error = function(e) NULL
  )
  loglik_nulo <- if (!is.null(mod_nulo)) as.numeric(logLik(mod_nulo)) else NA
  pseudo_r2   <- if (!is.na(loglik_nulo) && loglik_nulo != 0) 1 - (loglik_val / loglik_nulo) else NA

  .fmt_num <- function(v, digits = 2) {
    if (is.na(v)) return("-")
    formatC(v, format = "f", digits = digits, decimal.mark = ",")
  }
  .pad <- function(s, w, align = "center") {
    s   <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left")  return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }

  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("MÉTRICAS — MODELO NEGBIN IZ")
  } else {
    cat("\n── MÉTRICAS — MODELO NEGBIN IZ ──\n")
  }

  w_m  <- c(met = 34, val = 14)
  lin  <- paste0("  ", strrep("─", sum(w_m) + 1))

  cat(lin, "\n")
  cat("  ", .pad("Métrica", w_m["met"], "left"), .pad("Valor", w_m["val"], "center"), "\n", sep = "")
  cat(lin, "\n")
  cat("  ", .pad("Log-Verossimilhança",   w_m["met"], "left"), .pad(.fmt_num(loglik_val, 2),  w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("AIC",                   w_m["met"], "left"), .pad(.fmt_num(aic_val, 2),     w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("AICc",                  w_m["met"], "left"), .pad(.fmt_num(aicc_val, 2),    w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("BIC",                   w_m["met"], "left"), .pad(.fmt_num(bic_val, 2),     w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Parâmetro \u03b8 (theta)", w_m["met"], "left"), .pad(.fmt_num(theta_val, 3), w_m["val"], "center"), "\n", sep = "")
  if (!is.na(pseudo_r2)) {
    cat("  ", .pad("Pseudo-R² (McFadden)", w_m["met"], "left"), .pad(.fmt_num(pseudo_r2, 4), w_m["val"], "center"), "\n", sep = "")
  }
  cat(lin, "\n\n")

  if (exists(".print_rodape", mode = "function")) .print_rodape()

  invisible(list(loglik = loglik_val, aic = aic_val, aicc = aicc_val,
                 bic = bic_val, theta = theta_val, pseudo_r2 = pseudo_r2))
}


# Análise de resíduos do modelo NEGBIN IZ
.analise_residual_negbin_iz <- function(modelo, grafico = TRUE) {

  if (!inherits(modelo, "zeroinfl") || modelo$dist != "negbin") {
    stop("O modelo fornecido não é um modelo NEGBIN IZ válido.")
  }

  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("ANÁLISE DE RESÍDUOS — REGRESSÃO NEGBIN IZ")
  } else {
    cat("\n── ANÁLISE DE RESÍDUOS — REGRESSÃO NEGBIN IZ ──\n")
  }

  # Resíduos de Pearson
  pearson_res <- residuals(modelo, type = "pearson")
  y_fit       <- fitted(modelo)
  n           <- length(pearson_res)

  cor_ok   <- "\033[1;32m"
  cor_erro <- "\033[1;31m"
  cor_warn <- "\033[1;33m"
  reset    <- "\033[0m"

  # Média zero (t-test)
  if (exists(".print_topico", mode = "function")) .print_topico("MÉDIA ZERO DOS RESÍDUOS DE PEARSON")
  tt <- t.test(pearson_res, mu = 0)
  cat(sprintf("  t = %.3f, p = %s\n", tt$statistic,
              if (tt$p.value < 0.001) "< 0,001" else formatC(tt$p.value, format = "f", digits = 3, decimal.mark = ",")))
  if (tt$p.value > 0.05) {
    cat(sprintf("  %s[\u2713]%s Média dos resíduos de Pearson não difere significativamente de zero.\n\n", cor_ok, reset))
  } else {
    cat(sprintf("  %s[!]%s Média dos resíduos de Pearson difere de zero (p < 0,05).\n\n", cor_erro, reset))
  }

  # Normalidade
  if (exists(".print_topico", mode = "function")) .print_topico("NORMALIDADE DOS RESÍDUOS DE PEARSON")
  if (n <= 5000) {
    sw <- shapiro.test(pearson_res)
    cat(sprintf("  Shapiro-Wilk: W = %.4f, p = %s\n", sw$statistic,
                if (sw$p.value < 0.001) "< 0,001" else formatC(sw$p.value, format = "f", digits = 3, decimal.mark = ",")))
    if (sw$p.value > 0.05) {
      cat(sprintf("  %s[\u2713]%s Resíduos de Pearson não rejeitam normalidade.\n\n", cor_ok, reset))
    } else {
      cat(sprintf("  %s[!]%s Resíduos de Pearson não seguem distribuição normal (esperado para modelos de contagem).\n\n", cor_warn, reset))
    }
  } else {
    cat("  (n > 5000: teste de normalidade omitido)\n\n")
  }

  # Autocorrelação (Ljung-Box, lag = 10)
  if (exists(".print_topico", mode = "function")) .print_topico("AUTOCORRELAÇÃO (LJUNG-BOX, lag = 10)")
  lb <- Box.test(pearson_res, lag = 10, type = "Ljung-Box")
  cat(sprintf("  Estatística Q = %.3f, p = %s\n", lb$statistic,
              if (lb$p.value < 0.001) "< 0,001" else formatC(lb$p.value, format = "f", digits = 3, decimal.mark = ",")))
  if (lb$p.value > 0.05) {
    cat(sprintf("  %s[\u2713]%s Sem evidência de autocorrelação nos resíduos.\n\n", cor_ok, reset))
  } else {
    cat(sprintf("  %s[!]%s Autocorrelação detectada nos resíduos (p < 0,05).\n\n", cor_erro, reset))
  }

  if (exists(".print_rodape", mode = "function")) .print_rodape()

  # Gráfico
  if (grafico) {
    tryCatch({
      df_plot <- data.frame(
        Ajustados = y_fit,
        Residuos  = pearson_res
      )

      tema_base <- if (exists("meu_tema", mode = "function")) {
        meu_tema(grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2)
      } else {
        ggplot2::theme_minimal()
      }

      p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Ajustados, y = Residuos)) +
        ggplot2::geom_point(color = "#555555", alpha = 0.6) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#D90429") +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429",
                             linewidth = 0.8, formula = y ~ x) +
        ggplot2::labs(title = "Resíduos de Pearson vs Valores Preditos",
                      x = "Valores Preditos", y = "Resíduos de Pearson") +
        tema_base

      p2 <- ggplot2::ggplot(df_plot, ggplot2::aes(sample = Residuos)) +
        ggplot2::stat_qq(color = "#555555", alpha = 0.6) +
        ggplot2::stat_qq_line(color = "#D90429", linewidth = 0.8) +
        ggplot2::labs(title = "QQ-Plot dos Resíduos de Pearson",
                      x = "Quantis Teóricos", y = "Quantis Amostrais") +
        tema_base

      if (requireNamespace("patchwork", quietly = TRUE)) {
        print(p1 | p2)
      } else {
        print(p1); print(p2)
      }
    }, error = function(e) {
      message("[Aviso] Erro ao gerar painel de resíduos NEGBIN IZ: ", e$message)
    })
  }

  invisible(list(t_test = tt, ljung_box = lb))
}


# Seleção de modelos NEGBIN IZ (Best Subsets por AIC)
.selecao_modelos_negbin_iz <- function(formula, dados, top = 5) {

  if (!requireNamespace("pscl", quietly = TRUE)) {
    stop("O pacote 'pscl' é necessário. Instale com install.packages('pscl').")
  }

  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("SELEÇÃO DE MODELOS — NEGBIN IZ (BINOMIAL NEGATIVA INFLADA DE ZERO)")
  } else {
    cat("\n── SELEÇÃO DE MODELOS — NEGBIN IZ ──\n")
  }

  y_name     <- as.character(formula[[2]])
  preditores <- attr(terms(formula, data = dados), "term.labels")

  if (length(preditores) > 10) {
    cat(sprintf("  Muitos preditores (%d). Usando seleção heurística (top 10)...\n\n", length(preditores)))
    p_vals <- sapply(preditores, function(p) {
      f   <- as.formula(paste(y_name, "~", p))
      mod <- tryCatch(pscl::zeroinfl(f, data = dados, dist = "negbin"), error = function(e) NULL)
      if (is.null(mod)) return(1.0)
      sm  <- summary(mod)$coefficients$count
      if (nrow(sm) > 1) sm[2, 4] else 1.0
    })
    preditores <- names(sort(p_vals)[1:10])
  }

  cat(sprintf("  Variável Resposta: %s\n", y_name))
  cat(sprintf("  Preditores base:   %s\n", paste(preditores, collapse = ", ")))

  combinacoes <- list()
  for (i in seq_along(preditores)) {
    combinacoes <- c(combinacoes, combn(preditores, i, simplify = FALSE))
  }
  cat(sprintf("  Método utilizado:  Exaustivo (%d modelos avaliados)\n\n", length(combinacoes)))

  resultados <- data.frame(Modelo = character(), k = integer(), LogLik = numeric(),
                           AIC = numeric(), BIC = numeric(), Delta_AIC = numeric(),
                           stringsAsFactors = FALSE)

  for (vars in combinacoes) {
    f_str <- paste(y_name, "~", paste(vars, collapse = " + "))
    mod   <- suppressWarnings(tryCatch(
      pscl::zeroinfl(as.formula(f_str), data = dados, dist = "negbin"),
      error = function(e) NULL
    ))
    if (!is.null(mod)) {
      resultados <- rbind(resultados, data.frame(
        Modelo    = paste(vars, collapse = " + "),
        k         = length(vars) + 1L,  # +1 para o theta
        LogLik    = as.numeric(logLik(mod)),
        AIC       = AIC(mod),
        BIC       = BIC(mod),
        Delta_AIC = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
  }

  resultados       <- resultados[order(resultados$AIC), ]
  resultados$Delta_AIC <- resultados$AIC - min(resultados$AIC)
  df_top           <- head(resultados, top)
  df_top$Ranking   <- paste0(seq_len(nrow(df_top)), "º")
  df_top           <- df_top[, c("Ranking", "Modelo", "k", "AIC", "BIC", "LogLik", "Delta_AIC")]

  df_fmt           <- df_top
  df_fmt$AIC       <- sprintf("%.1f", df_top$AIC)
  df_fmt$BIC       <- sprintf("%.1f", df_top$BIC)
  df_fmt$LogLik    <- sprintf("%.1f", df_top$LogLik)
  df_fmt$Delta_AIC <- sprintf("%.2f", df_top$Delta_AIC)
  names(df_fmt)[names(df_fmt) == "Delta_AIC"] <- "\u0394AIC"

  if (exists(".print_topico", mode = "function")) {
    .print_topico(sprintf("TOP %d MODELOS (Ordenados por AIC)", nrow(df_fmt)))
  }

  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_fmt, align = c("left", "left", "center", "center", "center", "center", "center"))
  } else {
    print(df_fmt, row.names = FALSE)
  }

  cat("  * k = Número de parâmetros preditores (+1 do Theta).\n")
  cat("  * \u0394AIC = diferença em relação ao melhor modelo.\n\n")
  if (exists(".print_rodape", mode = "function")) .print_rodape()

  invisible(df_top)
}
