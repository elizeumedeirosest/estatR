#' @title Regressao de Poisson
#' @description Ajusta um modelo de regressao de Poisson para dados de contagem.
#' Apresenta estatisticas de ajuste, fator de sobredispersao, coeficientes,
#' Razao de Taxa de Incidencia (IRR) e interpretacao automatizada.
#' @param formula Formula do modelo (ex: Y ~ X1 + X2).
#' @param dados Data frame contendo as variaveis.
#' @param grafico Logico. Se TRUE, plota o grafico de Observado vs Predito.
#' @export
regressao_poisson <- function(formula, dados, grafico = TRUE) {
  
  # 1. Validacao e Ajuste
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  modelo <- glm(formula, data = dados, family = poisson(link = "log"))
  
  y_name <- as.character(formula[[2]])
  x_names <- attr(terms(formula), "term.labels")
  if (length(x_names) == 0) stop("O modelo deve ter pelo menos um preditor.")
  
  # 2. Extrair Resultados
  resumo <- summary(modelo)
  coefs <- resumo$coefficients
  
  # Intervalo de confianca e IRR (Incidence Rate Ratio = exp(beta))
  suppressMessages(ci <- confint(modelo))
  if (is.null(dim(ci))) {
    ci <- matrix(ci, ncol = 2, dimnames = list(names(coefs)[1], c("2.5 %", "97.5 %")))
  }
  
  # Deviance e LRT
  dev_nula <- modelo$null.deviance
  df_nula  <- modelo$df.null
  dev_res  <- modelo$deviance
  df_res   <- modelo$df.residual
  
  lrt_stat <- dev_nula - dev_res
  lrt_df   <- df_nula - df_res
  p_lrt    <- pchisq(lrt_stat, lrt_df, lower.tail = FALSE)
  
  # Pseudo-R2 McFadden
  loglik_mod  <- as.numeric(logLik(modelo))
  mod_nulo    <- glm(reformulate("1", y_name), data = dados, family = poisson(link = "log"))
  loglik_nulo <- as.numeric(logLik(mod_nulo))
  pseudo_r2   <- 1 - (loglik_mod / loglik_nulo)
  
  # Sobredispersao (Pearson Chi-2 / df)
  pearson_res <- residuals(modelo, type = "pearson")
  pearson_chi2 <- sum(pearson_res^2)
  dispersion_ratio <- pearson_chi2 / df_res
  
  # 3. Helpers de Formatacao
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
  .fmt_irr <- function(v) {
    if (is.na(v) || is.infinite(v)) return("-")
    if (abs(v) >= 1e4 || (abs(v) < 0.001 && v != 0)) return(formatC(v, format = "e", digits = 2))
    sprintf("%.3f", v)
  }
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  
  # 4. Impressao Console
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("REGRESS\u00c3O DE POISSON")
  } else {
    cat("\n\u2500\u2500 REGRESS\u00c3O DE POISSON \u2500\u2500\n")
  }
  
  # DIAGNOSTICO DO AJUSTE
  if (exists(".print_topico", mode = "function")) .print_topico("DIAGN\u00d3STICO DO AJUSTE")
  
  w_diag <- c(met = 38, val = 15)
  sep_diag <- paste0("  ", strrep("\u2500", sum(w_diag) + 3))
  cat(sep_diag, "\n")
  cat("  ", .pad("M\u00e9trica", w_diag["met"], "left"), .pad("Valor", w_diag["val"], "center"), "\n", sep = "")
  cat(sep_diag, "\n")
  cat("  ", .pad("Deviance Nula", w_diag["met"], "left"), .pad(sprintf("%.2f (df=%d)", dev_nula, df_nula), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Deviance Residual", w_diag["met"], "left"), .pad(sprintf("%.2f (df=%d)", dev_res, df_res), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Sobredispers\u00e3o (\u03c7\u00b2/df)", w_diag["met"], "left"), .pad(sprintf("%.3f", dispersion_ratio), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Teste LRT (\u03c7\u00b2)", w_diag["met"], "left"), .pad(sprintf("%.2f (%s)", lrt_stat, .fmt_p(p_lrt)), w_diag["val"], "center"), .ast(p_lrt), "\n", sep = "")
  cat("  ", .pad("Pseudo-R\u00b2 (McFadden)", w_diag["met"], "left"), .pad(sprintf("%.1f%%", pseudo_r2 * 100), w_diag["val"], "center"), "\n", sep = "")
  cat(sep_diag, "\n\n")
  
  # COEFICIENTES E IRR
  if (exists(".print_topico", mode = "function")) .print_topico("COEFICIENTES E RAZ\u00c3O DE TAXA DE INCID\u00caNCIA (IRR)")
  
  w_c <- c(var=14, est=12, err=11, z=9, pval=14, irr=10, ic_inf=12, ic_sup=12)
  sep_coef <- paste0("  ", strrep("\u2500", sum(w_c) + length(w_c) - 1))
  
  cat(sep_coef, "\n")
  cat("  ",
      .pad("Vari\u00e1vel", w_c["var"], "left"),
      .pad("Estimativa", w_c["est"], "center"),
      .pad("Erro Pad.", w_c["err"], "center"),
      .pad("Z", w_c["z"], "center"),
      .pad("p-valor", w_c["pval"], "center"),
      .pad("IRR", w_c["irr"], "center"),
      .pad("IC(IRR) Inf", w_c["ic_inf"], "center"),
      .pad("IC(IRR) Sup", w_c["ic_sup"], "center"),
      "\n", sep = "")
  cat(sep_coef, "\n")
  
  for (i in seq_len(nrow(coefs))) {
    nome_var <- rownames(coefs)[i]
    if (nome_var == "(Intercept)") nome_var <- "(Intercepto)"
    
    est  <- coefs[i, 1]
    err  <- coefs[i, 2]
    zval <- coefs[i, 3]
    pval <- coefs[i, 4]
    
    irr_val <- exp(est)
    
    cat("  ",
        .pad(nome_var,                           w_c["var"],    "left"),
        .pad(sprintf("%.4f", est),               w_c["est"],    "center"),
        .pad(sprintf("%.4f", err),               w_c["err"],    "center"),
        .pad(sprintf("%.2f", zval),              w_c["z"],      "center"),
        .pad(paste(.fmt_p(pval), .ast(pval)),    w_c["pval"],   "center"),
        .pad(.fmt_irr(irr_val),                  w_c["irr"],    "center"),
        .pad(.fmt_irr(exp(ci[i, 1])),            w_c["ic_inf"], "center"),
        .pad(.fmt_irr(exp(ci[i, 2])),            w_c["ic_sup"], "center"),
        "\n", sep = "")
  }
  cat(sep_coef, "\n")
  cat("  Notas:\n")
  cat("  Signific\u00e2ncia: *** p < 0,01   ** p < 0,05   * p < 0,10\n")
  cat("  IRR > 1: aumenta a contagem esperada. IRR < 1: reduz a contagem esperada.\n\n")
  
  # INTERPRETACAO
  if (exists(".print_topico", mode = "function")) .print_topico("INTERPRETA\u00c7\u00c3O")
  
  cat(sprintf("  O modelo de regress\u00e3o de Poisson foi ajustado para\n  estimar a contagem de '%s'.\n\n", y_name))
  
  if (p_lrt < 0.05) {
    cat(sprintf("  O modelo \u00e9 globalmente v\u00e1lido (LRT %s) e apresenta um\n  Pseudo-R\u00b2 de McFadden de %.1f%%, indicando a propor\u00e7\u00e3o da variabilidade\n  da contagem explicada pelo modelo.\n\n", 
                ifelse(p_lrt < 0.001, "p < 0,001", sprintf("p = %s", .fmt_p(p_lrt))), pseudo_r2 * 100))
  } else {
    cat(sprintf("  O modelo n\u00e3o se mostrou globalmente significativo (LRT p = %s).\n  Os preditores n\u00e3o explicam a variabilidade da contagem melhor que um modelo nulo.\n\n", .fmt_p(p_lrt)))
  }
  
  if (dispersion_ratio > 1.2) {
    cat(sprintf("  [!] ATEN\u00c7\u00c3O: O fator de sobredispers\u00e3o (\u03c7\u00b2/df) \u00e9 %.2f (maior que 1.2).\n  Isso indica que a vari\u00e2ncia dos dados \u00e9 maior que a m\u00e9dia.\n  O modelo de Poisson pode subestimar os erros padr\u00e3o.\n  Considere utilizar uma Regress\u00e3o Quase-Poisson ou Binomial Negativa.\n\n", dispersion_ratio))
  } else if (dispersion_ratio < 0.8) {
    cat(sprintf("  [!] ATEN\u00c7\u00c3O: O fator de subdispers\u00e3o (\u03c7\u00b2/df) \u00e9 %.2f (menor que 0.8).\n  A vari\u00e2ncia \u00e9 menor que a m\u00e9dia. Considere uma Regress\u00e3o Quase-Poisson.\n\n", dispersion_ratio))
  } else {
    cat(sprintf("  [\u2713] O fator de dispers\u00e3o (\u03c7\u00b2/df) \u00e9 %.2f, pr\u00f3ximo de 1.\n  N\u00e3o h\u00e1 ind\u00edcios fortes de viola\u00e7\u00e3o do pressuposto (m\u00e9dia = vari\u00e2ncia).\n\n", dispersion_ratio))
  }
  
  # Interpretacao dos coeficientes
  significantes <- 0
  for (i in 2:nrow(coefs)) {
    if (coefs[i, 4] < 0.10) {
      significantes <- significantes + 1
      nome_pred <- rownames(coefs)[i]
      pval_pred <- coefs[i, 4]
      irr_val   <- exp(coefs[i, 1])
      direcao   <- if (irr_val > 1) "aumenta" else "reduz"
      pct_irr   <- abs(irr_val - 1) * 100
      
      cat(sprintf("  O preditor '%s' apresentou efeito estatisticamente\n  significativo (p %s), mantendo os demais preditores constantes.\n  Um aumento de uma unidade em '%s' %s a contagem esperada de '%s'\n  em %.1f%% (IRR = %s; IC 95%%: [%s; %s]).\n\n",
                  nome_pred, ifelse(pval_pred < 0.001, "< 0,001", sprintf("= %s", .fmt_p(pval_pred))),
                  nome_pred, direcao, y_name,
                  pct_irr, .fmt_irr(irr_val),
                  .fmt_irr(exp(ci[i, 1])), .fmt_irr(exp(ci[i, 2]))))
    }
  }
  
  if (significantes == 0) {
    cat("  Nenhum preditor individual apresentou efeito estatisticamente significativo (p < 0,10).\n\n")
  }
  
  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(strrep("\u2500", sum(w_c) + length(w_c) - 1), "\n\n")
  }
  
  # 5. Grafico: Observado vs Predito
  if (grafico && exists("ggplot2", mode = "environment") && requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch({
      df_plot <- data.frame(
        Observado = dados[[y_name]],
        Predito = modelo$fitted.values
      )
      
      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Predito, y = Observado)) +
        ggplot2::geom_point(color = "#555555", alpha = 0.7, size = 2) +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "#D90429", linetype = "dashed", linewidth = 1) +
        ggplot2::labs(
          title = "Ajuste do Modelo de Poisson",
          subtitle = "Valores Observados vs Contagem Esperada (Predita)",
          x = "Contagem Predita",
          y = "Contagem Observada"
        )
        
      if (exists("tema_estatR", mode = "function")) {
        p <- p + tema_estatR(estilo = 2)
      } else {
        p <- p + ggplot2::theme_minimal()
      }
      
      print(p)
    }, error = function(e) {
      message("[Aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico: ", e$message)
    })
  }
  
  invisible(modelo)
}
