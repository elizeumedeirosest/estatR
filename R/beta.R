#' @title Regressao Beta
#' @description Ajusta um modelo de regressao Beta para variaveis respostas
#' que sao proporcoes ou taxas estritamente no intervalo (0, 1).
#' Apresenta estatisticas de ajuste, parametro de precisao (phi), coeficientes,
#' Odds Ratio (OR) para a media e interpretacao automatizada.
#' @param formula Formula do modelo (ex: Y ~ X1 + X2).
#' @param dados Data frame contendo as variaveis.
#' @param grafico Logico. Se TRUE, plota o grafico de Observado vs Predito.
#' @export
regressao_beta <- function(formula, dados, grafico = TRUE) {
  
  if (!requireNamespace("betareg", quietly = TRUE)) {
    stop("O pacote 'betareg' e necessario para ajustar modelos beta. Instale com install.packages('betareg').")
  }
  
  # 1. Validacao e Ajuste
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  
  y_name <- as.character(formula[[2]])
  y_val <- dados[[y_name]]
  
  if (any(y_val <= 0 | y_val >= 1, na.rm = TRUE)) {
    stop(sprintf("A regressao beta exige que a variavel resposta '%s' esteja ESTRITAMENTE no intervalo (0, 1).\nSeus dados contem valores 0 ou 1. Para corrigir, voce pode aplicar a transformacao de Smithson e Verkuilen:\ny_transf = (y * (n - 1) + 0.5) / n", y_name))
  }
  
  # Ajuste do modelo (link logit padrao)
  modelo <- tryCatch(
    betareg::betareg(formula, data = dados, link = "logit"),
    error = function(e) stop("Erro ao ajustar o modelo Beta: ", e$message)
  )
  
  x_names <- attr(terms(formula), "term.labels")
  if (length(x_names) == 0) stop("O modelo deve ter pelo menos um preditor.")
  
  # 2. Extrair Resultados
  resumo <- summary(modelo)
  coefs_mean <- resumo$coefficients$mean
  coefs_prec <- resumo$coefficients$precision
  
  phi_est <- coefs_prec[1, "Estimate"]
  phi_se  <- coefs_prec[1, "Std. Error"]
  
  # Intervalo de confianca da media
  suppressMessages(ci <- confint(modelo))
  # O confint retorna para media e precisao. Extrair apenas media:
  ci_mean <- ci[1:nrow(coefs_mean), , drop = FALSE]
  
  # Pseudo-R2
  pseudo_r2 <- resumo$pseudo.r.squared
  if (is.na(pseudo_r2)) pseudo_r2 <- 0 # Tratamento pra caso não consiga calcular
  
  # LRT
  loglik_mod <- as.numeric(logLik(modelo))
  mod_nulo <- suppressWarnings(betareg::betareg(stats::reformulate("1", y_name), data = dados, link = "logit"))
  loglik_nulo <- as.numeric(logLik(mod_nulo))
  
  lrt_stat <- 2 * (loglik_mod - loglik_nulo)
  lrt_df <- length(coef(modelo)) - length(coef(mod_nulo))
  if (lrt_df > 0) {
    p_lrt <- pchisq(lrt_stat, lrt_df, lower.tail = FALSE)
  } else {
    p_lrt <- NA
  }
  
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
  .fmt_or <- function(v) {
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
    .print_titulo("REGRESSÃO BETA (Proporções e Taxas)")
  } else {
    cat("\n── REGRESSÃO BETA ──\n")
  }
  
  # DIAGNOSTICO DO AJUSTE
  if (exists(".print_topico", mode = "function")) .print_topico("DIAGNÓSTICO DO AJUSTE")
  
  w_diag <- c(met = 38, val = 15)
  sep_diag <- paste0("  ", strrep("─", sum(w_diag) + 3))
  cat(sep_diag, "\n")
  cat("  ", .pad("Métrica", w_diag["met"], "left"), .pad("Valor", w_diag["val"], "center"), "\n", sep = "")
  cat(sep_diag, "\n")
  cat("  ", .pad("Log-Verossimilhança", w_diag["met"], "left"), .pad(sprintf("%.2f", loglik_mod), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Parâmetro de Precisão (Phi φ)", w_diag["met"], "left"), .pad(sprintf("%.2f", phi_est), w_diag["val"], "center"), "\n", sep = "")
  if (!is.na(p_lrt)) {
    cat("  ", .pad("Teste LRT (χ²)", w_diag["met"], "left"), .pad(sprintf("%.2f (%s)", lrt_stat, .fmt_p(p_lrt)), w_diag["val"], "center"), .ast(p_lrt), "\n", sep = "")
  }
  cat("  ", .pad("Pseudo-R²", w_diag["met"], "left"), .pad(sprintf("%.1f%%", pseudo_r2 * 100), w_diag["val"], "center"), "\n", sep = "")
  cat(sep_diag, "\n\n")
  
  # COEFICIENTES
  if (exists(".print_topico", mode = "function")) .print_topico("COEFICIENTES (Modelo da Média - Link Logit)")
  
  w_c <- c(var=14, est=12, err=11, z=9, pval=14, or=10, ic_inf=12, ic_sup=12)
  sep_coef <- paste0("  ", strrep("─", sum(w_c) + length(w_c) - 1))
  
  cat(sep_coef, "\n")
  cat("  ",
      .pad("Variável", w_c["var"], "left"),
      .pad("Estimativa", w_c["est"], "center"),
      .pad("Erro Pad.", w_c["err"], "center"),
      .pad("Z", w_c["z"], "center"),
      .pad("p-valor", w_c["pval"], "center"),
      .pad("OR*", w_c["or"], "center"),
      .pad("IC(OR) Inf", w_c["ic_inf"], "center"),
      .pad("IC(OR) Sup", w_c["ic_sup"], "center"),
      "\n", sep = "")
  cat(sep_coef, "\n")
  
  for (i in seq_len(nrow(coefs_mean))) {
    nome_var <- rownames(coefs_mean)[i]
    if (nome_var == "(Intercept)") nome_var <- "(Intercepto)"
    
    est  <- coefs_mean[i, 1]
    err  <- coefs_mean[i, 2]
    zval <- coefs_mean[i, 3]
    pval <- coefs_mean[i, 4]
    
    or_val <- exp(est)
    
    cat("  ",
        .pad(nome_var,                           w_c["var"],    "left"),
        .pad(sprintf("%.4f", est),               w_c["est"],    "center"),
        .pad(sprintf("%.4f", err),               w_c["err"],    "center"),
        .pad(sprintf("%.2f", zval),              w_c["z"],      "center"),
        .pad(paste(.fmt_p(pval), .ast(pval)),    w_c["pval"],   "center"),
        .pad(.fmt_or(or_val),                    w_c["or"],     "center"),
        .pad(.fmt_or(exp(ci_mean[i, 1])),        w_c["ic_inf"], "center"),
        .pad(.fmt_or(exp(ci_mean[i, 2])),        w_c["ic_sup"], "center"),
        "\n", sep = "")
  }
  cat(sep_coef, "\n")
  cat("  Notas:\n")
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n")
  cat("  *OR: Odds Ratio (Razão de Chances) da proporção esperada (exp(beta)).\n\n")
  
  # INTERPRETACAO
  if (exists(".print_topico", mode = "function")) .print_topico("INTERPRETAÇÃO")
  
  cat(sprintf("  O modelo de Regressão Beta foi ajustado para estimar a proporção\n  (ou taxa) da variável '%s' (limitada entre 0 e 1).\n\n", y_name))
  
  if (!is.na(p_lrt) && p_lrt < 0.05) {
    cat(sprintf("  O modelo é globalmente válido (LRT %s) e apresenta um\n  Pseudo-R² de %.1f%%, indicando a qualidade do ajustamento.\n\n", 
                ifelse(p_lrt < 0.001, "p < 0,001", sprintf("p = %s", .fmt_p(p_lrt))), pseudo_r2 * 100))
  }
  
  # Interpretacao dos coeficientes
  significantes <- 0
  for (i in 2:nrow(coefs_mean)) {
    if (coefs_mean[i, 4] < 0.10) {
      significantes <- significantes + 1
      nome_pred <- rownames(coefs_mean)[i]
      pval_pred <- coefs_mean[i, 4]
      or_val    <- exp(coefs_mean[i, 1])
      direcao   <- if (or_val > 1) "aumenta" else "reduz"
      
      cat(sprintf("  O preditor '%s' apresentou efeito estatisticamente\n  significativo (p %s). Um aumento de uma unidade em '%s'\n  %s as chances associadas à proporção média de '%s'\n  por um fator de %.2f (OR = %.2f; IC 95%%: [%.2f; %.2f]).\n\n",
                  nome_pred, ifelse(pval_pred < 0.001, "< 0,001", sprintf("= %s", .fmt_p(pval_pred))),
                  nome_pred, direcao, y_name,
                  or_val, or_val, exp(ci_mean[i, 1]), exp(ci_mean[i, 2])))
    }
  }
  
  if (significantes == 0) {
    cat("  Nenhum preditor individual apresentou efeito estatisticamente significativo (p < 0,10).\n\n")
  }
  
  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(strrep("─", sum(w_c) + length(w_c) - 1), "\n\n")
  }
  
  # 5. Grafico: Observado vs Predito
  if (grafico) {
    tryCatch({
      df_plot <- data.frame(
        Observado = y_val,
        Predito = predict(modelo, type = "response")
      )
      
      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Predito, y = Observado)) +
        ggplot2::geom_point(color = "#555555", alpha = 0.7, size = 2) +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "#D90429", linetype = "dashed", linewidth = 1) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(
          title = "Ajuste do Modelo Beta",
          subtitle = "Proporção Observada vs Proporção Predita",
          x = "Proporção Predita",
          y = "Proporção Observada"
        )
        
      if (exists("tema_estatR", mode = "function")) {
        p <- p + tema_estatR(estilo = 2)
      } else {
        p <- p + ggplot2::theme_minimal()
      }
      
      print(p)
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o gráfico: ", e$message)
    })
  }
  
  invisible(modelo)
}

#' @title Metricas do Modelo Beta
#' @description Calcula e exibe as principais metricas para um modelo betareg.
#' @param modelo Objeto de regressao beta (betareg)
#' @export
metricas_beta <- function(modelo) {
  if (!inherits(modelo, "betareg")) {
    stop("O modelo fornecido nao e um modelo de Regressao Beta ('betareg').")
  }
  
  resumo <- summary(modelo)
  pseudo_r2 <- resumo$pseudo.r.squared
  if (is.na(pseudo_r2)) pseudo_r2 <- 0
  
  loglik_mod <- as.numeric(logLik(modelo))
  mod_nulo <- suppressWarnings(tryCatch(
    betareg::betareg(stats::reformulate("1", as.character(formula(modelo)[[2]])), data = modelo$model, link = "logit"),
    error = function(e) NULL
  ))
  
  if (!is.null(mod_nulo)) {
    loglik_nulo <- as.numeric(logLik(mod_nulo))
    lrt_stat <- 2 * (loglik_mod - loglik_nulo)
    lrt_df <- length(coef(modelo)) - length(coef(mod_nulo))
    p_lrt <- pchisq(lrt_stat, lrt_df, lower.tail = FALSE)
  } else {
    lrt_stat <- NA; p_lrt <- NA
  }
  
  phi_est <- resumo$coefficients$precision[1, 1]
  
  aic_val <- AIC(modelo)
  bic_val <- BIC(modelo)
  k <- length(coef(modelo))
  n <- nobs(modelo)
  aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
  
  .fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("MÉTRICAS DO MODELO BETA")
  } else {
    cat("\n── MÉTRICAS DO MODELO BETA ──\n")
  }
  
  w_m <- c(met=38, val=15)
  linha <- paste0("  ", strrep("─", sum(w_m) + 1))
  
  cat(linha, "\n")
  cat("  ", .pad("Métrica", w_m["met"], "left"), .pad("Valor", w_m["val"], "center"), "\n", sep = "")
  cat(linha, "\n")
  cat("  ", .pad("Log-Verossimilhança", w_m["met"], "left"), .pad(sprintf("%.2f", loglik_mod), w_m["val"], "center"), "\n", sep = "")
  if (!is.na(p_lrt)) {
    cat("  ", .pad("Teste LRT (χ²)", w_m["met"], "left"), .pad(sprintf("%.2f (p %s)", lrt_stat, ifelse(p_lrt < 0.001, "<0,001", paste0("=", .fmt_p(p_lrt)))), w_m["val"], "center"), "\n", sep = "")
  }
  cat("  ", .pad("Parâmetro Precisão (Phi)", w_m["met"], "left"), .pad(sprintf("%.2f", phi_est), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Pseudo-R²", w_m["met"], "left"), .pad(sprintf("%.1f%%", pseudo_r2 * 100), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("AIC", w_m["met"], "left"), .pad(sprintf("%.2f", aic_val), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("AICc", w_m["met"], "left"), .pad(sprintf("%.2f", aicc_val), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("BIC", w_m["met"], "left"), .pad(sprintf("%.2f", bic_val), w_m["val"], "center"), "\n", sep = "")
  if (exists(".print_rodape", mode = "function")) .print_rodape() else cat(linha, "\n")
  
  invisible(list(loglik = loglik_mod, pseudo_r2 = pseudo_r2, phi = phi_est, aic = aic_val))
}

#' @title Analise Residual para Modelo Beta
#' @description Realiza o diagnostico de residuos para regressao Beta (residuos padronizados ponderados).
#' @param modelo Objeto betareg
#' @param grafico Logico. Plota o painel se TRUE.
#' @export
analise_residual_beta <- function(modelo, grafico = TRUE) {
  if (!inherits(modelo, "betareg")) stop("O modelo fornecido nao e 'betareg'.")
  
  if (exists(".print_titulo", mode = "function")) .print_titulo("ANÁLISE DE RESÍDUOS — REGRESSÃO BETA")
  
  # Usar residuos "sweighted2" (padronizados ponderados) recomendados para betareg
  res_pad <- residuals(modelo, type = "sweighted2")
  cooks_d <- cooks.distance(modelo)
  n <- nobs(modelo)
  corte_cook <- 4 / n
  infl_pts <- which(cooks_d > corte_cook)
  
  if (exists(".print_topico", mode = "function")) .print_topico("PONTOS INFLUENTES (Distância de Cook)")
  if (length(infl_pts) > 0) {
    cat(sprintf("  Observações com distância > %.4f (4/n):\n", corte_cook))
    for (i in infl_pts) {
      cat(sprintf("  Obs. %s: Cook = %.4f | Res. Padronizado = %.3f\n", names(cooks_d)[i], cooks_d[i], res_pad[i]))
    }
    cat("  [!] Avalie remover essas observações ou investigar outliers.\n\n")
  } else {
    cat("  [✓] Nenhuma observação excessivamente influente detectada.\n\n")
  }
  
  if (exists(".print_rodape", mode = "function")) .print_rodape()
  
  if (grafico) {
    tryCatch({
      df_plot <- data.frame(
        Ajustados = predict(modelo, type = "response"),
        Residuos  = res_pad,
        CooksD    = cooks_d,
        PredLin   = predict(modelo, type = "link"),
        Index     = seq_len(n)
      )
      
      p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Ajustados, y = Residuos)) +
        ggplot2::geom_point(color = "#555555", alpha = 0.6) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#D90429") +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429", linewidth = 1) +
        ggplot2::labs(title = "Resíduos vs Ajustados", x = "Proporção Ajustada", y = "Resíduo Padronizado")
        
      p2 <- ggplot2::ggplot(df_plot, ggplot2::aes(sample = Residuos)) +
        ggplot2::stat_qq(color = "#555555", alpha = 0.6) +
        ggplot2::stat_qq_line(color = "#D90429", linewidth = 1) +
        ggplot2::labs(title = "QQ-Plot (Resíduos)", x = "Quantis Teóricos", y = "Quantis Amostrais")
        
      p3 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Index, y = CooksD)) +
        ggplot2::geom_segment(ggplot2::aes(xend = Index, yend = 0), color = "#555555") +
        ggplot2::geom_point(color = "#555555") +
        ggplot2::geom_hline(yintercept = corte_cook, linetype = "dashed", color = "#D90429") +
        ggplot2::labs(title = "Distância de Cook", x = "Índice da Observação", y = "Distância")
        
      p4 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = PredLin, y = abs(Residuos))) +
        ggplot2::geom_point(color = "#555555", alpha = 0.6) +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429", linewidth = 1) +
        ggplot2::labs(title = "Escala-Localização", x = "Preditor Linear (Logit)", y = "|Resíduo Padronizado|")
        
      if (exists("tema_estatR", mode = "function")) {
        p1 <- p1 + tema_estatR(estilo = 2); p2 <- p2 + tema_estatR(estilo = 2)
        p3 <- p3 + tema_estatR(estilo = 2); p4 <- p4 + tema_estatR(estilo = 2)
      } else {
        p1 <- p1 + ggplot2::theme_minimal(); p2 <- p2 + ggplot2::theme_minimal()
        p3 <- p3 + ggplot2::theme_minimal(); p4 <- p4 + ggplot2::theme_minimal()
      }
      
      print((p1 | p2) / (p3 | p4))
    }, error = function(e) {
      message("[Aviso] Erro ao gerar painel de resíduos: ", e$message)
    })
  }
  invisible(list(cook = cooks_d))
}

#' @title Selecao de Modelos - Beta
#' @description Realiza Best Subsets para regressao Beta.
#' @param formula Formula do modelo
#' @param dados Data frame
#' @param top Numero de top modelos a retornar
#' @export
selecao_modelos_beta <- function(formula, dados, top = 5) {
  if (!requireNamespace("betareg", quietly = TRUE)) stop("Pacote 'betareg' e necessario.")
  
  if (exists(".print_titulo", mode = "function")) .print_titulo("SELEÇÃO DE MODELOS — REGRESSÃO BETA")
  
  y_name <- as.character(formula[[2]])
  preditores <- attr(terms(formula, data = dados), "term.labels")
  
  if (length(preditores) > 10) {
    cat(sprintf("  Muitos preditores (%d). Usando seleção heurística...\n\n", length(preditores)))
    p_vals <- sapply(preditores, function(p) {
       f <- as.formula(paste(y_name, "~", p))
       mod <- suppressWarnings(tryCatch(betareg::betareg(f, data = dados), error=function(e) NULL))
       if (is.null(mod)) return(1.0)
       coef(summary(mod))$mean[2, 4]
    })
    preditores <- names(sort(p_vals)[1:10])
  }
  
  cat(sprintf("  Variável Resposta: %s\n", y_name))
  cat(sprintf("  Preditores base:   %s\n", paste(preditores, collapse = ", ")))
  
  combinacoes <- list()
  for (i in 1:length(preditores)) {
    combinacoes <- c(combinacoes, combn(preditores, i, simplify = FALSE))
  }
  
  cat(sprintf("  Método utilizado:  Exaustivo (%d modelos avaliados)\n\n", length(combinacoes)))
  
  resultados <- data.frame(Modelo = character(), k = integer(), AIC = numeric(), AICc = numeric(), Pseudo_R2 = numeric(), stringsAsFactors = FALSE)
  
  n_obs <- nrow(dados)
  
  for (vars in combinacoes) {
    f_str <- paste(y_name, "~", paste(vars, collapse = " + "))
    mod <- suppressWarnings(tryCatch(betareg::betareg(as.formula(f_str), data = dados), error=function(e) NULL))
    
    if (!is.null(mod)) {
      k <- length(vars) + 1 # +1 pro parametro de precisao phi
      aic_val <- AIC(mod)
      aicc_val <- aic_val + (2 * (k + 1) * (k + 2)) / (n_obs - (k + 1) - 1)
      pseudo <- summary(mod)$pseudo.r.squared
      if (is.na(pseudo)) pseudo <- 0
      
      resultados <- rbind(resultados, data.frame(
        Modelo = paste(vars, collapse = " + "),
        k = k,
        AIC = aic_val,
        AICc = aicc_val,
        Pseudo_R2 = pseudo,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  resultados <- resultados[order(resultados$AIC), ]
  df_top <- head(resultados, top)
  df_top$Ranking <- paste0(1:nrow(df_top), "º")
  df_top <- df_top[, c("Ranking", "Modelo", "k", "AIC", "AICc", "Pseudo_R2")]
  
  df_fmt <- df_top
  names(df_fmt)[names(df_fmt) == "Pseudo_R2"] <- "Pseudo-R²"
  df_fmt$AIC <- sprintf("%.1f", df_top$AIC)
  df_fmt$AICc <- sprintf("%.1f", df_top$AICc)
  df_fmt[["Pseudo-R²"]] <- sprintf("%.1f%%", df_top$Pseudo_R2 * 100)
  
  if (exists(".print_topico", mode = "function")) .print_topico(sprintf("TOP %d MODELOS (Ordenados por AIC)", nrow(df_fmt)))
  
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_fmt, align = c("left", "left", "center", "center", "center", "center"))
  } else {
    print(df_fmt, row.names = FALSE)
  }
  
  cat("  * k = Número de parâmetros preditores (+1 do Phi).\n")
  cat("  * AIC menor = modelo mais parcimonioso.\n\n")
  if (exists(".print_rodape", mode = "function")) .print_rodape()
  
  invisible(df_top)
}
