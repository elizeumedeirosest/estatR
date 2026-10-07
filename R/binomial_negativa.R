#' @title Regressao Binomial Negativa
#' @description Ajusta um modelo de regressao Binomial Negativa para dados de contagem
#' com sobredispersao (variancia > media). Apresenta estatisticas de ajuste,
#' parametro de dispersao theta, coeficientes, Razao de Taxa de Incidencia (IRR)
#' e interpretacao automatizada.
#' @param formula Formula do modelo (ex: Y ~ X1 + X2).
#' @param dados Data frame contendo as variaveis.
#' @param grafico Logico. Se TRUE, plota o grafico de Observado vs Predito.
#' @export
regressao_binomial_negativa <- function(formula, dados, grafico = TRUE) {
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("O pacote 'MASS' e necessario para ajustar modelos binomiais negativos. Instale com install.packages('MASS').")
  }
  
  # 1. Validacao e Ajuste
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  
  # Capturando possiveis erros/warnings no ajuste (ex: nao convergencia)
  modelo <- tryCatch(
    MASS::glm.nb(formula, data = dados, link = log),
    error = function(e) stop("Erro ao ajustar o modelo Binomial Negativo: ", e$message)
  )
  
  y_name <- as.character(formula[[2]])
  x_names <- attr(terms(formula), "term.labels")
  if (length(x_names) == 0) stop("O modelo deve ter pelo menos um preditor.")
  
  # 2. Extrair Resultados
  resumo <- summary(modelo)
  coefs <- resumo$coefficients
  theta_val <- modelo$theta
  theta_se  <- modelo$SE.theta
  
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
  
  # Para o modelo nulo Binomial Negativo, usamos o mesmo glm.nb
  mod_nulo <- suppressWarnings(MASS::glm.nb(stats::reformulate("1", y_name), data = dados, link = log))
  loglik_nulo <- as.numeric(logLik(mod_nulo))
  pseudo_r2   <- 1 - (loglik_mod / loglik_nulo)
  
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
    .print_titulo("REGRESSÃO BINOMIAL NEGATIVA")
  } else {
    cat("\n── REGRESSÃO BINOMIAL NEGATIVA ──\n")
  }
  
  # DIAGNOSTICO DO AJUSTE
  if (exists(".print_topico", mode = "function")) .print_topico("DIAGNÓSTICO DO AJUSTE")
  
  w_diag <- c(met = 38, val = 15)
  sep_diag <- paste0("  ", strrep("─", sum(w_diag) + 3))
  cat(sep_diag, "\n")
  cat("  ", .pad("Métrica", w_diag["met"], "left"), .pad("Valor", w_diag["val"], "center"), "\n", sep = "")
  cat(sep_diag, "\n")
  cat("  ", .pad("Deviance Nula", w_diag["met"], "left"), .pad(sprintf("%.2f (df=%d)", dev_nula, df_nula), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Deviance Residual", w_diag["met"], "left"), .pad(sprintf("%.2f (df=%d)", dev_res, df_res), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Parâmetro Theta (θ)", w_diag["met"], "left"), .pad(sprintf("%.3f", theta_val), w_diag["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Teste LRT (χ²)", w_diag["met"], "left"), .pad(sprintf("%.2f (%s)", lrt_stat, .fmt_p(p_lrt)), w_diag["val"], "center"), .ast(p_lrt), "\n", sep = "")
  cat("  ", .pad("Pseudo-R² (McFadden)", w_diag["met"], "left"), .pad(sprintf("%.1f%%", pseudo_r2 * 100), w_diag["val"], "center"), "\n", sep = "")
  cat(sep_diag, "\n\n")
  
  # COEFICIENTES E IRR
  if (exists(".print_topico", mode = "function")) .print_topico("COEFICIENTES E RAZÃO DE TAXA DE INCIDÊNCIA (IRR)")
  
  w_c <- c(var=14, est=12, err=11, z=9, pval=14, irr=10, ic_inf=12, ic_sup=12)
  sep_coef <- paste0("  ", strrep("─", sum(w_c) + length(w_c) - 1))
  
  cat(sep_coef, "\n")
  cat("  ",
      .pad("Variável", w_c["var"], "left"),
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
  cat("  Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n")
  cat("  IRR > 1: aumenta a contagem esperada. IRR < 1: reduz a contagem esperada.\n\n")
  
  # INTERPRETACAO
  if (exists(".print_topico", mode = "function")) .print_topico("INTERPRETAÇÃO")
  
  cat(sprintf("  O modelo Binomial Negativo foi ajustado para estimar a\n  contagem de '%s', acomodando possível sobredispersão.\n\n", y_name))
  
  if (p_lrt < 0.05) {
    cat(sprintf("  O modelo é globalmente válido (LRT %s) e apresenta um\n  Pseudo-R² de McFadden de %.1f%%, indicando a proporção da variabilidade\n  explicada pelos preditores incluídos.\n\n", 
                ifelse(p_lrt < 0.001, "p < 0,001", sprintf("p = %s", .fmt_p(p_lrt))), pseudo_r2 * 100))
  } else {
    cat(sprintf("  O modelo não se mostrou globalmente significativo (LRT p = %s).\n  Os preditores não explicam a variabilidade da contagem melhor que um modelo nulo.\n\n", .fmt_p(p_lrt)))
  }
  
  cat(sprintf("  O parâmetro de dispersão (Theta = %.3f, EP = %.3f) reflete a \n  variabilidade adicional nos dados além da média (variância = μ + μ²/θ).\n\n", theta_val, theta_se))
  
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
      
      cat(sprintf("  O preditor '%s' apresentou efeito estatisticamente\n  significativo (p %s), mantendo os demais constantes.\n  Um aumento de uma unidade em '%s' %s a contagem esperada de '%s'\n  em %.1f%% (IRR = %s; IC 95%%: [%s; %s]).\n\n",
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
    cat(strrep("─", sum(w_c) + length(w_c) - 1), "\n\n")
  }
  
  # 5. Grafico: Observado vs Predito
  if (grafico) {
    tryCatch({
      df_plot <- data.frame(
        Observado = dados[[y_name]],
        Predito = modelo$fitted.values
      )
      
      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Predito, y = Observado)) +
        ggplot2::geom_point(color = "#555555", alpha = 0.7, size = 2) +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "#D90429", linetype = "dashed", linewidth = 1) +
        ggplot2::labs(
          title = "Ajuste do Modelo Binomial Negativo",
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
      message("[Aviso] Não foi possível gerar o gráfico: ", e$message)
    })
  }
  
  invisible(modelo)
}

#' @title Metricas do Modelo Binomial Negativo
#' @description Calcula e exibe as principais metricas para um modelo glm.nb.
#' @param modelo Objeto glm.nb (negbin)
#' @export
metricas_binomial_negativa <- function(modelo) {
  if (!inherits(modelo, "negbin")) {
    stop("O modelo fornecido nao e um modelo Binomial Negativo (glm.nb).")
  }
  
  # Deviance e df
  dev_nula <- modelo$null.deviance
  df_nula  <- modelo$df.null
  dev_res  <- modelo$deviance
  df_res   <- modelo$df.residual
  theta_val <- modelo$theta
  
  # LRT
  lrt_stat <- dev_nula - dev_res
  lrt_df   <- df_nula - df_res
  p_lrt    <- pchisq(lrt_stat, lrt_df, lower.tail = FALSE)
  
  # Pseudo-R2
  loglik_mod  <- as.numeric(logLik(modelo))
  mod_nulo    <- suppressWarnings(MASS::glm.nb(stats::reformulate("1", as.character(formula(modelo)[[2]])), 
                                               data = modelo$model, link = log))
  loglik_nulo <- as.numeric(logLik(mod_nulo))
  pseudo_r2   <- 1 - (loglik_mod / loglik_nulo)
  
  # Info Criterions
  aic_val <- AIC(modelo)
  bic_val <- BIC(modelo)
  k <- length(coef(modelo)) + 1 # +1 para o theta
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
    .print_titulo("MÉTRICAS DO MODELO BINOMIAL NEGATIVO")
  } else {
    cat("\n── MÉTRICAS DO MODELO BINOMIAL NEGATIVO ──\n")
  }
  
  w_m <- c(met=38, val=15)
  linha <- paste0("  ", strrep("─", sum(w_m) + 1))
  
  cat(linha, "\n")
  cat("  ", .pad("Métrica", w_m["met"], "left"), .pad("Valor", w_m["val"], "center"), "\n", sep = "")
  cat(linha, "\n")
  cat("  ", .pad("Deviance Nula", w_m["met"], "left"), .pad(sprintf("%.2f (df=%d)", dev_nula, df_nula), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Deviance Residual", w_m["met"], "left"), .pad(sprintf("%.2f (df=%d)", dev_res, df_res), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Teste LRT (χ²)", w_m["met"], "left"), .pad(sprintf("%.2f (p %s)", lrt_stat, ifelse(p_lrt < 0.001, "<0,001", paste0("=", .fmt_p(p_lrt)))), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Parâmetro Theta (θ)", w_m["met"], "left"), .pad(sprintf("%.3f", theta_val), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Pseudo-R² (McFadden)", w_m["met"], "left"), .pad(sprintf("%.1f%%", pseudo_r2 * 100), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("AIC", w_m["met"], "left"), .pad(sprintf("%.2f", aic_val), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("AICc", w_m["met"], "left"), .pad(sprintf("%.2f", aicc_val), w_m["val"], "center"), "\n", sep = "")
  cat("  ", .pad("BIC", w_m["met"], "left"), .pad(sprintf("%.2f", bic_val), w_m["val"], "center"), "\n", sep = "")
  if (exists(".print_rodape", mode = "function")) .print_rodape() else cat(linha, "\n")
  
  invisible(list(deviance = dev_res, pseudo_r2 = pseudo_r2, theta = theta_val, aic = aic_val))
}

#' @title Analise Residual para Modelo Binomial Negativo
#' @description Realiza o diagnostico de residuos para modelos de regressao Binomial Negativa.
#' @param modelo Objeto glm.nb (negbin)
#' @param grafico Logico. Plota o painel se TRUE.
#' @export
analise_residual_binomial_negativa <- function(modelo, grafico = TRUE) {
  if (exists(".print_titulo", mode = "function")) .print_titulo("ANÁLISE DE RESÍDUOS — BINOMIAL NEGATIVA")
  
  df_res <- modelo$df.residual
  pearson_res <- residuals(modelo, type = "pearson")
  
  cooks_d <- cooks.distance(modelo)
  n <- nobs(modelo)
  corte_cook <- 4 / n
  infl_pts <- which(cooks_d > corte_cook)
  
  if (exists(".print_topico", mode = "function")) .print_topico("PONTOS INFLUENTES (Distância de Cook)")
  if (length(infl_pts) > 0) {
    cat(sprintf("  Observações com distância > %.4f (4/n):\n", corte_cook))
    for (i in infl_pts) {
      cat(sprintf("  Obs. %s: Cook = %.4f | Res. Pearson = %.3f\n", names(cooks_d)[i], cooks_d[i], pearson_res[i]))
    }
    cat("  [!] Avalie remover essas observações ou investigar outliers.\n\n")
  } else {
    cat("  [✓] Nenhuma observação excessivamente influente detectada.\n\n")
  }
  
  if (exists(".print_rodape", mode = "function")) .print_rodape()
  
  if (grafico) {
    tryCatch({
      df_plot <- data.frame(
        Ajustados = modelo$fitted.values,
        Residuos  = pearson_res,
        CooksD    = cooks_d,
        PredLin   = modelo$linear.predictors,
        Index     = seq_len(n)
      )
      
      p1 <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Ajustados, y = Residuos)) +
        ggplot2::geom_point(color = "#555555", alpha = 0.6) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#D90429") +
        ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#D90429", linewidth = 1) +
        ggplot2::labs(title = "Resíduos vs Ajustados", x = "Valores Ajustados", y = "Resíduos Pearson")
        
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
        ggplot2::labs(title = "Escala-Localização", x = "Preditor Linear", y = "|Resíduo Pearson|")
        
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

#' @title Selecao de Modelos - Binomial Negativa
#' @description Realiza Best Subsets para regressao Binomial Negativa.
#' @param formula Formula do modelo
#' @param dados Data frame
#' @param top Numero de top modelos a retornar
#' @export
selecao_modelos_binomial_negativa <- function(formula, dados, top = 5) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("O pacote 'MASS' e necessario.")
  }
  
  if (exists(".print_titulo", mode = "function")) .print_titulo("SELEÇÃO DE MODELOS — BINOMIAL NEGATIVA")
  
  y_name <- as.character(formula[[2]])
  preditores <- attr(terms(formula, data = dados), "term.labels")
  
  if (length(preditores) > 10) {
    cat(sprintf("  Muitos preditores (%d). Usando seleção heurística...\n\n", length(preditores)))
    p_vals <- sapply(preditores, function(p) {
       f <- as.formula(paste(y_name, "~", p))
       mod <- suppressWarnings(tryCatch(MASS::glm.nb(f, data = dados, link = log), error=function(e) NULL))
       if (is.null(mod)) return(1.0)
       coef(summary(mod))[2, 4]
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
  
  mod_nulo <- suppressWarnings(tryCatch(MASS::glm.nb(stats::reformulate("1", y_name), data = dados, link = log), error=function(e) NULL))
  if (is.null(mod_nulo)) stop("Erro ao ajustar o modelo nulo Binomial Negativo.")
  loglik_nulo <- as.numeric(logLik(mod_nulo))
  n_obs <- nrow(dados)
  
  for (vars in combinacoes) {
    f_str <- paste(y_name, "~", paste(vars, collapse = " + "))
    mod <- suppressWarnings(tryCatch(MASS::glm.nb(as.formula(f_str), data = dados, link = log), error=function(e) NULL))
    
    if (!is.null(mod)) {
      k <- length(vars) + 1 # +1 do theta
      aic_val <- AIC(mod)
      aicc_val <- aic_val + (2 * (k + 1) * (k + 2)) / (n_obs - (k + 1) - 1)
      pseudo <- 1 - (as.numeric(logLik(mod)) / loglik_nulo)
      
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
  
  cat("  * k = Número de parâmetros preditores (+1 do Theta).\n")
  cat("  * AIC menor = modelo mais parcimonioso.\n\n")
  if (exists(".print_rodape", mode = "function")) .print_rodape()
  
  invisible(df_top)
}
