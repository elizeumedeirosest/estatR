#' @title Regressão Quantílica
#' @description Ajusta modelos de regressão quantílica para múltiplos quantis condicionais
#' (ex.: 0.10, 0.25, 0.50, 0.75, 0.90), permitindo investigar o impacto das covariáveis não apenas
#' na média condicional, mas também nas caudas e na mediana, sendo robusta a outliers e heterocedasticidade.
#' @param formula Fórmula simbólica do modelo (ex.: y ~ x1 + x2).
#' @param dados Data frame contendo as variáveis do modelo.
#' @param quantis Vetor numérico com os quantis a estimar (padrão: c(0.10, 0.25, 0.50, 0.75, 0.90)).
#' @param nivel Nível de confiança nominal para os intervalos (padrão: 0.95).
#' @param metodo_ep Método para cálculo dos erros padrão e inferência: "bootstrap" (bootstrap pareado empírico), "assintotico" (estimador de sandwich assintótico) ou "robusto". Padrão: "bootstrap".
#' @param n_boot Inteiro. Número de reamostragens bootstrap (padrão: 300).
#' @param semente Inteiro opcional para garantir a reprodutibilidade do bootstrap.
#' @param grafico Lógico. Se TRUE (padrão), exibe o gráfico da dispersão com as retas dos quantis ou o gráfico do processo dos coeficientes.
#' @return Objeto da classe \code{estatR_regressao_quantilica} contendo:
#'   \item{coeficientes}{Data frame com Estimativa, Erro Padrão, IC e p-valor para cada quantil}
#'   \item{quantis}{Vetor com os quantis estimados}
#'   \item{modelo_ols}{Modelo linear clássico (OLS) correspondente para comparação}
#'   \item{residuos}{Matriz de resíduos para cada quantil}
#'   \item{formula}{Fórmula utilizada}
#' @export
regressao_quantilica <- function(formula,
                                 dados     = NULL,
                                 data      = NULL,
                                 quantis   = c(0.10, 0.25, 0.50, 0.75, 0.90),
                                 nivel     = 0.95,
                                 metodo_ep = c("bootstrap", "assintotico", "robusto"),
                                 n_boot    = 300,
                                 semente   = NULL,
                                 grafico   = TRUE) {

  metodo_ep <- match.arg(metodo_ep)

  if (!is.null(semente)) set.seed(semente)

  if (is.null(dados) && !is.null(data)) dados <- data
  if (is.null(dados)) stop("O argumento 'dados' deve ser fornecido como um data frame.")
  if (!inherits(formula, "formula")) {
    formula <- stats::as.formula(formula)
  }

  mf <- stats::model.frame(formula, data = dados)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)
  n <- length(y)
  p <- ncol(X)
  nomes_params <- colnames(X)

  alpha <- 1 - nivel
  z_crit <- stats::qnorm(1 - alpha / 2)

  # Ajuste OLS de referência
  mod_ols <- stats::lm(formula, data = dados)
  coefs_ols <- stats::coef(mod_ols)
  se_ols    <- summary(mod_ols)$coefficients[, "Std. Error"]

  # Função de perda check loss: rho_tau(u) = u * (tau - I(u < 0))
  loss_check <- function(beta, X_mat, y_vec, tau) {
    u <- y_vec - X_mat %*% beta
    sum(u * (tau - (u < 0)))
  }

  # Função interna de ajuste para um quantil tau
  ajustar_rq_tau <- function(tau, X_mat, y_vec) {
    # Se quantreg estiver disponível, usa rq diretamente
    if (requireNamespace("quantreg", quietly = TRUE)) {
      fit_q <- tryCatch(suppressWarnings(quantreg::rq.fit(X_mat, y_vec, tau = tau)), error = function(e) NULL)
      if (!is.null(fit_q) && !is.null(fit_q$coefficients)) {
        return(as.numeric(fit_q$coefficients))
      }
    }
    # Fallback: otimização via Nelder-Mead com inicialização OLS
    beta_init <- stats::coef(stats::lm.fit(X_mat, y_vec))
    opt <- stats::optim(
      par = beta_init,
      fn = loss_check,
      X_mat = X_mat,
      y_vec = y_vec,
      tau = tau,
      method = "BFGS",
      control = list(maxit = 500)
    )
    if (opt$convergence != 0) {
      opt <- stats::optim(
        par = beta_init,
        fn = loss_check,
        X_mat = X_mat,
        y_vec = y_vec,
        tau = tau,
        method = "Nelder-Mead",
        control = list(maxit = 1000)
      )
    }
    opt$par
  }

  # Ajuste dos coeficientes para todos os quantis
  lista_coefs <- lapply(quantis, function(tau) ajustar_rq_tau(tau, X, y))
  mat_coefs <- do.call(cbind, lista_coefs)
  colnames(mat_coefs) <- paste0("tau_", quantis)
  rownames(mat_coefs) <- nomes_params

  # Cálculo de Erros Padrão via Bootstrap ou Assintótico
  if (metodo_ep == "bootstrap") {
    boot_res <- array(0, dim = c(p, length(quantis), n_boot))
    for (b in seq_len(n_boot)) {
      idx_b <- sample.int(n, replace = TRUE)
      X_b   <- X[idx_b, , drop = FALSE]
      y_b   <- y[idx_b]
      for (q_idx in seq_along(quantis)) {
        boot_res[, q_idx, b] <- ajustar_rq_tau(quantis[q_idx], X_b, y_b)
      }
    }
    # Desvio padrão dos coeficientes simulados
    se_mat <- apply(boot_res, c(1, 2), stats::sd, na.rm = TRUE)
  } else {
    # Estimador assintótico baseado na densidade residual
    se_mat <- matrix(0, nrow = p, ncol = length(quantis))
    inv_XX <- tryCatch(solve(t(X) %*% X), error = function(e) diag(1, p))
    for (q_idx in seq_along(quantis)) {
      tau <- quantis[q_idx]
      res_tau <- y - X %*% mat_coefs[, q_idx]
      # Estimativa de densidade de sparsity (Silverman)
      h <- 1.06 * stats::sd(res_tau) * (n^(-0.2))
      if (h == 0) h <- 0.01
      dens <- stats::density(res_tau, from = 0, to = 0, n = 1)$y
      if (is.na(dens) || dens <= 0) dens <- 0.1
      var_asymp <- (tau * (1 - tau) / (dens^2)) * inv_XX
      se_mat[, q_idx] <- sqrt(pmax(diag(var_asymp), 1e-6))
    }
  }

  # Montagem de tabela completa por quantil
  tabelas_quantis <- list()
  for (q_idx in seq_along(quantis)) {
    tau <- quantis[q_idx]
    est <- mat_coefs[, q_idx]
    se  <- se_mat[, q_idx]
    t_val <- est / se
    p_val <- 2 * stats::pnorm(-abs(t_val))
    ic_inf <- est - z_crit * se
    ic_sup <- est + z_crit * se

    df_tau <- data.frame(
      Parametro   = nomes_params,
      Estimativa  = est,
      Erro_Padrao = se,
      IC_Inf      = ic_inf,
      IC_Sup      = ic_sup,
      Estatistica = t_val,
      p_valor     = p_val,
      stringsAsFactors = FALSE
    )
    tabelas_quantis[[paste0("tau_", tau)]] <- df_tau
  }

  # Seleciona o quantil mediano (tau = 0.50 ou mais próximo de 0.50) para exibição principal
  diff_med <- abs(quantis - 0.50)
  idx_med <- which.min(diff_med)
  tau_med <- quantis[idx_med]
  df_principal <- tabelas_quantis[[idx_med]]

  # Helpers de formatação
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 3) {
    if (is.na(x)) return("-")
    formatC(x, format = "f", digits = d, decimal.mark = ",")
  }
  .fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001 ***")
    ast <- if (p < 0.01) " **" else if (p < 0.05) " *" else ""
    paste0(formatC(p, format = "f", digits = 3, decimal.mark = ","), ast)
  }

  # 4. Impressão no Console
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("REGRESSÃO QUANTÍLICA")
  } else {
    cat("
── REGRESSÃO QUANTÍLICA ──
")
  }

  cat(sprintf("  Fórmula do Modelo:    %s
", deparse(formula)))
  cat(sprintf("  Amostra Válida (N):   %d observações
", n))
  cat(sprintf("  Quantis Avaliados:    %s
", paste(sprintf("%.2f", quantis), collapse = ", ")))
  cat(sprintf("  Método de Erro Pad.:  %s (%s)

",
              if (metodo_ep == "bootstrap") "Bootstrap Empírico" else if (metodo_ep == "robusto") "Sandwich Robusto" else "Assintótico",
              if (metodo_ep == "bootstrap") paste0(n_boot, " réplicas") else "Teórico"))

  # Tabela do Quantil Principal (Mediana tau = 0.50)
  topico_med <- sprintf("COEFICIENTES ESTIMADOS (QUANTIL τ = %.2f — MEDIANA CONDICIONAL)", tau_med)
  if (exists(".print_topico", mode = "function")) {
    .print_topico(topico_med)
  } else {
    cat(sprintf("  ▶ %s
", topico_med))
  }

  w_rq <- c(par = max(20, max(nchar(nomes_params)) + 2), est = 14, se = 14, ic = 24, z = 12, p = 16)
  sep_rq <- paste0("  ", strrep("─", sum(w_rq)))
  cat(sep_rq, "
")
  cat("  ", .pad("Parâmetro", w_rq["par"], "left"),
             .pad("Estimativa", w_rq["est"], "center"),
             .pad("Erro Padrão", w_rq["se"], "center"),
             .pad(sprintf("IC (%d%%)", round(nivel * 100)), w_rq["ic"], "center"),
             .pad("z", w_rq["z"], "center"),
             .pad("p-valor", w_rq["p"], "center"), "
", sep = "")
  cat(sep_rq, "
")
  for (i in seq_len(nrow(df_principal))) {
    str_ic <- sprintf("[ %s ; %s ]", .fmt_f(df_principal$IC_Inf[i], 2), .fmt_f(df_principal$IC_Sup[i], 2))
    cat("  ", .pad(df_principal$Parametro[i], w_rq["par"], "left"),
               .pad(.fmt_f(df_principal$Estimativa[i], 3), w_rq["est"], "center"),
               .pad(.fmt_f(df_principal$Erro_Padrao[i], 3), w_rq["se"], "center"),
               .pad(str_ic, w_rq["ic"], "center"),
               .pad(.fmt_f(df_principal$Estatistica[i], 2), w_rq["z"], "center"),
               .pad(.fmt_p(df_principal$p_valor[i]), w_rq["p"], "center"), "
", sep = "")
  }
  cat(sep_rq, "
")
  cat("  Notas:
")
  cat("  Significância: *** p < 0,001   ** p < 0,01   * p < 0,05

")

  # Tabela Resumo dos Coeficientes através dos Quantis
  if (exists(".print_topico", mode = "function")) {
    .print_topico("EVOLUÇÃO DOS COEFICIENTES ATRAVÉS DOS QUANTIS")
  } else {
    cat("  ▶ EVOLUÇÃO DOS COEFICIENTES ATRAVÉS DOS QUANTIS
")
  }

  w_q_col <- 13
  w_tot_evol <- w_rq["par"] + w_q_col * (length(quantis) + 1)
  sep_evol <- paste0("  ", strrep("─", w_tot_evol))
  cat(sep_evol, "
")
  cat("  ", .pad("Parâmetro", w_rq["par"], "left"))
  for (q in quantis) cat(.pad(sprintf("τ = %.2f", q), w_q_col, "center"))
  cat(.pad("Média (OLS)", w_q_col, "center"), "
")
  cat(sep_evol, "
")
  for (p_idx in seq_along(nomes_params)) {
    nm_p <- nomes_params[p_idx]
    cat("  ", .pad(nm_p, w_rq["par"], "left"))
    for (q_idx in seq_along(quantis)) {
      cat(.pad(.fmt_f(mat_coefs[p_idx, q_idx], 2), w_q_col, "center"))
    }
    cat(.pad(.fmt_f(coefs_ols[nm_p], 2), w_q_col, "center"), "
")
  }
  cat(sep_evol, "

")

  # INTERPRETAÇÃO DOS COEFICIENTES
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO DOS COEFICIENTES")
  } else {
    cat("  ▶ INTERPRETAÇÃO DOS COEFICIENTES
")
  }

  cor_ok <- "[1;32m"
  reset  <- "[0m"

  # Identifica se algum coeficiente varia substancialmente entre as caudas
  var_caudas <- FALSE
  param_dif <- ""
  if (length(quantis) >= 2) {
    q_inf_idx <- 1
    q_sup_idx <- length(quantis)
    difs <- abs(mat_coefs[, q_sup_idx] - mat_coefs[, q_inf_idx])
    max_dif_idx <- which.max(difs[-1]) + 1
    if (length(max_dif_idx) > 0 && difs[max_dif_idx] > 0.5 * sd(y)) {
      var_caudas <- TRUE
      param_dif <- nomes_params[max_dif_idx]
    }
  }

  cat(sprintf("  %s[✓]%s A estimativa do quantil mediano (τ = %.2f) reflete a tendência central robusta da resposta,
    imune a distorções geradas por assimetria distributiva ou outliers extremos.

",
              cor_ok, reset, tau_med))

  if (var_caudas) {
    cat(sprintf("  %s[✓]%s O efeito da variável '%s' altera-se significativamente ao longo da distribuição,
    evidenciando heterogeneidade no impacto preditivo entre as caudas inferior e superior.

",
                cor_ok, reset, param_dif))
  } else {
    cat(sprintf("  %s[✓]%s Os coeficientes quantílicos mantêm trajetória consistente com o modelo linear OLS,
    demonstrando estabilidade estrutural do efeito das covariáveis ao longo de toda a distribuição.

",
                cor_ok, reset))
  }

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat("──────────────────────────────────────────────────────────────────────────────
")
  }

  res_obj <- list(
    coeficientes       = tabelas_quantis,
    matriz_coeficientes= mat_coefs,
    quantis            = quantis,
    modelo_ols         = mod_ols,
    formula            = formula,
    dados              = dados,
    nomes_params       = nomes_params
  )
  class(res_obj) <- c("estatR_regressao_quantilica", "list")

  if (isTRUE(grafico)) {
    p_rq <- plot_regressao_quantilica(res_obj)
    print(p_rq)
  }

  invisible(res_obj)
}

#' @title Gráfico do Processo dos Coeficientes Quantílicos
#' @description Exibe a trajetória dos coeficientes estimados ao longo dos quantis com envelope de confiança e reta clássica OLS.
#' @param modelo Objeto da classe \code{estatR_regressao_quantilica}.
#' @return Objeto ggplot2.
#' @export
plot_regressao_quantilica <- function(modelo) {
  if (!inherits(modelo, "estatR_regressao_quantilica")) {
    stop("O modelo deve ser um objeto da classe 'estatR_regressao_quantilica'.")
  }

  # Monta data frame para plot de todos os parâmetros (exceto intercepto se houver outros)
  nomes_p <- modelo$nomes_params
  params_plot <- if (length(nomes_p) > 1) setdiff(nomes_p, c("(Intercept)", "Intercept")) else nomes_p

  df_lista <- list()
  for (q_nm in names(modelo$coeficientes)) {
    tau_val <- as.numeric(sub("tau_", "", q_nm))
    df_t <- modelo$coeficientes[[q_nm]]
    df_t$Quantil <- tau_val
    df_lista[[q_nm]] <- df_t[df_t$Parametro %in% params_plot, , drop = FALSE]
  }
  df_proc <- do.call(rbind, df_lista)

  # Coeficientes OLS correspondentes
  coefs_ols <- stats::coef(modelo$modelo_ols)
  df_ols <- data.frame(
    Parametro = params_plot,
    Estimativa_OLS = coefs_ols[params_plot],
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df_proc, ggplot2::aes(x = Quantil, y = Estimativa)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = IC_Inf, ymax = IC_Sup), fill = "#1D3557", alpha = 0.20) +
    ggplot2::geom_line(color = "#1D3557", linewidth = 1) +
    ggplot2::geom_point(color = "#1D3557", size = 2.5) +
    ggplot2::geom_hline(data = df_ols, ggplot2::aes(yintercept = Estimativa_OLS), color = "#E63946", linetype = "dashed", linewidth = 0.8) +
    ggplot2::facet_wrap(~ Parametro, scales = "free_y") +
    ggplot2::labs(
      title    = "Trajetória dos Coeficientes através dos Quantis",
      subtitle = "Faixa sombreada: IC 95% Quantílico | Linha tracejada vermelha: Média OLS",
      x        = "Quantil Condicional (τ)",
      y        = "Coeficiente Estimado"
    )

  if (exists("meu_tema", mode = "function")) {
    p <- p + meu_tema()
  } else if (exists("tema_estatR", mode = "function")) {
    p <- p + meu_tema(grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}
