# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: CONTROLE ESTATÍSTICO DA QUALIDADE (CEQ / SIX SIGMA)
# ─────────────────────────────────────────────────────────────────────────────

# Tabela de Constantes de Shewhart para n = 2:25
.shewhart_consts <- list(
  A2 = c(1.880, 1.023, 0.729, 0.577, 0.483, 0.419, 0.373, 0.337, 0.308,
         0.285, 0.266, 0.249, 0.235, 0.223, 0.212, 0.203, 0.194, 0.187,
         0.180, 0.173, 0.167, 0.162, 0.157, 0.153),
  A3 = c(2.659, 1.954, 1.628, 1.427, 1.287, 1.182, 1.099, 1.032, 0.975,
         0.927, 0.886, 0.850, 0.817, 0.789, 0.763, 0.739, 0.718, 0.698,
         0.680, 0.663, 0.647, 0.633, 0.619, 0.606),
  D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223, 0.256, 0.283, 0.307,
         0.328, 0.347, 0.363, 0.378, 0.391, 0.403, 0.415, 0.425, 0.434, 0.443, 0.451, 0.459),
  D4 = c(3.267, 2.574, 2.282, 2.114, 2.004, 1.924, 1.864, 1.816, 1.777,
         1.744, 1.717, 1.693, 1.672, 1.653, 1.637, 1.622, 1.608, 1.597,
         1.585, 1.575, 1.566, 1.557, 1.548, 1.541),
  B3 = c(0, 0, 0, 0, 0.030, 0.118, 0.185, 0.239, 0.284, 0.321, 0.354, 0.382,
         0.406, 0.428, 0.448, 0.466, 0.482, 0.497, 0.510, 0.523, 0.534, 0.545, 0.555, 0.565),
  B4 = c(3.267, 2.568, 2.266, 2.089, 1.970, 1.882, 1.815, 1.761, 1.716,
         1.679, 1.646, 1.618, 1.594, 1.572, 1.552, 1.534, 1.518, 1.503,
         1.490, 1.477, 1.466, 1.455, 1.445, 1.435),
  d2 = c(1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.970, 3.078,
         3.173, 3.258, 3.336, 3.407, 3.472, 3.532, 3.588, 3.640, 3.689,
         3.735, 3.778, 3.819, 3.858, 3.895, 3.931),
  c4 = c(0.7979, 0.8862, 0.9213, 0.9400, 0.9515, 0.9594, 0.9650, 0.9693,
         0.9727, 0.9754, 0.9776, 0.9794, 0.9810, 0.9823, 0.9835, 0.9845,
         0.9854, 0.9862, 0.9869, 0.9876, 0.9882, 0.9887, 0.9892, 0.9896)
)

.get_const <- function(nome, n) {
  if (n < 2) return(1)
  idx <- min(max(n - 1, 1), 24)
  .shewhart_consts[[nome]][idx]
}


#' @title Gráficos de Controle Estatístico de Processo (Cartas de Shewhart)
#' @description Constrói cartas de controle de Shewhart para monitoramento de estabilidade
#'   e controle de processos contínuos e por atributos, com detecção formal de causas
#'   especiais (regras de Western Electric / Nelson) e visualização gráfica com grade dupla.
#'
#' @param dados Vetor numérico ou matriz/data.frame com as medições de processo.
#' @param tipo Tipo de carta de controle:
#'   \itemize{
#'     \item \code{"imr"}: Carta de Valores Individuais e Amplitude Móvel (I-MR).
#'     \item \code{"xbar_r"}: Carta de Média e Amplitude (\eqn{\bar{X}}-R).
#'     \item \code{"xbar_s"}: Carta de Média e Desvio-Padrão (\eqn{\bar{X}}-S).
#'     \item \code{"p"}: Proporção de não-conformes (tamanho amostral variável ou constante).
#'     \item \code{"np"}: Número de unidades defeituosas (tamanho amostral constante).
#'     \item \code{"c"}: Contagem total de defeitos por amostra (tamanho constante).
#'     \item \code{"u"}: Taxa de defeitos por unidade de inspeção (tamanho variável).
#'   }
#' @param subgrupo Vetor indicando a identificação dos subgrupos racionais (opcional).
#' @param tamanho_amostra Tamanho dos subgrupos \eqn{n} para cartas por atributos (\code{p}, \code{np}, \code{u}).
#' @param regras Lógico. Se \code{TRUE} (padrão), testa as regras de causas especiais e destaca pontos fora de controle.
#' @param grafico Lógico. Se \code{TRUE} (padrão), renderiza os gráficos de controle no padrão do estatR.
#'
#' @return Retorna um objeto da classe \code{estatR_ceq} contendo:
#' \itemize{
#'   \item \code{limites}: Tabela de Limites Estatísticos (LSC, Linha Média, LIC).
#'   \item \code{pontos}: Tabela com os valores observados de cada subgrupo e status.
#'   \item \code{alarmes}: Subgrupos que violaram regras de controle.
#'   \item \code{estavel}: Lógico indicando se o processo está sob controle estatístico.
#'   \item \code{grafico}: Objeto gráfico \code{ggplot2}.
#' }
#' @export
grafico_controle <- function(dados,
                             tipo            = c("imr", "xbar_r", "xbar_s", "p", "np", "c", "u"),
                             subgrupo        = NULL,
                             tamanho_amostra = NULL,
                             regras          = TRUE,
                             grafico         = TRUE) {

  tipo <- match.arg(tipo)

  # ── Despachar por tipo de carta ──────────────────────────────────────────────
  if (tipo == "imr") {
    res <- .carta_imr(dados, regras = regras, grafico = grafico)
  } else if (tipo %in% c("xbar_r", "xbar_s")) {
    res <- .carta_xbar(dados, tipo = tipo, subgrupo = subgrupo, regras = regras, grafico = grafico)
  } else {
    res <- .carta_atributos(dados, tipo = tipo, tamanho_amostra = tamanho_amostra, regras = regras, grafico = grafico)
  }

  return(invisible(res))
}


# ── Carta I-MR (Indivíduos e Amplitude Móvel) ──────────────────────────────────
.carta_imr <- function(x, regras, grafico) {
  if (is.data.frame(x)) x <- x[[1]]
  x <- as.numeric(stats::na.omit(x))
  N <- length(x)
  if (N < 4) stop("A carta I-MR necessita de pelo menos 4 medições.", call. = FALSE)

  # Amplitude móvel de ordem 2
  mr <- c(NA, abs(diff(x)))
  mr_val <- mr[-1]
  mr_bar <- mean(mr_val)
  x_bar  <- mean(x)

  # Limites
  d2 <- .get_const("d2", 2) # 1.128
  D4 <- .get_const("D4", 2) # 3.267
  sigma_est <- mr_bar / d2

  lsc_x <- x_bar + 3 * sigma_est
  lic_x <- x_bar - 3 * sigma_est
  lsc_mr <- D4 * mr_bar
  lic_mr <- 0

  # Regras de alarme
  fora_x <- (x > lsc_x) | (x < lic_x)
  fora_mr <- rep(FALSE, N)
  fora_mr[-1] <- (mr_val > lsc_mr) | (mr_val < lic_mr)

  # Violações adicionais: 8 seguidos do mesmo lado da média
  desloc_x <- rep(FALSE, N)
  if (regras && N >= 8) {
    sinais <- sign(x - x_bar)
    for (i in 8:N) {
      if (all(sinais[(i-7):i] == 1) || all(sinais[(i-7):i] == -1)) {
        desloc_x[(i-7):i] <- TRUE
      }
    }
  }

  status_x <- ifelse(fora_x, "Fora 3s", ifelse(desloc_x, "Deslocamento", "Normal"))
  alarmes_idx <- which(fora_x | desloc_x | fora_mr)
  estavel <- length(alarmes_idx) == 0

  tb_lim <- data.frame(
    Carta       = c("Individuos (I)", "Amplitude Movel (MR)"),
    LIC         = c(lic_x, lic_mr),
    Linha_Media = c(x_bar, mr_bar),
    LSC         = c(lsc_x, lsc_mr),
    Sigma_Curto = c(sigma_est, NA)
  )

  # ── Saída console ────────────────────────────────────────────────────────────
  .print_titulo("CONTROLE ESTAT\u00cdSTICO DE PROCESSO: CARTA I-MR")
  cat(sprintf("  Medidas Analisadas: %d  |  Estimativa de \u03c3 (Curto Prazo): %.4f\n\n", N, sigma_est))

  .print_topico("LIMITES DE CONTROLE ESTAT\u00cdSTICO (3\u03c3)")
  w_lim <- c(24, 16, 16, 16)
  cat(.th_sep(sum(w_lim)), "\n")
  cat("  ", .th_pad("Carta", w_lim[1], "left"), .th_pad("LIC (Inferior)", w_lim[2], "center"),
      .th_pad("Linha M\u00e9dia", w_lim[3], "center"), .th_pad("LSC (Superior)", w_lim[4], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_lim)), "\n")

  cat("  ", .th_pad("Valores Individuais (I)", w_lim[1], "left"),
      .th_pad(formatC(lic_x, format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
      .th_pad(formatC(x_bar, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
      .th_pad(formatC(lsc_x, format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")

  cat("  ", .th_pad("Amplitude M\u00f3vel (MR)", w_lim[1], "left"),
      .th_pad(formatC(lic_mr, format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
      .th_pad(formatC(mr_bar, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
      .th_pad(formatC(lsc_mr, format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_lim)), "\n\n")

  .print_topico("DIAGN\u00d3STICO DE ESTABILIDADE E ALARMES")
  if (estavel) {
    cat("  \u2713 Processo SOB CONTROLE ESTAT\u00cdSTICO. Nenhuma causa especial detectada.\n")
  } else {
    cat(sprintf("  \u26a0 ATEN\u00c7\u00c3O: Foram detectados %d pontos indicando INSTABILIDADE / Causas Especiais:\n", length(alarmes_idx)))
    for (idx in alarmes_idx) {
      motivo <- c()
      if (fora_x[idx]) motivo <- c(motivo, "Ponto al\u00e9m de 3\u03c3 em I")
      if (desloc_x[idx]) motivo <- c(motivo, "8 pontos consecutivos do mesmo lado da m\u00e9dia")
      if (fora_mr[idx]) motivo <- c(motivo, "Amplitude M\u00f3vel excedeu LSC")
      cat(sprintf("    * Ponto %d (Valor = %.3f): %s\n", idx, x[idx], paste(motivo, collapse = " | ")))
    }
  }
  .print_rodape()

  # ── Gráfico ggplot2 ──────────────────────────────────────────────────────────
  p_combined <- NULL
  if (grafico) {
    df_plot_i <- data.frame(
      Ponto  = 1:N,
      Valor  = x,
      Status = ifelse(fora_x | desloc_x, "Alarme", "Normal")
    )
    df_plot_mr <- data.frame(
      Ponto  = 2:N,
      Valor  = mr[-1],
      Status = ifelse(fora_mr[-1], "Alarme", "Normal")
    )

    p1 <- ggplot2::ggplot(df_plot_i, ggplot2::aes(x = Ponto, y = Valor)) +
      ggplot2::geom_line(color = "#2B5C8F", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = Status), size = 2.5) +
      ggplot2::scale_color_manual(values = c("Normal" = "#2B5C8F", "Alarme" = "#D90429")) +
      ggplot2::geom_hline(yintercept = x_bar, color = "#222222", linewidth = 0.9) +
      ggplot2::geom_hline(yintercept = c(lic_x, lsc_x), color = "#D90429", linetype = "dashed", linewidth = 0.8) +
      ggplot2::annotate("text", x = N * 1.01, y = lsc_x, label = sprintf("LSC = %.2f", lsc_x), color = "#D90429", size = 3.5, hjust = 0) +
      ggplot2::annotate("text", x = N * 1.01, y = x_bar, label = sprintf("LM = %.2f", x_bar), color = "#222222", size = 3.5, hjust = 0) +
      ggplot2::annotate("text", x = N * 1.01, y = lic_x, label = sprintf("LIC = %.2f", lic_x), color = "#D90429", size = 3.5, hjust = 0) +
      ggplot2::labs(title = "Carta de Valores Individuais (I)", x = "Amostra / Observa\u00e7\u00e3o", y = "Valor Individual") +
      ggplot2::theme(legend.position = "none")

    p2 <- ggplot2::ggplot(df_plot_mr, ggplot2::aes(x = Ponto, y = Valor)) +
      ggplot2::geom_line(color = "#3470a3", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = Status), size = 2.3) +
      ggplot2::scale_color_manual(values = c("Normal" = "#3470a3", "Alarme" = "#D90429")) +
      ggplot2::geom_hline(yintercept = mr_bar, color = "#222222", linewidth = 0.9) +
      ggplot2::geom_hline(yintercept = lsc_mr, color = "#D90429", linetype = "dashed", linewidth = 0.8) +
      ggplot2::annotate("text", x = N * 1.01, y = lsc_mr, label = sprintf("LSC = %.2f", lsc_mr), color = "#D90429", size = 3.5, hjust = 0) +
      ggplot2::annotate("text", x = N * 1.01, y = mr_bar, label = sprintf("LM = %.2f", mr_bar), color = "#222222", size = 3.5, hjust = 0) +
      ggplot2::labs(title = "Carta de Amplitude M\u00f3vel (MR)", x = "Amostra / Observa\u00e7\u00e3o", y = "Amplitude M\u00f3vel (MR)") +
      ggplot2::theme(legend.position = "none")

    if (exists("meu_tema", mode = "function")) {
      p1 <- p1 + meu_tema(grade = "dupla")
      p2 <- p2 + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p1 <- p1 + tema_estatR(estilo = 1, grade = "dupla")
      p2 <- p2 + tema_estatR(estilo = 1, grade = "dupla")
    } else {
      p1 <- p1 + ggplot2::theme_minimal()
      p2 <- p2 + ggplot2::theme_minimal()
    }

    if (requireNamespace("patchwork", quietly = TRUE)) {
      p_combined <- patchwork::wrap_plots(p1, p2, ncol = 1)
      suppressMessages(print(p_combined))
    } else {
      p_combined <- p1
      suppressMessages(print(p1))
    }
  }

  res <- list(
    limites = tb_lim,
    pontos  = data.frame(Ponto = 1:N, Valor = x, MR = mr, Status = status_x),
    alarmes = alarmes_idx,
    estavel = estavel,
    grafico = p_combined
  )
  class(res) <- c("estatR_ceq", "list")
  return(res)
}


# ── Cartas Xbar-R e Xbar-S ─────────────────────────────────────────────────────
.carta_xbar <- function(dados, tipo, subgrupo, regras, grafico) {
  # Organizar matriz de subgrupos
  if (is.matrix(dados) || is.data.frame(dados)) {
    mat <- as.matrix(dados)
    k <- nrow(mat)
    n <- ncol(mat)
    subg_nomes <- rownames(dados)
    if (is.null(subg_nomes)) subg_nomes <- paste0("Subg_", 1:k)
  } else {
    if (is.null(subgrupo)) stop("Forne\u00e7a o vetor 'subgrupo' para agrupar as medi\u00e7\u00f5es.", call. = FALSE)
    df_temp <- data.frame(y = dados, g = as.factor(subgrupo))
    lst <- split(df_temp$y, df_temp$g)
    k <- length(lst)
    n <- unique(sapply(lst, length))
    if (length(n) > 1) stop("Para cartas Xbar-R/S, os subgrupos devem ter o mesmo tamanho n.", call. = FALSE)
    mat <- do.call(rbind, lst)
    subg_nomes <- names(lst)
  }

  if (n < 2) stop("O tamanho de subgrupo (n) deve ser >= 2.", call. = FALSE)

  medias <- apply(mat, 1, mean)
  x_dbar <- mean(medias)

  if (tipo == "xbar_r") {
    amplitudes <- apply(mat, 1, function(row) max(row) - min(row))
    r_bar <- mean(amplitudes)
    A2 <- .get_const("A2", n)
    D3 <- .get_const("D3", n)
    D4 <- .get_const("D4", n)
    d2 <- .get_const("d2", n)
    sigma_est <- r_bar / d2

    lsc_x <- x_dbar + A2 * r_bar
    lic_x <- x_dbar - A2 * r_bar
    lsc_sec <- D4 * r_bar
    lic_sec <- D3 * r_bar
    sec_vals <- amplitudes
    sec_bar  <- r_bar
    nome_sec <- "Amplitude (R)"
  } else {
    desvios <- apply(mat, 1, stats::sd)
    s_bar <- mean(desvios)
    A3 <- .get_const("A3", n)
    B3 <- .get_const("B3", n)
    B4 <- .get_const("B4", n)
    c4 <- .get_const("c4", n)
    sigma_est <- s_bar / c4

    lsc_x <- x_dbar + A3 * s_bar
    lic_x <- x_dbar - A3 * s_bar
    lsc_sec <- B4 * s_bar
    lic_sec <- B3 * s_bar
    sec_vals <- desvios
    sec_bar  <- s_bar
    nome_sec <- "Desvio-Padr\u00e3o (S)"
  }

  fora_x <- (medias > lsc_x) | (medias < lic_x)
  fora_sec <- (sec_vals > lsc_sec) | (sec_vals < lic_sec)
  alarmes_idx <- which(fora_x | fora_sec)
  estavel <- length(alarmes_idx) == 0

  tb_lim <- data.frame(
    Carta       = c("Media (X-Barra)", nome_sec),
    LIC         = c(lic_x, lic_sec),
    Linha_Media = c(x_dbar, sec_bar),
    LSC         = c(lsc_x, lsc_sec),
    Sigma_Curto = c(sigma_est, NA)
  )

  .print_titulo(sprintf("CONTROLE ESTAT\u00cdSTICO DE PROCESSO: CARTA %s", toupper(tipo)))
  cat(sprintf("  Subgrupos: %d  |  Tamanho por Subgrupo (n): %d  |  \u03c3 Estimado: %.4f\n\n", k, n, sigma_est))

  .print_topico("LIMITES DE CONTROLE ESTAT\u00cdSTICO (3\u03c3)")
  w_lim <- c(24, 16, 16, 16)
  cat(.th_sep(sum(w_lim)), "\n")
  cat("  ", .th_pad("Carta", w_lim[1], "left"), .th_pad("LIC (Inferior)", w_lim[2], "center"),
      .th_pad("Linha M\u00e9dia", w_lim[3], "center"), .th_pad("LSC (Superior)", w_lim[4], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_lim)), "\n")

  cat("  ", .th_pad("M\u00e9dia (X-Barra)", w_lim[1], "left"),
      .th_pad(formatC(lic_x, format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
      .th_pad(formatC(x_dbar, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
      .th_pad(formatC(lsc_x, format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")

  cat("  ", .th_pad(nome_sec, w_lim[1], "left"),
      .th_pad(formatC(lic_sec, format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
      .th_pad(formatC(sec_bar, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
      .th_pad(formatC(lsc_sec, format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_lim)), "\n\n")

  .print_topico("DIAGN\u00d3STICO DE ESTABILIDADE E ALARMES")
  if (estavel) {
    cat("  \u2713 Processo SOB CONTROLE ESTAT\u00cdSTICO. Nenhum subgrupo fora de controle.\n")
  } else {
    cat(sprintf("  \u26a0 ATEN\u00c7\u00c3O: Foram detectados %d subgrupos indicando INSTABILIDADE:\n", length(alarmes_idx)))
    for (idx in alarmes_idx) {
      motivo <- c()
      if (fora_x[idx]) motivo <- c(motivo, "M\u00e9dia excedeu limites 3\u03c3")
      if (fora_sec[idx]) motivo <- c(motivo, sprintf("%s excedeu limites 3\u03c3", nome_sec))
      cat(sprintf("    * Subgrupo %s: %s\n", subg_nomes[idx], paste(motivo, collapse = " | ")))
    }
  }
  .print_rodape()

  p_combined <- NULL
  if (grafico) {
    df_p1 <- data.frame(Ponto = 1:k, Valor = medias, Status = ifelse(fora_x, "Alarme", "Normal"))
    df_p2 <- data.frame(Ponto = 1:k, Valor = sec_vals, Status = ifelse(fora_sec, "Alarme", "Normal"))

    p1 <- ggplot2::ggplot(df_p1, ggplot2::aes(x = Ponto, y = Valor)) +
      ggplot2::geom_line(color = "#2B5C8F", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = Status), size = 2.5) +
      ggplot2::scale_color_manual(values = c("Normal" = "#2B5C8F", "Alarme" = "#D90429")) +
      ggplot2::geom_hline(yintercept = x_dbar, color = "#222222", linewidth = 0.9) +
      ggplot2::geom_hline(yintercept = c(lic_x, lsc_x), color = "#D90429", linetype = "dashed", linewidth = 0.8) +
      ggplot2::labs(title = sprintf("Carta de M\u00e9dias (\u00afX) [n = %d]", n), x = "Subgrupo", y = "M\u00e9dia (\u00afX)") +
      ggplot2::theme(legend.position = "none")

    p2 <- ggplot2::ggplot(df_p2, ggplot2::aes(x = Ponto, y = Valor)) +
      ggplot2::geom_line(color = "#3470a3", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = Status), size = 2.5) +
      ggplot2::scale_color_manual(values = c("Normal" = "#3470a3", "Alarme" = "#D90429")) +
      ggplot2::geom_hline(yintercept = sec_bar, color = "#222222", linewidth = 0.9) +
      ggplot2::geom_hline(yintercept = c(lic_sec, lsc_sec), color = "#D90429", linetype = "dashed", linewidth = 0.8) +
      ggplot2::labs(title = sprintf("Carta de Dispers\u00e3o: %s", nome_sec), x = "Subgrupo", y = nome_sec) +
      ggplot2::theme(legend.position = "none")

    if (exists("meu_tema", mode = "function")) {
      p1 <- p1 + meu_tema(grade = "dupla")
      p2 <- p2 + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p1 <- p1 + tema_estatR(estilo = 1, grade = "dupla")
      p2 <- p2 + tema_estatR(estilo = 1, grade = "dupla")
    } else {
      p1 <- p1 + ggplot2::theme_minimal()
      p2 <- p2 + ggplot2::theme_minimal()
    }

    if (requireNamespace("patchwork", quietly = TRUE)) {
      p_combined <- patchwork::wrap_plots(p1, p2, ncol = 1)
      suppressMessages(print(p_combined))
    } else {
      p_combined <- p1
      suppressMessages(print(p1))
    }
  }

  res <- list(
    limites = tb_lim,
    pontos  = data.frame(Subgrupo = subg_nomes, Media = medias, Dispersao = sec_vals),
    alarmes = alarmes_idx,
    estavel = estavel,
    grafico = p_combined
  )
  class(res) <- c("estatR_ceq", "list")
  return(res)
}


# ── Cartas por Atributos (p, np, c, u) ─────────────────────────────────────────
.carta_atributos <- function(dados, tipo, tamanho_amostra, regras, grafico) {
  k <- length(dados)

  if (tipo %in% c("p", "u")) {
    if (is.null(tamanho_amostra)) stop("Para cartas 'p' e 'u', forne\u00e7a 'tamanho_amostra'.", call. = FALSE)
    n_vec <- if (length(tamanho_amostra) == 1) rep(tamanho_amostra, k) else tamanho_amostra
    if (length(n_vec) != k) stop("'tamanho_amostra' deve ter o mesmo comprimento dos dados.", call. = FALSE)

    if (tipo == "p") {
      # Fração defeituosa
      p_bar <- sum(dados) / sum(n_vec)
      taxa_obs <- dados / n_vec
      lsc <- p_bar + 3 * sqrt(p_bar * (1 - p_bar) / n_vec)
      lic <- pmax(0, p_bar - 3 * sqrt(p_bar * (1 - p_bar) / n_vec))
      lm  <- p_bar
      tit <- "Carta p (Fra\u00e7\u00e3o N\u00e3o-Conforme)"
      ylab <- "Fra\u00e7\u00e3o Defeituosa (p)"
      y_plot <- taxa_obs
    } else {
      # Taxa u
      u_bar <- sum(dados) / sum(n_vec)
      taxa_obs <- dados / n_vec
      lsc <- u_bar + 3 * sqrt(u_bar / n_vec)
      lic <- pmax(0, u_bar - 3 * sqrt(u_bar / n_vec))
      lm  <- u_bar
      tit <- "Carta u (Taxa de Defeitos por Unidade)"
      ylab <- "Defeitos por Unidade (u)"
      y_plot <- taxa_obs
    }
  } else if (tipo == "np") {
    if (is.null(tamanho_amostra) || length(tamanho_amostra) != 1) {
      stop("A carta 'np' exige tamanho de amostra constante 'tamanho_amostra'.", call. = FALSE)
    }
    n <- tamanho_amostra
    p_bar <- sum(dados) / (k * n)
    lm <- n * p_bar
    lsc <- lm + 3 * sqrt(n * p_bar * (1 - p_bar))
    lic <- max(0, lm - 3 * sqrt(n * p_bar * (1 - p_bar)))
    tit <- sprintf("Carta np (N\u00famero de Defeituosos) [n = %d]", n)
    ylab <- "Contagem de N\u00e3o-Conformes (np)"
    y_plot <- dados
  } else {
    # tipo == "c"
    c_bar <- mean(dados)
    lm <- c_bar
    lsc <- c_bar + 3 * sqrt(c_bar)
    lic <- max(0, c_bar - 3 * sqrt(c_bar))
    tit <- "Carta c (Contagem Total de Defeitos)"
    ylab <- "N\u00famero de Defeitos (c)"
    y_plot <- dados
  }

  fora <- (y_plot > lsc) | (y_plot < lic)
  alarmes_idx <- which(fora)
  estavel <- length(alarmes_idx) == 0

  .print_titulo(sprintf("CONTROLE ESTAT\u00cdSTICO DE PROCESSO: %s", toupper(tit)))
  cat(sprintf("  Amostras Inspecionadas: %d  |  Linha M\u00e9dia: %.4f\n\n", k, mean(lm)))

  .print_topico("LIMITES DE CONTROLE ESTAT\u00cdSTICO (3\u03c3)")
  cat(sprintf("  LIC: %.4f  |  Linha M\u00e9dia: %.4f  |  LSC: %.4f\n\n", min(lic), mean(lm), max(lsc)))

  .print_topico("DIAGN\u00d3STICO DE ESTABILIDADE E ALARMES")
  if (estavel) {
    cat("  \u2713 Processo SOB CONTROLE ESTAT\u00cdSTICO. Nenhuma viola\u00e7\u00e3o detectada.\n")
  } else {
    cat(sprintf("  \u26a0 ATEN\u00c7\u00c3O: %d amostras violaram os limites de controle:\n", length(alarmes_idx)))
    for (idx in alarmes_idx) {
      cat(sprintf("    * Amostra %d (Valor = %.3f) fora dos limites [LIC = %.3f, LSC = %.3f]\n",
                  idx, y_plot[idx], lic[min(idx, length(lic))], lsc[min(idx, length(lsc))]))
    }
  }
  .print_rodape()

  p_out <- NULL
  if (grafico) {
    df_p <- data.frame(
      Amostra = 1:k,
      Valor   = y_plot,
      LSC     = lsc,
      LM      = lm,
      LIC     = lic,
      Status  = ifelse(fora, "Alarme", "Normal")
    )

    p_out <- ggplot2::ggplot(df_p, ggplot2::aes(x = Amostra, y = Valor)) +
      ggplot2::geom_line(color = "#2B5C8F", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = Status), size = 2.5) +
      ggplot2::scale_color_manual(values = c("Normal" = "#2B5C8F", "Alarme" = "#D90429")) +
      ggplot2::geom_line(ggplot2::aes(y = LM), color = "#222222", linewidth = 0.9) +
      ggplot2::geom_line(ggplot2::aes(y = LSC), color = "#D90429", linetype = "dashed", linewidth = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = LIC), color = "#D90429", linetype = "dashed", linewidth = 0.8) +
      ggplot2::labs(title = tit, x = "Amostra / Subgrupo", y = ylab) +
      ggplot2::theme(legend.position = "none")

    if (exists("meu_tema", mode = "function")) {
      p_out <- p_out + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p_out <- p_out + tema_estatR(estilo = 1, grade = "dupla")
    } else {
      p_out <- p_out + ggplot2::theme_minimal()
    }
    suppressMessages(print(p_out))
  }

  res <- list(
    limites = data.frame(Amostra = 1:k, LIC = lic, LM = lm, LSC = lsc),
    pontos  = data.frame(Amostra = 1:k, Valor = y_plot, Status = ifelse(fora, "Alarme", "Normal")),
    alarmes = alarmes_idx,
    estavel = estavel,
    grafico = p_out
  )
  class(res) <- c("estatR_ceq", "list")
  return(res)
}


# ── Capacidade de Processo (Cp, Cpk, Pp, Ppk, PPM, Nível Sigma) ───────────────

#' @title Análise de Capacidade do Processo (Six Sigma & Índices Cp / Cpk)
#' @description Avalia a capacidade de um processo manufatureiro ou operacional de atender
#'   aos limites de especificação da engenharia (LIE e LSE). Calcula índices de potencial
#'   e desempenho de curto e longo prazo (\eqn{C_p}, \eqn{C_{pk}}, \eqn{C_{pm}}, \eqn{P_p},
#'   \eqn{P_{pk}}), PPM esperado e nível Sigma do processo.
#'
#' @param dados Vetor numérico com as medições de processo.
#' @param lie Limite Inferior de Especificação (LSL / LIE). Opcional se houver LSE.
#' @param lse Limite Superior de Especificação (USL / LSE). Opcional se houver LIE.
#' @param alvo Valor nominal ou alvo ideal do processo (opcional, para índice \eqn{C_{pm}}).
#' @param confianca Nível de confiança para o intervalo (padrão 0.95).
#' @param grafico Lógico. Se \code{TRUE} (padrão), plota o histograma com curva normal teórica,
#'   linhas de especificação e painel de resumo.
#'
#' @return Retorna um objeto da classe \code{estatR_capacidade} contendo:
#' \itemize{
#'   \item \code{indices}: Tabela com os índices calculados (Cp, Cpk, Cpl, Cpu, Cpm, Pp, Ppk).
#'   \item \code{ppm}: Tabela com partes por milhão (PPM) esperado e observado.
#'   \item \code{sigma_nivel}: Nível Sigma do processo (\eqn{Z}).
#'   \item \code{grafico}: Objeto gráfico \code{ggplot2}.
#' }
#' @export
capacidade_processo <- function(dados,
                                lie       = NULL,
                                lse       = NULL,
                                alvo      = NULL,
                                confianca = 0.95,
                                grafico   = TRUE) {

  if (is.data.frame(dados)) dados <- dados[[1]]
  x <- as.numeric(stats::na.omit(dados))
  n <- length(x)

  if (n < 10) stop("A análise de capacidade exige pelo menos 10 medições.", call. = FALSE)
  if (is.null(lie) && is.null(lse)) {
    stop("Informe pelo menos um limite de especificação ('lie' ou 'lse').", call. = FALSE)
  }
  if (!is.null(lie) && !is.null(lse) && lie >= lse) {
    stop("O 'lie' deve ser estritamente menor que o 'lse'.", call. = FALSE)
  }

  media_obs <- mean(x)
  s_longo   <- stats::sd(x)

  # Estimativa de curto prazo via amplitude móvel média
  mr <- abs(diff(x))
  d2 <- .get_const("d2", 2)
  s_curto <- mean(mr) / d2

  # ── Cálculo dos Índices ──────────────────────────────────────────────────────
  cp <- NA_real_; cpl <- NA_real_; cpu <- NA_real_; cpk <- NA_real_
  pp <- NA_real_; ppl <- NA_real_; ppu <- NA_real_; ppk <- NA_real_
  cpm <- NA_real_

  if (!is.null(lie) && !is.null(lse)) {
    cp  <- (lse - lie) / (6 * s_curto)
    pp  <- (lse - lie) / (6 * s_longo)
  }
  if (!is.null(lie)) {
    cpl <- (media_obs - lie) / (3 * s_curto)
    ppl <- (media_obs - lie) / (3 * s_longo)
  }
  if (!is.null(lse)) {
    cpu <- (lse - media_obs) / (3 * s_curto)
    ppu <- (lse - media_obs) / (3 * s_longo)
  }

  if (!is.null(lie) && !is.null(lse)) {
    cpk <- min(cpl, cpu)
    ppk <- min(ppl, ppu)
  } else if (!is.null(lie)) {
    cpk <- cpl; ppk <- ppl
  } else {
    cpk <- cpu; ppk <- ppu
  }

  if (!is.null(alvo) && !is.null(lie) && !is.null(lse)) {
    tau <- sqrt(s_curto^2 + (media_obs - alvo)^2)
    cpm <- (lse - lie) / (6 * tau)
  }

  # ── PPM Esperado e Observado ─────────────────────────────────────────────────
  p_inf_esp <- if (!is.null(lie)) stats::pnorm(lie, mean = media_obs, sd = s_curto) else 0
  p_sup_esp <- if (!is.null(lse)) (1 - stats::pnorm(lse, mean = media_obs, sd = s_curto)) else 0
  ppm_esp_total <- (p_inf_esp + p_sup_esp) * 1e6

  n_inf_obs <- if (!is.null(lie)) sum(x < lie) else 0
  n_sup_obs <- if (!is.null(lse)) sum(x > lse) else 0
  ppm_obs_total <- ((n_inf_obs + n_sup_obs) / n) * 1e6

  # Nível Sigma do processo (com 1.5 sigma shift habitual)
  z_bench <- 3 * cpk
  nivel_sigma <- z_bench + 1.5

  # Classificação qualitativa
  classif <- if (!is.na(cpk) && cpk >= 1.67) {
    "Excelente (Seis Sigma / Capacidade de Classe Mundial)"
  } else if (!is.na(cpk) && cpk >= 1.33) {
    "Adequado / Capaz (Atende aos requisitos industriais padrão)"
  } else if (!is.na(cpk) && cpk >= 1.00) {
    "Marginalmente Capaz (Risco elevado de defeituosos caso haja desvio)"
  } else {
    "Incapaz (Processo gerando produtos fora de especificação)"
  }

  tb_indices <- data.frame(
    Indice    = c("Cp (Potencial)", "Cpk (Capacidade Real)", "Cpm (Taguchi)", "Pp (Desempenho)", "Ppk (Desempenho Real)"),
    Estimativa = c(cp, cpk, cpm, pp, ppk),
    stringsAsFactors = FALSE
  )
  tb_indices <- tb_indices[!is.na(tb_indices$Estimativa), ]

  # ── Saída console ────────────────────────────────────────────────────────────
  .print_titulo("AN\u00c1LISE DE CAPACIDADE DO PROCESSO (SIX SIGMA)")

  cat(sprintf("  Amostras (n): %d  |  M\u00e9dia: %.4f  |  Alvo (Target): %s\n",
              n, media_obs, if (!is.null(alvo)) formatC(alvo, format = "f", digits = 2) else "N/A"))
  cat(sprintf("  Especifica\u00e7\u00f5es: LIE = %s  |  LSE = %s\n",
              if (!is.null(lie)) formatC(lie, format = "f", digits = 2) else "-Inf",
              if (!is.null(lse)) formatC(lse, format = "f", digits = 2) else "+Inf"))
  cat(sprintf("  Dispers\u00e3o: \u03c3 Curto Prazo = %.4f  |  s Longo Prazo = %.4f\n\n", s_curto, s_longo))

  .print_topico("\u00cdNDICES DE CAPACIDADE E DESEMPENHO")
  w_cap <- c(28, 16)
  cat(.th_sep(sum(w_cap)), "\n")
  cat("  ", .th_pad("\u00cdndice", w_cap[1], "left"), .th_pad("Valor", w_cap[2], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_cap)), "\n")

  for (i in 1:nrow(tb_indices)) {
    cat("  ", .th_pad(tb_indices$Indice[i], w_cap[1], "left"),
        .th_pad(formatC(tb_indices$Estimativa[i], format = "f", digits = 3, decimal.mark = ","), w_cap[2], "center"), "\n", sep = "")
  }
  cat(.th_sep(sum(w_cap)), "\n\n")

  .print_topico("DESEMPENHO DO PROCESSO (PPM E N\u00cdVEL SIGMA)")
  cat(sprintf("  PPM Esperado (Te\u00f3rico):   %.1f defeitos por milh\u00e3o\n", ppm_esp_total))
  cat(sprintf("  PPM Observado na Amostra: %.1f defeitos por milh\u00e3o (%d pe\u00e7as n\u00e3o-conformes)\n",
              ppm_obs_total, n_inf_obs + n_sup_obs))
  cat(sprintf("  N\u00edvel Sigma (\u03c3-level):      %.2f \u03c3\n\n", nivel_sigma))

  .print_topico("VEREDITO")
  cat(sprintf("  Diagn\u00f3stico: %s\n", classif))
  .print_rodape()

  # ── Gráfico ggplot2 ──────────────────────────────────────────────────────────
  p_cap <- NULL
  if (grafico) {
    x_min <- min(c(x, if (!is.null(lie)) lie else min(x)))
    x_max <- max(c(x, if (!is.null(lse)) lse else max(x)))
    rng <- x_max - x_min
    seq_x <- seq(x_min - 0.1 * rng, x_max + 0.1 * rng, length.out = 400)
    df_curva <- data.frame(x = seq_x, y = stats::dnorm(seq_x, mean = media_obs, sd = s_curto))

    p_cap <- ggplot2::ggplot(data.frame(x = x), ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 20, fill = "#85b4d4", color = "#555555", alpha = 0.6) +
      ggplot2::geom_line(data = df_curva, ggplot2::aes(x = x, y = y), color = "#2B5C8F", linewidth = 1.2)

    if (!is.null(lie)) {
      p_cap <- p_cap +
        ggplot2::geom_vline(xintercept = lie, color = "#D90429", linewidth = 1.1, linetype = "dashed") +
        ggplot2::annotate("text", x = lie, y = max(df_curva$y) * 0.9, label = sprintf(" LIE = %.2f", lie), color = "#D90429", fontface = "bold", hjust = 1.05)
    }
    if (!is.null(lse)) {
      p_cap <- p_cap +
        ggplot2::geom_vline(xintercept = lse, color = "#D90429", linewidth = 1.1, linetype = "dashed") +
        ggplot2::annotate("text", x = lse, y = max(df_curva$y) * 0.9, label = sprintf("LSE = %.2f ", lse), color = "#D90429", fontface = "bold", hjust = -0.05)
    }
    if (!is.null(alvo)) {
      p_cap <- p_cap +
        ggplot2::geom_vline(xintercept = alvo, color = "#008000", linewidth = 1.1, linetype = "dotted") +
        ggplot2::annotate("text", x = alvo, y = max(df_curva$y) * 1.02, label = sprintf("Alvo = %.2f", alvo), color = "#008000", fontface = "bold", hjust = 0.5)
    }

    p_cap <- p_cap +
      ggplot2::labs(
        title    = "Capacidade do Processo",
        subtitle = sprintf("Cp = %.2f  |  Cpk = %.2f  |  PPM Esp = %.0f  |  N\u00edvel = %.2f\u03c3",
                           if (!is.na(cp)) cp else cpk, cpk, ppm_esp_total, nivel_sigma),
        x = "Caracter\u00edstica de Qualidade",
        y = "Densidade"
      )

    if (exists("meu_tema", mode = "function")) {
      p_cap <- p_cap + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p_cap <- p_cap + tema_estatR(estilo = 1, grade = "dupla")
    } else {
      p_cap <- p_cap + ggplot2::theme_minimal()
    }

    suppressMessages(print(p_cap))
  }

  res <- list(
    indices     = tb_indices,
    ppm         = data.frame(Tipo = c("Esperado", "Observado"), PPM = c(ppm_esp_total, ppm_obs_total)),
    sigma_nivel = nivel_sigma,
    grafico     = p_cap
  )
  class(res) <- c("estatR_capacidade", "list")
  return(invisible(res))
}
