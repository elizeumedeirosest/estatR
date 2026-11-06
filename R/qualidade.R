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

# ── Helper de renderização gráfica unificada para Cartas de Shewhart ─────────
.montar_grafico_shewhart <- function(pontos_x, valores, lc, lic, lsc, titulo, xlab, ylab, fora) {
  df <- data.frame(
    Ponto = pontos_x,
    Valor = valores,
    Fora  = fora
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Ponto, y = Valor))

  if (length(lc) == 1) {
    p <- p +
      ggplot2::geom_hline(ggplot2::aes(yintercept = lc, color = "Linha Central", linetype = "Linha Central"), linewidth = 0.9) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = lsc, color = "Limite de Controle", linetype = "Limite de Controle"), linewidth = 0.9) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = lic, color = "Limite de Controle", linetype = "Limite de Controle"), linewidth = 0.9)
  } else {
    df$LC  <- lc
    df$LSC <- lsc
    df$LIC <- lic
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = LC, color = "Linha Central", linetype = "Linha Central"), linewidth = 0.9) +
      ggplot2::geom_line(ggplot2::aes(y = LSC, color = "Limite de Controle", linetype = "Limite de Controle"), linewidth = 0.9) +
      ggplot2::geom_line(ggplot2::aes(y = LIC, color = "Limite de Controle", linetype = "Limite de Controle"), linewidth = 0.9)
  }

  # Linha e pontos de dados
  p <- p +
    ggplot2::geom_line(ggplot2::aes(color = "Dados"), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(color = "Dados"), size = 3.0)

  # Destaque de pontos fora de controle (anel/halo circular)
  df_fora <- subset(df, Fora)
  tem_fora <- nrow(df_fora) > 0
  if (tem_fora) {
    p <- p +
      ggplot2::geom_point(data = df_fora, ggplot2::aes(color = "Valor Isolado"), shape = 1, size = 6.8, stroke = 1.3) +
      ggplot2::geom_point(data = df_fora, color = "#D90429", size = 3.8)
  }

  cores <- c(
    "Dados"              = "#1E5B94",
    "Valor Isolado"      = "#D90429",
    "Limite de Controle" = "#D90429",
    "Linha Central"      = "#222222"
  )

  linetypes <- c(
    "Dados"              = "solid",
    "Valor Isolado"      = "blank",
    "Limite de Controle" = "dashed",
    "Linha Central"      = "solid"
  )

  niveis_legenda <- if (tem_fora) {
    c("Dados", "Valor Isolado", "Limite de Controle", "Linha Central")
  } else {
    c("Dados", "Limite de Controle", "Linha Central")
  }

  shapes_override <- if (tem_fora) c(16, 1, NA, NA) else c(16, NA, NA)
  sizes_override  <- if (tem_fora) c(3.0, 4.5, NA, NA) else c(3.0, NA, NA)
  lt_override     <- if (tem_fora) c("solid", "blank", "dashed", "solid") else c("solid", "dashed", "solid")
  lw_override     <- if (tem_fora) c(0.8, 0, 0.9, 0.9) else c(0.8, 0.9, 0.9)

  p <- p +
    ggplot2::scale_color_manual(name = "", values = cores, breaks = niveis_legenda) +
    ggplot2::scale_linetype_manual(name = "", values = linetypes, breaks = niveis_legenda) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(
        shape = shapes_override,
        size = sizes_override,
        linetype = lt_override,
        linewidth = lw_override
      ))
    ) +
    ggplot2::labs(title = titulo, x = xlab, y = ylab)

  if (exists("meu_tema", mode = "function")) {
    p <- p + meu_tema(grade = "dupla")
  } else if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 1, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p <- p + ggplot2::theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box = "horizontal",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 9.5, color = "#333333"),
    plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0, margin = ggplot2::margin(b = 6))
  )

  p
}


#' @title Gráficos de Controle Estatístico de Processo (Cartas de Shewhart)
#' @description Constrói cartas de controle de Shewhart para monitoramento de estabilidade
#'   e controle de processos contínuos e por atributos, com limites 3\eqn{\sigma}, detecção
#'   de causas especiais e visualização gráfica padronizada com legenda no topo.
#'
#' @param dados Vetor numérico ou matriz/data.frame com as medições de processo.
#' @param tipo Tipo de carta de controle:
#'   \itemize{
#'     \item \code{"imr"}: Carta de Valores Individuais e Amplitude Móvel (I-MR conjunta).
#'     \item \code{"x"} ou \code{"i"}: Carta de Valores Individuais (X / I isolada).
#'     \item \code{"xbar_r"}: Carta de Média e Amplitude (\eqn{\bar{X}}-R conjunta).
#'     \item \code{"xbar_s"}: Carta de Média e Desvio-Padrão (\eqn{\bar{X}}-S conjunta).
#'     \item \code{"xbar"}: Carta de Médias (\eqn{\bar{X}} isolada).
#'     \item \code{"r"}: Carta de Amplitudes (R isolada).
#'     \item \code{"s"}: Carta de Desvios-Padrão (S isolada).
#'     \item \code{"p"}: Proporção de não-conformes (tamanho amostral variável ou constante).
#'     \item \code{"np"}: Número de unidades defeituosas (tamanho amostral constante).
#'     \item \code{"c"}: Contagem total de defeitos por amostra (tamanho constante).
#'     \item \code{"u"}: Taxa de defeitos por unidade de inspeção (tamanho variável).
#'   }
#' @param subgrupo Vetor indicando a identificação dos subgrupos racionais (opcional).
#' @param tamanho_amostra Tamanho dos subgrupos \eqn{n} para cartas por atributos (\code{p}, \code{np}, \code{u}).
#' @param regras Lógico. Se \code{TRUE} (padrão), identifica e destaca pontos além dos limites 3\eqn{\sigma}.
#' @param grafico Lógico. Se \code{TRUE} (padrão), renderiza os gráficos de controle no padrão do estatR.
#'
#' @return Retorna um objeto da classe \code{estatR_ceq} contendo:
#' \itemize{
#'   \item \code{limites}: Tabela de Limites Estatísticos (LSC, Linha Média, LIC).
#'   \item \code{pontos}: Tabela com os valores observados de cada subgrupo e status.
#'   \item \code{alarmes}: Subgrupos que violaram os limites 3\eqn{\sigma}.
#'   \item \code{estavel}: Lógico indicando se o processo está sob controle estatístico.
#'   \item \code{grafico}: Objeto gráfico \code{ggplot2}.
#' }
#' @export
grafico_controle <- function(dados,
                             tipo            = c("xbar_r", "xbar_s", "imr", "x", "xbar", "r", "s", "p", "np", "c", "u", "i"),
                             subgrupo        = NULL,
                             tamanho_amostra = NULL,
                             regras          = TRUE,
                             grafico         = TRUE) {

  tipo <- match.arg(tipo)
  if (tipo == "i") tipo <- "x"

  # ── Despachar por tipo de carta ──────────────────────────────────────────────
  if (tipo %in% c("imr", "x")) {
    res <- .carta_imr_internal(dados, isolada = (tipo == "x"), regras = regras, grafico = grafico)
  } else if (tipo %in% c("xbar_r", "xbar_s", "xbar", "r", "s")) {
    res <- .carta_xbar_internal(dados, tipo = tipo, subgrupo = subgrupo, regras = regras, grafico = grafico)
  } else {
    res <- .carta_atributos(dados, tipo = tipo, tamanho_amostra = tamanho_amostra, regras = regras, grafico = grafico)
  }

  return(invisible(res))
}


# ── Carta I-MR e Carta X Isolada ───────────────────────────────────────────────
.carta_imr_internal <- function(x, isolada = FALSE, regras = TRUE, grafico = TRUE) {
  if (is.data.frame(x)) x <- x[[1]]
  x <- as.numeric(stats::na.omit(x))
  N <- length(x)
  if (N < 4) stop("A carta necessita de pelo menos 4 medições.", call. = FALSE)

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

  fora_x <- (x > lsc_x) | (x < lic_x)
  fora_mr <- rep(FALSE, N)
  fora_mr[-1] <- (mr_val > lsc_mr) | (mr_val < lic_mr)

  if (isolada) {
    alarmes_idx <- which(fora_x)
    estavel <- length(alarmes_idx) == 0
    tb_lim <- data.frame(
      Carta       = "Valores Individuais (X)",
      LIC         = lic_x,
      Linha_Media = x_bar,
      LSC         = lsc_x,
      Sigma_Curto = sigma_est
    )
  } else {
    alarmes_idx <- which(fora_x | fora_mr)
    estavel <- length(alarmes_idx) == 0
    tb_lim <- data.frame(
      Carta       = c("Indivíduos (I)", "Amplitude Móvel (MR)"),
      LIC         = c(lic_x, lic_mr),
      Linha_Media = c(x_bar, mr_bar),
      LSC         = c(lsc_x, lsc_mr),
      Sigma_Curto = c(sigma_est, NA)
    )
  }

  # ── Saída console ────────────────────────────────────────────────────────────
  titulo_bloco <- if (isolada) "CONTROLE ESTATÍSTICO DE PROCESSO: CARTA X (VALORES INDIVIDUAIS)" else "CONTROLE ESTATÍSTICO DE PROCESSO: CARTA I-MR"
  .print_titulo(titulo_bloco)
  cat(sprintf("  Medidas Analisadas: %d  |  σ̂: %s\n\n", N, formatC(sigma_est, format = "f", digits = 4, decimal.mark = ",")))

  .print_topico("LIMITES DE CONTROLE ESTATÍSTICO (3σ)")
  w_lim <- c(24, 16, 16, 16)
  cat(.th_sep(sum(w_lim)), "\n")
  cat("  ", .th_pad("Carta", w_lim[1], "left"), .th_pad("LIC (Inferior)", w_lim[2], "center"),
      .th_pad("Linha Média", w_lim[3], "center"), .th_pad("LSC (Superior)", w_lim[4], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_lim)), "\n")

  if (isolada) {
    cat("  ", .th_pad("Valores Individuais (X)", w_lim[1], "left"),
        .th_pad(formatC(lic_x, format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
        .th_pad(formatC(x_bar, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
        .th_pad(formatC(lsc_x, format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")
  } else {
    cat("  ", .th_pad("Valores Individuais (I)", w_lim[1], "left"),
        .th_pad(formatC(lic_x, format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
        .th_pad(formatC(x_bar, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
        .th_pad(formatC(lsc_x, format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")
    cat("  ", .th_pad("Amplitude Móvel (MR)", w_lim[1], "left"),
        .th_pad(formatC(lic_mr, format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
        .th_pad(formatC(mr_bar, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
        .th_pad(formatC(lsc_mr, format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")
  }
  cat(.th_sep(sum(w_lim)), "\n\n")

  .print_topico("INTERPRETAÇÃO E ESTABILIDADE DO PROCESSO")
  if (estavel) {
    cat("  Processo sob controle estatístico (estável): todos os pontos estão contidos dentro dos limites 3σ.\n\n")
  } else {
    cat(sprintf("  Processo fora de controle estatístico (instável): %d ponto(s) além dos limites 3σ de controle:\n", length(alarmes_idx)))
    for (idx in alarmes_idx) {
      motivos <- c()
      if (fora_x[idx]) motivos <- c(motivos, sprintf("Valor (%.3f) excedeu limites de X", x[idx]))
      if (!isolada && idx > 1 && fora_mr[idx]) motivos <- c(motivos, sprintf("Amplitude Móvel (%.3f) excedeu LSC", mr[idx]))
      cat(sprintf("    * Ponto %d: %s\n", idx, paste(motivos, collapse = " | ")))
    }
    cat("\n")
  }
  .print_rodape()

  # ── Gráficos ggplot2 ─────────────────────────────────────────────────────────
  p_combined <- NULL
  if (grafico) {
    p1 <- .montar_grafico_shewhart(
      pontos_x = 1:N,
      valores  = x,
      lc       = x_bar,
      lic      = lic_x,
      lsc      = lsc_x,
      titulo   = if (isolada) "Carta de Valores Individuais (X)" else "Carta de Valores Individuais (I)",
      xlab     = "Amostra / Observação",
      ylab     = "Valor Individual",
      fora     = fora_x
    )

    if (isolada) {
      p_combined <- p1
      suppressMessages(print(p1))
    } else {
      p2 <- .montar_grafico_shewhart(
        pontos_x = 2:N,
        valores  = mr_val,
        lc       = mr_bar,
        lic      = lic_mr,
        lsc      = lsc_mr,
        titulo   = "Carta de Amplitude Móvel (MR)",
        xlab     = "Amostra / Observação",
        ylab     = "Amplitude Móvel (MR)",
        fora     = fora_mr[-1]
      )

      if (requireNamespace("patchwork", quietly = TRUE)) {
        p_combined <- patchwork::wrap_plots(p1, p2, ncol = 1)
        suppressMessages(print(p_combined))
      } else {
        p_combined <- p1
        suppressMessages(print(p1))
      }
    }
  }

  res <- list(
    limites = tb_lim,
    pontos  = data.frame(Ponto = 1:N, Valor = x, MR = mr, Status = ifelse(fora_x, "Alarme", "Normal")),
    alarmes = alarmes_idx,
    estavel = estavel,
    grafico = p_combined
  )
  class(res) <- c("estatR_ceq", "list")
  return(res)
}


# ── Cartas Xbar-R, Xbar-S e Cartas Isoladas (Xbar, R, S) ──────────────────────
.carta_xbar_internal <- function(dados, tipo, subgrupo, regras, grafico) {
  # Organizar matriz de subgrupos
  if (is.matrix(dados) || is.data.frame(dados)) {
    mat <- as.matrix(dados)
    k <- nrow(mat)
    n <- ncol(mat)
    subg_nomes <- rownames(dados)
    if (is.null(subg_nomes)) subg_nomes <- paste0("Subg_", 1:k)
  } else {
    if (is.null(subgrupo)) stop("Forneça o vetor 'subgrupo' para agrupar as medições.", call. = FALSE)
    df_temp <- data.frame(y = dados, g = as.factor(subgrupo))
    lst <- split(df_temp$y, df_temp$g)
    k <- length(lst)
    n <- unique(sapply(lst, length))
    if (length(n) > 1) stop("Para cartas baseadas em subgrupos, todos os subgrupos devem ter tamanho n constante.", call. = FALSE)
    mat <- do.call(rbind, lst)
    subg_nomes <- names(lst)
  }

  if (n < 2) stop("O tamanho de subgrupo (n) deve ser >= 2.", call. = FALSE)

  medias <- apply(mat, 1, mean)
  x_dbar <- mean(medias)

  usa_desvio <- tipo %in% c("xbar_s", "s")

  if (!usa_desvio) {
    # Baseado em Amplitude (R)
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
    # Baseado em Desvio-Padrão (S)
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
    nome_sec <- "Desvio-Padrão (S)"
  }

  fora_x   <- (medias > lsc_x) | (medias < lic_x)
  fora_sec <- (sec_vals > lsc_sec) | (sec_vals < lic_sec)

  # Determinar modo de exibição (isolada vs conjunta)
  if (tipo == "xbar") {
    alarmes_idx <- which(fora_x)
    estavel <- length(alarmes_idx) == 0
    tb_lim <- data.frame(
      Carta       = "Média (X̄)",
      LIC         = lic_x,
      Linha_Media = x_dbar,
      LSC         = lsc_x
    )
    titulo_console <- "CONTROLE ESTATÍSTICO DE PROCESSO: CARTA X̄ (MÉDIAS)"
  } else if (tipo == "r") {
    alarmes_idx <- which(fora_sec)
    estavel <- length(alarmes_idx) == 0
    tb_lim <- data.frame(
      Carta       = "Amplitude (R)",
      LIC         = lic_sec,
      Linha_Media = sec_bar,
      LSC         = lsc_sec
    )
    titulo_console <- "CONTROLE ESTATÍSTICO DE PROCESSO: CARTA R (AMPLITUDE)"
  } else if (tipo == "s") {
    alarmes_idx <- which(fora_sec)
    estavel <- length(alarmes_idx) == 0
    tb_lim <- data.frame(
      Carta       = "Desvio-Padrão (S)",
      LIC         = lic_sec,
      Linha_Media = sec_bar,
      LSC         = lsc_sec
    )
    titulo_console <- "CONTROLE ESTATÍSTICO DE PROCESSO: CARTA S (DESVIO-PADRÃO)"
  } else {
    # Conjunta (xbar_r ou xbar_s)
    alarmes_idx <- which(fora_x | fora_sec)
    estavel <- length(alarmes_idx) == 0
    tb_lim <- data.frame(
      Carta       = c("Média (X̄)", nome_sec),
      LIC         = c(lic_x, lic_sec),
      Linha_Media = c(x_dbar, sec_bar),
      LSC         = c(lsc_x, lsc_sec)
    )
    titulo_console <- if (tipo == "xbar_r") "CONTROLE ESTATÍSTICO DE PROCESSO: CARTA X̄-R (MÉDIAS E AMPLITUDES)" else "CONTROLE ESTATÍSTICO DE PROCESSO: CARTA X̄-S (MÉDIAS E DESVIOS)"
  }

  # ── Saída console ────────────────────────────────────────────────────────────
  .print_titulo(titulo_console)
  cat(sprintf("  Subgrupos: %d  |  Tamanho por Subgrupo (n): %d  |  σ̂: %s\n\n",
              k, n, formatC(sigma_est, format = "f", digits = 4, decimal.mark = ",")))

  .print_topico("LIMITES DE CONTROLE ESTATÍSTICO (3σ)")
  w_lim <- c(24, 16, 16, 16)
  cat(.th_sep(sum(w_lim)), "\n")
  cat("  ", .th_pad("Carta", w_lim[1], "left"), .th_pad("LIC (Inferior)", w_lim[2], "center"),
      .th_pad("Linha Média", w_lim[3], "center"), .th_pad("LSC (Superior)", w_lim[4], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_lim)), "\n")

  for (i in 1:nrow(tb_lim)) {
    cat("  ", .th_pad(tb_lim$Carta[i], w_lim[1], "left"),
        .th_pad(formatC(tb_lim$LIC[i], format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
        .th_pad(formatC(tb_lim$Linha_Media[i], format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
        .th_pad(formatC(tb_lim$LSC[i], format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")
  }
  cat(.th_sep(sum(w_lim)), "\n\n")

  .print_topico("INTERPRETAÇÃO E ESTABILIDADE DO PROCESSO")
  if (estavel) {
    cat("  Processo sob controle estatístico (estável): nenhum subgrupo violou os limites 3σ.\n\n")
  } else {
    cat(sprintf("  Processo fora de controle estatístico (instável): %d subgrupo(s) violaram os limites 3σ:\n", length(alarmes_idx)))
    for (idx in alarmes_idx) {
      motivo <- c()
      if (fora_x[idx] && tipo != "r" && tipo != "s") motivo <- c(motivo, "Média excedeu limites 3σ")
      if (fora_sec[idx] && tipo != "xbar") motivo <- c(motivo, sprintf("%s excedeu limites 3σ", nome_sec))
      cat(sprintf("    * Subgrupo %s: %s\n", subg_nomes[idx], paste(motivo, collapse = " | ")))
    }
    cat("\n")
  }
  .print_rodape()

  p_combined <- NULL
  if (grafico) {
    p_x <- .montar_grafico_shewhart(
      pontos_x = 1:k,
      valores  = medias,
      lc       = x_dbar,
      lic      = lic_x,
      lsc      = lsc_x,
      titulo   = bquote("Carta de Médias (" ~ bar(X) ~ ")"),
      xlab     = "Subgrupo",
      ylab     = expression(bar(X)),
      fora     = fora_x
    )

    p_sec <- .montar_grafico_shewhart(
      pontos_x = 1:k,
      valores  = sec_vals,
      lc       = sec_bar,
      lic      = lic_sec,
      lsc      = lsc_sec,
      titulo   = sprintf("Carta de Dispersão: %s", nome_sec),
      xlab     = "Subgrupo",
      ylab     = nome_sec,
      fora     = fora_sec
    )

    if (tipo == "xbar") {
      p_combined <- p_x
      suppressMessages(print(p_x))
    } else if (tipo %in% c("r", "s")) {
      p_combined <- p_sec
      suppressMessages(print(p_sec))
    } else {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        p_combined <- patchwork::wrap_plots(p_x, p_sec, ncol = 1)
        suppressMessages(print(p_combined))
      } else {
        p_combined <- p_x
        suppressMessages(print(p_x))
      }
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
    if (is.null(tamanho_amostra)) stop("Para cartas 'p' e 'u', forneça 'tamanho_amostra'.", call. = FALSE)
    n_vec <- if (length(tamanho_amostra) == 1) rep(tamanho_amostra, k) else tamanho_amostra
    if (length(n_vec) != k) stop("'tamanho_amostra' deve ter o mesmo comprimento dos dados.", call. = FALSE)

    if (tipo == "p") {
      p_bar <- sum(dados) / sum(n_vec)
      taxa_obs <- dados / n_vec
      lsc <- p_bar + 3 * sqrt(p_bar * (1 - p_bar) / n_vec)
      lic <- pmax(0, p_bar - 3 * sqrt(p_bar * (1 - p_bar) / n_vec))
      lm  <- p_bar
      tit <- "Carta p (Fração Não-Conforme)"
      ylab <- "Fração Defeituosa (p)"
      y_plot <- taxa_obs
    } else {
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
    tit <- "Carta np (Número de Defeituosos)"
    ylab <- "Contagem de Não-Conformes (np)"
    y_plot <- dados
  } else {
    # tipo == "c"
    c_bar <- mean(dados)
    lm <- c_bar
    lsc <- c_bar + 3 * sqrt(c_bar)
    lic <- max(0, c_bar - 3 * sqrt(c_bar))
    tit <- "Carta c (Contagem de Defeitos)"
    ylab <- "Contagem de Defeitos (c)"
    y_plot <- dados
  }

  fora <- (y_plot > lsc) | (y_plot < lic)
  alarmes_idx <- which(fora)
  estavel <- length(alarmes_idx) == 0

  # ── Saída console ────────────────────────────────────────────────────────────
  .print_titulo(sprintf("CONTROLE ESTATÍSTICO DE PROCESSO: %s", toupper(tit)))
  cat(sprintf("  Subgrupos / Amostras Inspecionadas: %d\n\n", k))

  .print_topico("LIMITES DE CONTROLE ESTATÍSTICO (3σ)")
  w_lim <- c(24, 16, 16, 16)
  cat(.th_sep(sum(w_lim)), "\n")
  cat("  ", .th_pad("Parâmetro", w_lim[1], "left"), .th_pad("LIC (Inferior)", w_lim[2], "center"),
      .th_pad("Linha Média", w_lim[3], "center"), .th_pad("LSC (Superior)", w_lim[4], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_lim)), "\n")

  if (length(lsc) == 1) {
    cat("  ", .th_pad("Limites Constantes", w_lim[1], "left"),
        .th_pad(formatC(lic, format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
        .th_pad(formatC(lm, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
        .th_pad(formatC(lsc, format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")
  } else {
    cat("  ", .th_pad("Limites Variáveis (médios)", w_lim[1], "left"),
        .th_pad(formatC(mean(lic), format = "f", digits = 3, decimal.mark = ","), w_lim[2], "center"),
        .th_pad(formatC(lm, format = "f", digits = 3, decimal.mark = ","), w_lim[3], "center"),
        .th_pad(formatC(mean(lsc), format = "f", digits = 3, decimal.mark = ","), w_lim[4], "center"), "\n", sep = "")
  }
  cat(.th_sep(sum(w_lim)), "\n\n")

  .print_topico("INTERPRETAÇÃO E ESTABILIDADE DO PROCESSO")
  if (estavel) {
    cat("  Processo sob controle estatístico (estável): nenhuma amostra violou os limites 3σ.\n\n")
  } else {
    cat(sprintf("  Processo fora de controle estatístico (instável): %d amostra(s) além dos limites 3σ de controle:\n", length(alarmes_idx)))
    for (idx in alarmes_idx) {
      cat(sprintf("    * Amostra %d (Valor = %s): além dos limites 3σ\n",
                  idx, formatC(y_plot[idx], format = "f", digits = 3, decimal.mark = ",")))
    }
    cat("\n")
  }
  .print_rodape()

  p_out <- NULL
  if (grafico) {
    p_out <- .montar_grafico_shewhart(
      pontos_x = 1:k,
      valores  = y_plot,
      lc       = lm,
      lic      = lic,
      lsc      = lsc,
      titulo   = tit,
      xlab     = "Amostra / Subgrupo",
      ylab     = ylab,
      fora     = fora
    )
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


# ── Funções Diretas e Aliases Exportados ──────────────────────────────────────

#' @title Carta de Controle X (Valores Individuais)
#' @description Constrói a carta de controle para valores individuais isolada.
#' @inheritParams grafico_controle
#' @export
carta_x <- function(dados, regras = TRUE, grafico = TRUE) {
  grafico_controle(dados, tipo = "x", regras = regras, grafico = grafico)
}

#' @title Carta de Controle I (Valores Individuais)
#' @description Alias para \code{carta_x}.
#' @inheritParams grafico_controle
#' @export
carta_i <- carta_x

#' @title Carta de Controle Xbar (Médias)
#' @description Constrói a carta de controle de médias (\eqn{\bar{X}}) isolada.
#' @inheritParams grafico_controle
#' @export
carta_xbar <- function(dados, subgrupo = NULL, regras = TRUE, grafico = TRUE) {
  grafico_controle(dados, tipo = "xbar", subgrupo = subgrupo, regras = regras, grafico = grafico)
}

#' @title Carta de Controle R (Amplitude)
#' @description Constrói a carta de controle de amplitudes (R) isolada.
#' @inheritParams grafico_controle
#' @export
carta_r <- function(dados, subgrupo = NULL, regras = TRUE, grafico = TRUE) {
  grafico_controle(dados, tipo = "r", subgrupo = subgrupo, regras = regras, grafico = grafico)
}

#' @title Carta de Controle S (Desvio-Padrão)
#' @description Constrói a carta de controle de desvios-padrão (S) isolada.
#' @inheritParams grafico_controle
#' @export
carta_s <- function(dados, subgrupo = NULL, regras = TRUE, grafico = TRUE) {
  grafico_controle(dados, tipo = "s", subgrupo = subgrupo, regras = regras, grafico = grafico)
}

#' @title Carta de Controle c (Contagem de Defeitos)
#' @description Constrói a carta de controle c para contagem de defeitos por amostra de tamanho constante.
#' @inheritParams grafico_controle
#' @export
carta_c <- function(dados, regras = TRUE, grafico = TRUE) {
  grafico_controle(dados, tipo = "c", regras = regras, grafico = grafico)
}

#' @title Carta de Controle np (Número de Defeituosos)
#' @description Constrói a carta de controle np para número de itens não-conformes em amostras de tamanho constante.
#' @inheritParams grafico_controle
#' @export
carta_np <- function(dados, tamanho_amostra, regras = TRUE, grafico = TRUE) {
  grafico_controle(dados, tipo = "np", tamanho_amostra = tamanho_amostra, regras = regras, grafico = grafico)
}

#' @title Carta de Controle p (Fração Defeituosa)
#' @description Constrói a carta de controle p para fração de não-conformes.
#' @inheritParams grafico_controle
#' @export
carta_p <- function(dados, tamanho_amostra, regras = TRUE, grafico = TRUE) {
  grafico_controle(dados, tipo = "p", tamanho_amostra = tamanho_amostra, regras = regras, grafico = grafico)
}

#' @title Carta de Controle u (Taxa de Defeitos por Unidade)
#' @description Constrói a carta de controle u para taxa de defeitos por unidade de inspeção.
#' @inheritParams grafico_controle
#' @export
carta_u <- function(dados, tamanho_amostra, regras = TRUE, grafico = TRUE) {
  grafico_controle(dados, tipo = "u", tamanho_amostra = tamanho_amostra, regras = regras, grafico = grafico)
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
  .print_titulo("ANÁLISE DE CAPACIDADE DO PROCESSO (SIX SIGMA)")

  cat(sprintf("  Amostras (n): %d  |  Média: %.4f  |  Alvo (Target): %s\n",
              n, media_obs, if (!is.null(alvo)) formatC(alvo, format = "f", digits = 2) else "N/A"))
  cat(sprintf("  Especificações: LIE = %s  |  LSE = %s\n",
              if (!is.null(lie)) formatC(lie, format = "f", digits = 2) else "-Inf",
              if (!is.null(lse)) formatC(lse, format = "f", digits = 2) else "+Inf"))
  cat(sprintf("  Dispersão: σ Curto Prazo = %.4f  |  s Longo Prazo = %.4f\n\n", s_curto, s_longo))

  .print_topico("ÍNDICES DE CAPACIDADE E DESEMPENHO")
  w_cap <- c(28, 16)
  cat(.th_sep(sum(w_cap)), "\n")
  cat("  ", .th_pad("Índice", w_cap[1], "left"), .th_pad("Valor", w_cap[2], "center"), "\n", sep = "")
  cat(.th_sep(sum(w_cap)), "\n")

  for (i in 1:nrow(tb_indices)) {
    cat("  ", .th_pad(tb_indices$Indice[i], w_cap[1], "left"),
        .th_pad(formatC(tb_indices$Estimativa[i], format = "f", digits = 3, decimal.mark = ","), w_cap[2], "center"), "\n", sep = "")
  }
  cat(.th_sep(sum(w_cap)), "\n\n")

  .print_topico("DESEMPENHO DO PROCESSO (PPM E NÍVEL SIGMA)")
  cat(sprintf("  PPM Esperado (Teórico):   %.1f defeitos por milhão\n", ppm_esp_total))
  cat(sprintf("  PPM Observado na Amostra: %.1f defeitos por milhão (%d peças não-conformes)\n",
              ppm_obs_total, n_inf_obs + n_sup_obs))
  cat(sprintf("  Nível Sigma (σ-level):      %.2f σ\n\n", nivel_sigma))

  .print_topico("VEREDITO")
  cat(sprintf("  Diagnóstico: %s\n", classif))
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
        subtitle = sprintf("Cp = %.2f  |  Cpk = %.2f  |  PPM Esp = %.0f  |  Nível = %.2fσ",
                           if (!is.na(cp)) cp else cpk, cpk, ppm_esp_total, nivel_sigma),
        x = "Característica de Qualidade",
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

# ── 5. CARTA EWMA (MÉDIA MÓVEL EXPONENCIALMENTE PONDERADA) ───────────────────

#' @title Carta de Controle EWMA
#' @description Constrói carta de controle EWMA (Exponentially Weighted Moving Average)
#'   para monitoramento de processos e detecção precoce de pequenos deslocamentos de média
#'   (\eqn{\le 1{,}5\sigma}). Suporta limites dinâmicos no tempo ou assintóticos, com interface
#'   gráfica no padrão do ecossistema e destaque para observações fora de controle.
#'
#' @param dados Vetor numérico com as observações ou data frame quando utilizado via pipe (\code{|>}).
#' @param lambda Parâmetro de suavização/peso do histórico (\eqn{0 < \lambda \le 1}, padrão \code{0.20}).
#' @param L Multiplicador da largura dos limites em número de desvios-padrão (padrão \code{3.0}).
#' @param media Valor da média histórica de referência (\eqn{\mu_0}). Se \code{NULL}, utiliza a média amostral.
#' @param sigma Valor do desvio-padrão histórico do processo. Se \code{NULL}, estima via amplitude móvel média (\eqn{\overline{MR} / d_2}).
#' @param limites_variaveis Lógico. Se TRUE (padrão), utiliza limites dinâmicos no tempo que se expandem até a assíntota.
#' @param grafico Lógico. Se TRUE (padrão), exibe o gráfico no padrão visual estatR.
#'
#' @return Objeto invisível da classe \code{"estatR_carta_ewma"} contendo \code{z}, \code{lc}, \code{lic},
#'   \code{lsc}, \code{fora}, \code{n_fora} e o gráfico ggplot2.
#' @export
carta_ewma <- function(dados, lambda = 0.2, L = 3, media = NULL, sigma = NULL,
                       limites_variaveis = TRUE, grafico = TRUE) {

  if (is.data.frame(dados)) {
    num_cols <- which(sapply(dados, is.numeric))
    if (length(num_cols) == 0) stop("O data frame não contém colunas numéricas.", call. = FALSE)
    x <- as.numeric(dados[[num_cols[1]]])
  } else {
    x <- as.numeric(dados)
  }
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) stop("São necessárias pelo menos 4 observações para a carta EWMA.", call. = FALSE)

  mu0 <- if (!is.null(media)) as.numeric(media) else mean(x)

  if (!is.null(sigma)) {
    sig0 <- as.numeric(sigma)
  } else {
    mr <- abs(diff(x))
    sig0 <- mean(mr) / 1.128
  }

  z <- numeric(n)
  z_anterior <- mu0
  for (i in seq_len(n)) {
    z[i] <- lambda * x[i] + (1 - lambda) * z_anterior
    z_anterior <- z[i]
  }

  idx <- seq_len(n)
  if (limites_variaveis) {
    termo_tempo <- sqrt((lambda / (2 - lambda)) * (1 - (1 - lambda)^(2 * idx)))
    lic <- mu0 - L * sig0 * termo_tempo
    lsc <- mu0 + L * sig0 * termo_tempo
  } else {
    termo_assint <- sqrt(lambda / (2 - lambda))
    lic <- rep(mu0 - L * sig0 * termo_assint, n)
    lsc <- rep(mu0 + L * sig0 * termo_assint, n)
  }

  fora <- (z > lsc) | (z < lic)
  n_fora <- sum(fora)

  cat("\n── CARTA DE CONTROLE EWMA (MÉDIA MÓVEL EXPONENCIALMENTE PONDERADA) ────────\n\n")
  cat(sprintf("  Tamanho da Série      : %d observações\n", n))
  cat(sprintf("  Parâmetro de Peso (\u03bb) : %s\n", formatC(lambda, format = "f", digits = 2, decimal.mark = ",")))
  cat(sprintf("  Multiplicador Sigma (L): %s\n", formatC(L, format = "f", digits = 2, decimal.mark = ",")))
  cat(sprintf("  Média do Processo (\u03bc) : %s\n", formatC(mu0, format = "f", digits = 3, decimal.mark = ",")))
  cat(sprintf("  Desvio-Padrão (\u03c3)     : %s\n", formatC(sig0, format = "f", digits = 3, decimal.mark = ",")))
  cat(sprintf("  Limites Assintóticos  : LIC = %s  |  LSC = %s\n\n",
              formatC(mu0 - L * sig0 * sqrt(lambda / (2 - lambda)), format = "f", digits = 3, decimal.mark = ","),
              formatC(mu0 + L * sig0 * sqrt(lambda / (2 - lambda)), format = "f", digits = 3, decimal.mark = ",")))

  cat(sprintf("  Pontos Fora de Controle : %d de %d (%.1f%%)\n\n", n_fora, n, (n_fora / n) * 100))

  cat("  ▶ INTERPRETAÇÃO DO PROCESSO\n")
  if (n_fora == 0) {
    cat("  O processo opera sob controle estatístico em relação a pequenos desvios de média.\n")
    cat("  Nenhum ponto excedeu os limites dinâmicos da carta EWMA.\n")
  } else {
    posicoes <- which(fora)
    cat(sprintf("  Detectado(s) %d ponto(s) fora dos limites de controle (amostra(s): %s).\n",
                n_fora, paste(head(posicoes, 10), collapse = ", ")))
    cat("  A carta EWMA acusa uma tendência ou deslocamento sistemático na média do processo.\n")
  }
  cat("───────────────────────────────────────────────────────────────────────────\n\n")

  p <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    p <- .montar_grafico_shewhart(
      pontos_x = idx,
      valores  = z,
      lc       = rep(mu0, n),
      lic      = lic,
      lsc      = lsc,
      titulo   = "Carta de Controle EWMA",
      xlab     = "Ordem de Coleta / Subgrupo",
      ylab     = "Estatística EWMA (z)",
      fora     = fora
    )
    suppressMessages(print(p))
  }

  res <- list(
    z = z,
    lc = rep(mu0, n),
    lic = lic,
    lsc = lsc,
    fora = fora,
    n_fora = n_fora,
    lambda = lambda,
    L = L,
    grafico = p
  )
  class(res) <- c("estatR_carta_ewma", "list")
  invisible(res)
}


# ── 6. CARTA CUSUM (SOMAS CUMULATIVAS TABULARES) ─────────────────────────────

#' @title Carta de Controle CUSUM Tabular
#' @description Constrói carta de controle CUSUM (Cumulative Sum) bilateral tabular
#'   para monitoramento contínuo e detecção rápida de desvios na média do processo.
#'
#' @param dados Vetor numérico com as observações ou data frame quando utilizado via pipe (\code{|>}).
#' @param k Parâmetro de tolerância/referência em unidades de sigma (padrão \code{0.50}, ótimo para detectar desvio de \eqn{1\sigma}).
#' @param h Limite de decisão em unidades de sigma (padrão \code{4.50}).
#' @param media Valor da média histórica de referência (\eqn{\mu_0}). Se \code{NULL}, utiliza a média amostral.
#' @param sigma Valor do desvio-padrão histórico. Se \code{NULL}, estima via amplitude móvel média (\eqn{\overline{MR} / d_2}).
#' @param grafico Lógico. Se TRUE (padrão), plota as trajetórias de \eqn{C^+} e \eqn{C^-} com o limite \eqn{H}.
#'
#' @return Objeto invisível da classe \code{"estatR_carta_cusum"} contendo \code{c_pos}, \code{c_neg},
#'   \code{H}, \code{K}, \code{n_fora} e o gráfico ggplot2.
#' @export
carta_cusum <- function(dados, k = 0.5, h = 4.5, media = NULL, sigma = NULL, grafico = TRUE) {

  if (is.data.frame(dados)) {
    num_cols <- which(sapply(dados, is.numeric))
    if (length(num_cols) == 0) stop("O data frame não contém colunas numéricas.", call. = FALSE)
    x <- as.numeric(dados[[num_cols[1]]])
  } else {
    x <- as.numeric(dados)
  }
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) stop("São necessárias pelo menos 4 observações para a carta CUSUM.", call. = FALSE)

  mu0 <- if (!is.null(media)) as.numeric(media) else mean(x)

  if (!is.null(sigma)) {
    sig0 <- as.numeric(sigma)
  } else {
    mr <- abs(diff(x))
    sig0 <- mean(mr) / 1.128
  }

  K <- k * sig0
  H <- h * sig0

  c_pos <- numeric(n)
  c_neg <- numeric(n)
  cp_ant <- 0
  cn_ant <- 0

  for (i in seq_len(n)) {
    c_pos[i] <- max(0, x[i] - (mu0 + K) + cp_ant)
    c_neg[i] <- max(0, (mu0 - K) - x[i] + cn_ant)
    cp_ant <- c_pos[i]
    cn_ant <- c_neg[i]
  }

  fora_pos <- c_pos > H
  fora_neg <- c_neg > H
  fora_geral <- fora_pos | fora_neg
  n_fora <- sum(fora_geral)

  cat("\n── CARTA DE CONTROLE CUSUM (SOMAS CUMULATIVAS TABULARES) ──────────────────\n\n")
  cat(sprintf("  Tamanho da Série        : %d observações\n", n))
  cat(sprintf("  Parâmetro de Tolerância : k = %s  (K = %s)\n",
              formatC(k, format = "f", digits = 2, decimal.mark = ","),
              formatC(K, format = "f", digits = 3, decimal.mark = ",")))
  cat(sprintf("  Limite de Decisão       : h = %s  (H = %s)\n",
              formatC(h, format = "f", digits = 2, decimal.mark = ","),
              formatC(H, format = "f", digits = 3, decimal.mark = ",")))
  cat(sprintf("  Média do Processo (\u03bc)   : %s\n", formatC(mu0, format = "f", digits = 3, decimal.mark = ",")))
  cat(sprintf("  Desvio-Padrão (\u03c3)       : %s\n\n", formatC(sig0, format = "f", digits = 3, decimal.mark = ",")))

  cat(sprintf("  Pontos em Alarme (C+ ou C- > H): %d de %d (%.1f%%)\n\n",
              n_fora, n, (n_fora / n) * 100))

  cat("  ▶ INTERPRETAÇÃO DO PROCESSO\n")
  if (n_fora == 0) {
    cat("  O processo opera sob controle estatístico. As somas cumulativas permanecem dentro da faixa H.\n")
  } else {
    pos_p <- which(fora_pos)
    pos_n <- which(fora_neg)
    if (length(pos_p) > 0) {
      cat(sprintf("  Deslocamento positivo detectado na média: %d ponto(s) (amostra(s): %s).\n",
                  length(pos_p), paste(head(pos_p, 8), collapse = ", ")))
    }
    if (length(pos_n) > 0) {
      cat(sprintf("  Deslocamento negativo detectado na média: %d ponto(s) (amostra(s): %s).\n",
                  length(pos_n), paste(head(pos_n, 8), collapse = ", ")))
    }
  }
  cat("───────────────────────────────────────────────────────────────────────────\n\n")

  p <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    idx <- seq_len(n)
    df_p <- data.frame(Ponto = idx, CUSUM = c_pos, Tipo = "C+ (Superior)", Fora = fora_pos)
    df_n <- data.frame(Ponto = idx, CUSUM = c_neg, Tipo = "C- (Inferior)", Fora = fora_neg)
    df_plot <- rbind(df_p, df_n)

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Ponto, y = CUSUM, color = Tipo)) +
      ggplot2::geom_hline(yintercept = H, linetype = "dashed", color = "#D90429", linewidth = 0.9) +
      ggplot2::geom_hline(yintercept = 0, color = "#222222", linewidth = 0.8) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2.8) +
      ggplot2::scale_color_manual(values = c("C+ (Superior)" = "#005EA5", "C- (Inferior)" = "#E66101"))

    df_fora_c <- subset(df_plot, Fora)
    if (nrow(df_fora_c) > 0) {
      p <- p +
        ggplot2::geom_point(data = df_fora_c, color = "#D90429", shape = 1, size = 6.2, stroke = 1.3, show.legend = FALSE) +
        ggplot2::geom_point(data = df_fora_c, color = "#D90429", size = 3.2, show.legend = FALSE)
    }

    p <- p +
      ggplot2::labs(
        title = "Carta de Controle CUSUM Tabular",
        subtitle = sprintf("H = %s (Limite de Decisão)  |  K = %s  |  Alarme: %d ponto(s)",
                           formatC(H, format = "f", digits = 2, decimal.mark = ","),
                           formatC(K, format = "f", digits = 2, decimal.mark = ","),
                           n_fora),
        x = "Ordem de Coleta / Subgrupo",
        y = "Soma Cumulativa (CUSUM)"
      )

    if (exists("meu_tema", mode = "function")) {
      p <- p + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p <- p + tema_estatR(estilo = 2, grade = "dupla")
    } else {
      p <- p + ggplot2::theme_minimal()
    }
    p <- p + ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
    suppressMessages(print(p))
  }

  res <- list(
    c_pos = c_pos,
    c_neg = c_neg,
    H = H,
    K = K,
    fora_pos = fora_pos,
    fora_neg = fora_neg,
    n_fora = n_fora,
    grafico = p
  )
  class(res) <- c("estatR_carta_cusum", "list")
  invisible(res)
}


# ── 7. ESTUDO GAGE R&R (ANÁLISE DE SISTEMAS DE MEDIÇÃO — MSA) ────────────────

#' @title Análise do Sistema de Medição (Gage R&R)
#' @description Realiza estudo de Repetibilidade e Reprodutibilidade (Gage R&R) pelo método
#'   ANOVA de dois fatores com interação (modelo com efeitos aleatórios), calculando
#'   componentes de variância, percentuais de contribuição, variação do estudo (\%SV),
#'   número de categorias distintas (\eqn{ndc}) e veredito segundo critérios AIAG.
#'
#' @aliases gage_rr
#' @param dados Data frame contendo os dados do estudo ou fórmula do tipo \code{medicao ~ peca + operador}.
#' @param peca Nome da coluna de peças quando \code{dados} for data frame, ou data frame quando usado com fórmula.
#' @param operador Nome da coluna de operadores/avaliadores.
#' @param medicao Nome da coluna com as respostas numéricas medidas.
#' @param tolerancia Tolerância opcional do processo (LSE - LIE) para cálculo de \%Tolerância.
#' @param confianca Nível de confiança (padrão \code{0.95}).
#' @param grafico Lógico. Se TRUE (padrão), plota os componentes de variabilidade.
#'
#' @return Objeto invisível da classe \code{"estatR_gage_rr"} contendo a tabela ANOVA,
#'   componentes de variância, \code{ndc}, percentuais e o gráfico ggplot2.
#' @export
estudo_gage_rr <- function(dados, peca = NULL, operador = NULL, medicao = NULL,
                           tolerancia = NULL, confianca = 0.95, grafico = TRUE) {

  s_peca <- deparse(substitute(peca))
  s_op   <- deparse(substitute(operador))
  s_med  <- deparse(substitute(medicao))

  if (inherits(dados, "formula")) {
    formula_in <- dados
    dados_df <- peca
    mf <- stats::model.frame(formula_in, data = dados_df)
    med_vec <- stats::model.response(mf)
    peca_vec <- mf[[2]]
    op_vec <- mf[[3]]
    nome_peca <- names(mf)[2]
    nome_op   <- names(mf)[3]
    nome_med  <- names(mf)[1]
  } else {
    dados_df <- as.data.frame(dados)
    col_peca <- if (s_peca %in% names(dados_df)) s_peca else if (is.character(peca) && peca %in% names(dados_df)) peca else names(dados_df)[1]
    col_op   <- if (s_op %in% names(dados_df)) s_op else if (is.character(operador) && operador %in% names(dados_df)) operador else names(dados_df)[2]
    col_med  <- if (s_med %in% names(dados_df)) s_med else if (is.character(medicao) && medicao %in% names(dados_df)) medicao else names(dados_df)[3]

    peca_vec <- dados_df[[col_peca]]
    op_vec   <- dados_df[[col_op]]
    med_vec  <- as.numeric(dados_df[[col_med]])
    nome_peca <- col_peca
    nome_op   <- col_op
    nome_med  <- col_med
  }

  peca_fac <- as.factor(peca_vec)
  op_fac   <- as.factor(op_vec)
  p_num <- nlevels(peca_fac)
  o_num <- nlevels(op_fac)
  N_tot <- length(med_vec)
  r_num <- N_tot / (p_num * o_num)

  if (p_num < 2) stop("O estudo Gage R&R requer pelo menos 2 peças.", call. = FALSE)
  if (o_num < 2) stop("O estudo Gage R&R requer pelo menos 2 operadores.", call. = FALSE)

  fit <- stats::aov(med_vec ~ peca_fac * op_fac)
  sm  <- summary(fit)[[1]]
  rownames(sm) <- trimws(rownames(sm))

  ms_peca <- sm["peca_fac", "Mean Sq"]
  ms_op   <- sm["op_fac", "Mean Sq"]
  ms_int  <- sm["peca_fac:op_fac", "Mean Sq"]
  ms_res  <- sm["Residuals", "Mean Sq"]
  p_int   <- sm["peca_fac:op_fac", "Pr(>F)"]

  var_repetibilidade <- ms_res
  var_interacao <- if (!is.na(p_int) && p_int < 0.25) max(0, (ms_int - ms_res) / r_num) else 0
  var_operador  <- max(0, (ms_op - if (var_interacao > 0) ms_int else ms_res) / (p_num * r_num))
  var_reprodutibilidade <- var_operador + var_interacao
  var_grr <- var_repetibilidade + var_reprodutibilidade
  var_peca <- max(0, (ms_peca - if (var_interacao > 0) ms_int else ms_res) / (o_num * r_num))
  var_total <- var_grr + var_peca

  dp_repetibilidade <- sqrt(var_repetibilidade)
  dp_reprodutibilidade <- sqrt(var_reprodutibilidade)
  dp_grr <- sqrt(var_grr)
  dp_peca <- sqrt(var_peca)
  dp_total <- sqrt(var_total)

  pct_contrib_grr  <- (var_grr / var_total) * 100
  pct_contrib_rep  <- (var_repetibilidade / var_total) * 100
  pct_contrib_op   <- (var_reprodutibilidade / var_total) * 100
  pct_contrib_peca <- (var_peca / var_total) * 100

  pct_sv_grr  <- (dp_grr / dp_total) * 100
  pct_sv_rep  <- (dp_repetibilidade / dp_total) * 100
  pct_sv_op   <- (dp_reprodutibilidade / dp_total) * 100
  pct_sv_peca <- (dp_peca / dp_total) * 100

  ndc <- max(1, floor(1.41 * (dp_peca / dp_grr)))

  df_comp <- data.frame(
    Componente = c("Gage R&R Total", "  Repetibilidade (EV)", "  Reprodutibilidade (AV)", "Peça a Peça (PV)", "Variação Total (TV)"),
    Variancia = c(var_grr, var_repetibilidade, var_reprodutibilidade, var_peca, var_total),
    DP = c(dp_grr, dp_repetibilidade, dp_reprodutibilidade, dp_peca, dp_total),
    Pct_Contrib = c(pct_contrib_grr, pct_contrib_rep, pct_contrib_op, pct_contrib_peca, 100),
    Pct_Estudo = c(pct_sv_grr, pct_sv_rep, pct_sv_op, pct_sv_peca, 100),
    stringsAsFactors = FALSE
  )

  cat("\n── ANÁLISE DO SISTEMA DE MEDIÇÃO (GAGE R&R — MÉTODO ANOVA) ────────────────\n\n")
  cat(sprintf("  Medição Avaliada : %s\n", nome_med))
  cat(sprintf("  Peças            : %d peças (%s)\n", p_num, nome_peca))
  cat(sprintf("  Operadores       : %d operadores (%s)\n", o_num, nome_op))
  cat(sprintf("  Réplicas         : %.0f réplicas por peça/operador\n", r_num))
  cat(sprintf("  Total de Ensaios : %d medições\n\n", N_tot))

  cat("  COMPONENTES DE VARIÂNCIA E CAPACIDADE DO SISTEMA DE MEDIÇÃO:\n")
  sep_linha <- paste0("  ", strrep("\u2500", 75))
  cat(sep_linha, "\n")
  
  pad_txt <- function(s, w) {
    s <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    paste0(s, strrep(" ", sp))
  }
  
  wc_m <- c(comp = 27, var = 13, dp = 13, pc = 12, pe = 12)
  hdr_m <- paste0(pad_txt("Fonte de Variação", wc_m["comp"]),
                  pad_txt("Variância", wc_m["var"]),
                  pad_txt("Desvio-Padrão", wc_m["dp"]),
                  pad_txt("%Contribuição", wc_m["pc"]),
                  pad_txt("%Var.Estudo", wc_m["pe"]))
  cat("  ", hdr_m, "\n", sep = "")
  cat(sep_linha, "\n")

  for (i in seq_len(nrow(df_comp))) {
    row <- df_comp[i, ]
    cat("  ",
        pad_txt(row$Componente, wc_m["comp"]),
        pad_txt(formatC(row$Variancia, format = "f", digits = 4, decimal.mark = ","), wc_m["var"]),
        pad_txt(formatC(row$DP, format = "f", digits = 4, decimal.mark = ","), wc_m["dp"]),
        pad_txt(paste0(formatC(row$Pct_Contrib, format = "f", digits = 2, decimal.mark = ","), "%"), wc_m["pc"]),
        pad_txt(paste0(formatC(row$Pct_Estudo, format = "f", digits = 2, decimal.mark = ","), "%"), wc_m["pe"]),
        "\n", sep = "")
  }
  cat(sep_linha, "\n\n")

  cat(sprintf("  Número de Categorias Distintas (ndc): %d\n", ndc))
  cat("  (Critério AIAG: ndc \u2265 5 indica sistema com resolução suficiente)\n\n")

  cat("  ▶ VEREDITO DO SISTEMA DE MEDIÇÃO (CRITÉRIOS AIAG)\n")
  if (pct_sv_grr < 10) {
    cat(sprintf("  O sistema de medição é EXCELENTE / ACEITÁVEL (%%GRR = %.2f%% < 10%%).\n", pct_sv_grr))
    cat("  A variação do instrumento e dos avaliadores é desprezível frente à variação das peças.\n")
  } else if (pct_sv_grr <= 30) {
    cat(sprintf("  O sistema de medição é MARGINAL / CONDICIONAL (%%GRR = %.2f%% entre 10%% e 30%%).\n", pct_sv_grr))
    cat("  Pode ser aceito dependendo da criticidade da aplicação e dos custos de melhoria.\n")
  } else {
    cat(sprintf("  O sistema de medição é INACEITÁVEL (%%GRR = %.2f%% > 30%%).\n", pct_sv_grr))
    cat("  O sistema necessita de melhorias urgentes (recalibração, treinamento ou troca de dispositivo).\n")
  }
  cat("───────────────────────────────────────────────────────────────────────────\n\n")

  p <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot_m <- data.frame(
      Componente = rep(c("Gage R&R", "Repetibilidade", "Reprodutibilidade", "Peça a Peça"), 2),
      Percentual = c(pct_contrib_grr, pct_contrib_rep, pct_contrib_op, pct_contrib_peca,
                     pct_sv_grr, pct_sv_rep, pct_sv_op, pct_sv_peca),
      Metrica    = rep(c("% Contribuição da Variância", "% Variação do Estudo (%SV)"), each = 4)
    )

    p <- ggplot2::ggplot(df_plot_m, ggplot2::aes(x = Componente, y = Percentual, fill = Metrica)) +
      ggplot2::geom_col(position = "dodge", width = 0.65) +
      ggplot2::geom_hline(yintercept = c(10, 30), linetype = "dashed", color = c("#008000", "#D90429"), linewidth = 0.7) +
      ggplot2::scale_fill_manual(values = c("% Contribuição da Variância" = "#1E5B94", "% Variação do Estudo (%SV)" = "#E66101")) +
      ggplot2::labs(
        title = "Componentes de Variação do Sistema de Medição (Gage R&R)",
        subtitle = sprintf("%%GRR = %.2f%%  |  ndc = %d categorias  |  Linhas de referência AIAG em 10%% e 30%%", pct_sv_grr, ndc),
        x = NULL,
        y = "Percentual (%)"
      )

    if (exists("meu_tema", mode = "function")) {
      p <- p + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p <- p + tema_estatR(estilo = 2, grade = "dupla")
    } else {
      p <- p + ggplot2::theme_minimal()
    }
    p <- p + ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
    suppressMessages(print(p))
  }

  res <- list(
    anova = sm,
    componentes = df_comp,
    ndc = ndc,
    pct_grr = pct_sv_grr,
    pct_contrib = pct_contrib_grr,
    grafico = p
  )
  class(res) <- c("estatR_gage_rr", "list")
  invisible(res)
}

#' @rdname estudo_gage_rr
#' @export
gage_rr <- estudo_gage_rr
