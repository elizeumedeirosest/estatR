# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# M\u00d3DULO: TESTE DE EQUIVAL\u00caNCIA (TOST)
# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500

#' @title Teste de Equival\u00eancia (TOST)
#' @description Realiza o procedimento TOST (Two One-Sided Tests) para avaliar
#'   se dois grupos s\u00e3o estatisticamente equivalentes dentro de uma margem
#'   pr\u00e9-definida. Adequado para estudos de biossimilaridade, n\u00e3o-inferioridade
#'   e compara\u00e7\u00e3o de m\u00e9todos.
#' @param x Vetor num\u00e9rico do grupo 1, ou f\u00f3rmula no formato \code{y ~ grupo}.
#' @param y Vetor num\u00e9rico do grupo 2 (opcional se \code{x} for f\u00f3rmula).
#' @param dados Data frame (opcional se \code{x} for f\u00f3rmula).
#' @param margem Margem de equival\u00eancia \u0394 (obrigat\u00f3rio). Define o intervalo
#'   \code{[-\u0394, +\u0394]} dentro do qual a diferen\u00e7a deve estar para equivaler.
#' @param confianca N\u00edvel de confian\u00e7a (padr\u00e3o \code{0.95}).
#' @param grafico L\u00f3gico. Se \code{TRUE}, exibe gr\u00e1fico de equival\u00eancia.
#' @return Retorna invisivelmente uma lista com a diferen\u00e7a observada, IC do
#'   TOST, p-valores dos dois testes e o veredito de equival\u00eancia.
#' @export
teste_equivalencia <- function(x,
                               y         = NULL,
                               dados     = NULL,
                               margem,
                               confianca = 0.95,
                               grafico   = TRUE) {

  # ── Parsear entradas ─────────────────────────────────────────────────────────
  if (inherits(x, "formula")) {
    if (is.null(dados))
      stop("Forne\u00e7a 'dados' ao usar f\u00f3rmula.", call. = FALSE)
    resp_nome  <- as.character(x[[2]])
    grupo_nome <- as.character(x[[3]])
    if (!resp_nome  %in% names(dados))
      stop(sprintf("Coluna '%s' n\u00e3o encontrada.", resp_nome),  call. = FALSE)
    if (!grupo_nome %in% names(dados))
      stop(sprintf("Coluna '%s' n\u00e3o encontrada.", grupo_nome), call. = FALSE)

    grupos_fator <- as.factor(dados[[grupo_nome]])
    niveis       <- levels(grupos_fator)
    if (length(niveis) != 2)
      stop("A vari\u00e1vel de grupo deve ter exatamente 2 n\u00edveis.", call. = FALSE)

    x_vec    <- dados[[resp_nome]][grupos_fator == niveis[1]]
    y_vec    <- dados[[resp_nome]][grupos_fator == niveis[2]]
    nome_g1  <- niveis[1]
    nome_g2  <- niveis[2]
  } else {
    if (!is.numeric(x))
      stop("'x' deve ser num\u00e9rico ou uma f\u00f3rmula.", call. = FALSE)
    if (is.null(y) || !is.numeric(y))
      stop("Forne\u00e7a 'y' como vetor num\u00e9rico ao n\u00e3o usar f\u00f3rmula.", call. = FALSE)
    x_vec   <- x
    y_vec   <- y
    nome_g1 <- deparse(substitute(x))
    nome_g2 <- deparse(substitute(y))
    # Limpar nomes (remover $coluna se necess\u00e1rio)
    nome_g1 <- sub(".*\\$", "", nome_g1)
    nome_g2 <- sub(".*\\$", "", nome_g2)
  }

  x_vec <- x_vec[!is.na(x_vec)]
  y_vec <- y_vec[!is.na(y_vec)]

  if (length(x_vec) < 2 || length(y_vec) < 2)
    stop("Cada grupo deve ter pelo menos 2 observa\u00e7\u00f5es.", call. = FALSE)
  if (missing(margem) || !is.numeric(margem) || margem <= 0)
    stop("Forne\u00e7a uma margem de equival\u00eancia positiva em 'margem'.", call. = FALSE)

  alpha <- 1 - confianca

  # ── Estat\u00edsticas descritivas ───────────────────────────────────────────────
  n1     <- length(x_vec)
  n2     <- length(y_vec)
  m1     <- mean(x_vec)
  m2     <- mean(y_vec)
  dp1    <- sd(x_vec)
  dp2    <- sd(y_vec)
  ep1    <- dp1 / sqrt(n1)
  ep2    <- dp2 / sqrt(n2)

  # ── C\u00e1lculo TOST ─────────────────────────────────────────────────────────────
  diff_obs <- m1 - m2
  EP_diff  <- sqrt(dp1^2 / n1 + dp2^2 / n2)

  # Graus de liberdade de Welch (Satterthwaite)
  df_welch <- (dp1^2 / n1 + dp2^2 / n2)^2 /
              ((dp1^2 / n1)^2 / (n1 - 1) + (dp2^2 / n2)^2 / (n2 - 1))

  # Testes unicaudais TOST
  t1   <- (diff_obs - (-margem)) / EP_diff   # H0: diff <= -Delta
  t2   <- (diff_obs - (+margem)) / EP_diff   # H0: diff >= +Delta

  p1   <- pt(t1, df = df_welch, lower.tail = FALSE)   # P(T > t1)
  p2   <- pt(t2, df = df_welch, lower.tail = TRUE)    # P(T < t2)

  p_tost <- max(p1, p2)

  # IC com n\u00edvel (1 - 2*alpha) = confianca - alpha = ex: 90% para alpha=0.05
  t_ic   <- qt(1 - alpha, df = df_welch)           # quantil para IC (1-2*alpha)
  ic_inf <- diff_obs - t_ic * EP_diff
  ic_sup <- diff_obs + t_ic * EP_diff

  ic_pct <- (1 - 2 * alpha) * 100

  # Veredito
  equivalente <- (ic_inf > -margem) && (ic_sup < margem)

  # ── Imprimir sa\u00edda ──────────────────────────────────────────────────────────
  w_sep <- 70

  .print_titulo("TESTE DE EQUIVAL\u00caNCIA (TOST)")

  cat(sprintf(
    "  Grupos: %s vs %s\n",
    nome_g1, nome_g2
  ))
  cat(sprintf(
    "  Margem de equival\u00eancia: \u0394 = \u00b1%.2f\n\n",
    margem
  ))

  cat("  H0: A diferen\u00e7a EST\u00c1 FORA da margem de equival\u00eancia (|diff| \u2265 \u0394)\n")
  cat("  H1: A diferen\u00e7a EST\u00c1 DENTRO da margem (|diff| < \u0394) \u2014 s\u00e3o equivalentes\n\n")

  # ── Bloco: Descritivas ───────────────────────────────────────────────────────
  .print_topico("ESTAT\u00cdSTICAS DESCRITIVAS")

  wc_d <- c(16, 8, 10, 10, 10)
  hdr_d <- paste0(
    "  ",
    .th_pad("Grupo",   wc_d[1], "left"),
    .th_pad("N",       wc_d[2], "center"),
    .th_pad("M\u00e9dia",  wc_d[3], "center"),
    .th_pad("DP",      wc_d[4], "center"),
    .th_pad("EP",      wc_d[5], "center")
  )
  sep_d <- paste0("  ", strrep("\u2500", nchar(hdr_d) - 2))
  cat(sep_d, "\n")
  cat(hdr_d, "\n")
  cat(sep_d, "\n")

  for (i in 1:2) {
    gn  <- c(nome_g1, nome_g2)[i]
    ni  <- c(n1, n2)[i]
    mi  <- c(m1, m2)[i]
    dpi <- c(dp1, dp2)[i]
    epi <- c(ep1, ep2)[i]
    cat(
      "  ",
      .th_pad(gn,                                                              wc_d[1], "left"),
      .th_pad(ni,                                                              wc_d[2], "center"),
      .th_pad(formatC(mi,  format = "f", digits = 2, decimal.mark = ","),     wc_d[3], "center"),
      .th_pad(formatC(dpi, format = "f", digits = 2, decimal.mark = ","),     wc_d[4], "center"),
      .th_pad(formatC(epi, format = "f", digits = 2, decimal.mark = ","),     wc_d[5], "center"),
      "\n", sep = ""
    )
  }
  cat(sep_d, "\n")
  cat(sprintf(
    "  Diferen\u00e7a observada (%s \u2212 %s):  %.2f\n\n",
    nome_g1, nome_g2,
    diff_obs
  ))

  # ── Bloco: Teste TOST ────────────────────────────────────────────────────────
  .print_topico("TESTE TOST")

  wc_t <- c(38, 22)
  hdr_t <- paste0(
    "  ",
    .th_pad("M\u00e9trica",  wc_t[1], "left"),
    .th_pad("Valor",   wc_t[2], "center")
  )
  sep_t <- paste0("  ", strrep("\u2500", nchar(hdr_t) - 2))
  cat(sep_t, "\n")
  cat(hdr_t, "\n")
  cat(sep_t, "\n")

  cat(
    "  ",
    .th_pad(sprintf("Teste 1 (H0: diff \u2264 -\u0394)    t = %s",
                    formatC(t1, format = "f", digits = 2)), wc_t[1], "left"),
    .th_pad(sprintf("p = %s", formatC(p1, format = "f", digits = 3, decimal.mark = ",")),
            wc_t[2], "center"),
    "\n", sep = ""
  )
  cat(
    "  ",
    .th_pad(sprintf("Teste 2 (H0: diff \u2265 +\u0394)    t = %s",
                    formatC(t2, format = "f", digits = 2)), wc_t[1], "left"),
    .th_pad(sprintf("p = %s", formatC(p2, format = "f", digits = 3, decimal.mark = ",")),
            wc_t[2], "center"),
    "\n", sep = ""
  )
  cat(
    "  ",
    .th_pad("p-valor TOST (max)",  wc_t[1], "left"),
    .th_pad(paste(
      formatC(p_tost, format = "f", digits = 3, decimal.mark = ","),
      .th_asterisk(p_tost)
    ), wc_t[2], "center"),
    "\n", sep = ""
  )
  cat(
    "  ",
    .th_pad(sprintf("IC %.0f%% da diferen\u00e7a", ic_pct), wc_t[1], "left"),
    .th_pad(sprintf("[%s; %s]",
                    formatC(ic_inf, format = "f", digits = 2, decimal.mark = ","),
                    formatC(ic_sup, format = "f", digits = 2, decimal.mark = ",")),
            wc_t[2], "center"),
    "\n", sep = ""
  )
  cat(sep_t, "\n\n")

  # ── Veredito ─────────────────────────────────────────────────────────────────
  if (equivalente) {
    cat(paste0("  ", "\033[1m", "\u25b6 VEREDITO: \u2713 EQUIVALENTES", "\033[0m", "\n"))
    cat(sprintf(
      "  Os grupos s\u00e3o estatisticamente equivalentes ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n",
      alpha * 100
    ))
    cat(sprintf(
      "  O IC %.0f%% da diferen\u00e7a [%s; %s] est\u00e1 inteiramente contido na margem [-%.2f; +%.2f].\n\n",
      ic_pct,
      formatC(ic_inf, format = "f", digits = 2, decimal.mark = ","),
      formatC(ic_sup, format = "f", digits = 2, decimal.mark = ","),
      margem, margem
    ))
  } else {
    cat(paste0("  ", "\033[1m", "\u25b6 VEREDITO: \u2717 N\u00c3O EQUIVALENTES", "\033[0m", "\n"))
    cat(sprintf(
      "  Os grupos n\u00e3o s\u00e3o estatisticamente equivalentes ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n",
      alpha * 100
    ))
    cat(sprintf(
      "  O IC %.0f%% da diferen\u00e7a [%s; %s] n\u00e3o est\u00e1 inteiramente contido na margem [-%.2f; +%.2f].\n\n",
      ic_pct,
      formatC(ic_inf, format = "f", digits = 2, decimal.mark = ","),
      formatC(ic_sup, format = "f", digits = 2, decimal.mark = ","),
      margem, margem
    ))
  }

  .print_rodape()

  # ── Gr\u00e1fico de equival\u00eancia ────────────────────────────────────────────────
  if (grafico) {
    tryCatch({
      df_graf <- data.frame(
        x    = 0,
        diff = diff_obs,
        inf  = ic_inf,
        sup  = ic_sup
      )

      cor_ponto <- if (equivalente) "#377EB8" else "#E31A1C"
      label_verd <- if (equivalente) "Equivalentes" else "N\u00e3o equivalentes"

      tema_fn <- if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }

      p_graf <- ggplot2::ggplot(df_graf, ggplot2::aes(x = diff, y = 0)) +
        # Zona de equival\u00eancia
        ggplot2::annotate(
          "rect",
          xmin = -margem, xmax = margem,
          ymin = -0.5, ymax = 0.5,
          fill = "gray80", alpha = 0.5
        ) +
        # Linhas de margem
        ggplot2::geom_vline(xintercept = c(-margem, margem),
                            colour = "gray50", linetype = "dashed", linewidth = 0.7) +
        # Linha de refer\u00eancia em zero
        ggplot2::geom_vline(xintercept = 0,
                            colour = "firebrick", linetype = "dashed", linewidth = 0.7) +
        # Barra de IC
        ggplot2::geom_segment(
          ggplot2::aes(x = inf, xend = sup, y = 0, yend = 0),
          colour = cor_ponto, linewidth = 1.5
        ) +
        # Tampas do IC
        ggplot2::geom_point(
          data = data.frame(x = c(ic_inf, ic_sup), y = 0),
          ggplot2::aes(x = x, y = y),
          shape = 3, size = 4, colour = cor_ponto
        ) +
        # Ponto central (diferen\u00e7a)
        ggplot2::geom_point(
          ggplot2::aes(x = diff, y = 0),
          size = 4, colour = cor_ponto
        ) +
        # Anotac\u00f5es de margem
        ggplot2::annotate("text", x = -margem, y = 0.45,
                          label = sprintf("-\u0394 = -%.2f", margem),
                          hjust = 1.1, size = 3.5, colour = "gray40") +
        ggplot2::annotate("text", x =  margem, y = 0.45,
                          label = sprintf("+\u0394 = +%.2f", margem),
                          hjust = -0.1, size = 3.5, colour = "gray40") +
        ggplot2::labs(
          title    = "Gr\u00e1fico de Equival\u00eancia (TOST)",
          subtitle = sprintf(
            "IC %.0f%% da diferen\u00e7a: [%.2f; %.2f] | Veredito: %s",
            ic_pct, ic_inf, ic_sup, label_verd
          ),
          x = sprintf("Diferen\u00e7a (%s \u2212 %s)", nome_g1, nome_g2),
          y = NULL
        ) +
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank()
        ) +
        tema_fn

      print(p_graf)
    }, error = function(e) {
      message("  [aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico: ", e$message)
    })
  }

  # ── Retorno invis\u00edvel ────────────────────────────────────────────────────────
  invisible(list(
    n1          = n1,
    n2          = n2,
    media1      = m1,
    media2      = m2,
    diff        = diff_obs,
    EP_diff     = EP_diff,
    df_welch    = df_welch,
    t1          = t1,
    t2          = t2,
    p1          = p1,
    p2          = p2,
    p_tost      = p_tost,
    ic_inf      = ic_inf,
    ic_sup      = ic_sup,
    ic_pct      = ic_pct,
    margem      = margem,
    equivalente = equivalente
  ))
}
