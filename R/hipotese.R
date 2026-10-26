# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: TESTES DE HIPÓTESE
# ─────────────────────────────────────────────────────────────────────────────

# ── Helpers internos ──────────────────────────────────────────────────────────

.th_fmt_p_interp <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0,001")
  paste0("p = ", formatC(p, format = "f", digits = 3, decimal.mark = ","))
}
.th_fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  formatC(p, format = "f", digits = 3, decimal.mark = ",")
}

.th_asterisk <- function(p) {
  if (is.na(p) || p >= 0.10) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  return("*")
}

.th_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
  paste0(s, strrep(" ", sp))
}

.th_sep <- function(w) paste0("  ", strrep("\u2500", w))

.th_nota_sig <- function(p) {
  if (!is.na(p) && p < 0.10)
    cat("  Signific\u00e2ncia: *** p < 0.01   ** p < 0.05   * p < 0.10\n")
}

.th_hipotese_labels <- function(hipotese, param = "\u03bc", ref = NULL, g1 = NULL, g2 = NULL) {
  op_h1 <- switch(hipotese,
    bilateral = "\u2260",
    maior     = ">",
    menor     = "<"
  )
  if (!is.null(ref)) {
    h0_str <- sprintf("H0: %s = %s", param, ref)
    h1_str <- sprintf("H1: %s %s %s (%s)", param, op_h1, ref, hipotese)
  } else {
    h0_str <- sprintf("H0: %s(%s) = %s(%s)", param, g1, param, g2)
    h1_str <- sprintf("H1: %s(%s) %s %s(%s) (%s)", param, g1, op_h1, param, g2, hipotese)
  }
  list(h0 = h0_str, h1 = h1_str)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teste T para Uma Amostra
#' @description Testa se a média de uma variável numérica é igual a um valor de referência.
#' @param x Vetor numérico.
#' @param mu Valor de referência para a média populacional (μ₀).
#' @param hipotese Direção do teste: \code{"bilateral"} (padrão), \code{"maior"} ou \code{"menor"}.
#' @param confianca Nível de confiança para o IC (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @param grafico Lógico. Se TRUE, exibe o gráfico da distribuição t com a região crítica.
#' @return Retorna invisivelmente o objeto do teste (\code{htest}).
#' @export
teste_t_uma_amostra <- function(x, mu, hipotese = c("bilateral", "maior", "menor"),
                                 confianca = 0.95, decimais = 2, grafico = TRUE) {
  hipotese <- match.arg(hipotese)
  if (!is.numeric(x)) stop("'x' deve ser um vetor num\u00e9rico.")

  alt_r <- switch(hipotese, bilateral = "two.sided", maior = "greater", menor = "less")

  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)
  if (identical(var_nome, "x")) var_nome <- "x"

  x_c    <- x[!is.na(x)]
  n      <- length(x_c)
  media  <- mean(x_c)
  dp     <- sd(x_c)
  ep     <- dp / sqrt(n)

  res    <- t.test(x_c, mu = mu, alternative = alt_r, conf.level = confianca)
  t_stat <- res$statistic
  gl     <- res$parameter
  p_val  <- res$p.value
  ic     <- res$conf.int
  p_ast  <- .th_asterisk(p_val)

  lbs    <- .th_hipotese_labels(hipotese, param = "\u03bc", ref = mu)
  w_sep  <- 70

  .print_titulo("TESTE T PARA UMA AMOSTRA")
  cat(sprintf("  Vari\u00e1vel: %s   |   μ0: %s\n\n", var_nome, mu))
  cat(sprintf("  %s\n", lbs$h0))
  cat(sprintf("  %s\n\n", lbs$h1))

  # Tabela de resultados
  wc <- c(med = 18, t = 15, gl = 12, ic = 20, p = 16)
  hdr <- paste0(
    .th_pad("M\u00e9dia amostral", wc["med"], "center"),
    .th_pad("Estat\u00edstica t",  wc["t"],   "center"),
    .th_pad("Graus Lib.",    wc["gl"],  "center"),
    .th_pad(sprintf("IC (%.0f%%)", confianca * 100), wc["ic"], "center"),
    .th_pad("p-valor",       wc["p"],   "center")
  )
  cat(.th_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.th_sep(nchar(hdr)), "\n")

  ic_str   <- sprintf("[%s; %s]",
                      formatC(ic[1], format = "f", digits = decimais, decimal.mark = ","),
                      formatC(ic[2], format = "f", digits = decimais, decimal.mark = ","))
  p_str    <- paste(.th_fmt_p(p_val), p_ast)
  media_str <- formatC(media, format = "f", digits = decimais, decimal.mark = ",")
  t_str    <- formatC(t_stat, format = "f", digits = 3, decimal.mark = ",")
  gl_str   <- formatC(gl,     format = "f", digits = 0, decimal.mark = ",")

  cat("  ",
      .th_pad(media_str, wc["med"], "center"),
      .th_pad(t_str,     wc["t"],   "center"),
      .th_pad(gl_str,    wc["gl"],  "center"),
      .th_pad(ic_str,    wc["ic"],  "center"),
      .th_pad(p_str,     wc["p"],   "center"),
      "\n", sep = "")
  cat(.th_sep(nchar(hdr)), "\n")
  .th_nota_sig(p_val)
  cat("\n")

  # Interpretação
  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_val < (1 - confianca)) {
    cat(sprintf("  H\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  A m\u00e9dia de %s (%s) \u00e9 estatisticamente diferente\n",
                var_nome, media_str))
    cat(sprintf("  de %s ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n",
                mu, (1 - confianca) * 100))
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  A m\u00e9dia de %s (%s) n\u00e3o difere significativamente\n",
                var_nome, media_str))
    cat(sprintf("  de %s ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n",
                mu, (1 - confianca) * 100))
  }
  .print_rodape()

  # Gráfico
  if (grafico) {
    tryCatch({
      if (exists("meu_tema")) {
        alpha_nivel <- 1 - confianca
        tc <- switch(hipotese,
          bilateral = qt(1 - alpha_nivel / 2, df = gl),
          maior     = qt(1 - alpha_nivel, df = gl),
          menor     = qt(alpha_nivel, df = gl)
        )

        x_seq  <- seq(-4, 4, length.out = 500)
        df_curv <- data.frame(x = x_seq, y = dt(x_seq, df = gl))

        p_plot <- ggplot2::ggplot(df_curv, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_line(linewidth = 0.9, color = "#333333")

        if (hipotese == "bilateral") {
          df_esq <- df_curv[df_curv$x <= -abs(tc), ]
          df_dir <- df_curv[df_curv$x >= abs(tc), ]
          p_plot <- p_plot +
            ggplot2::geom_area(data = df_esq, fill = "red1", alpha = 0.35) +
            ggplot2::geom_area(data = df_dir, fill = "red1", alpha = 0.35)
        } else if (hipotese == "maior") {
          df_crit <- df_curv[df_curv$x >= tc, ]
          p_plot <- p_plot + ggplot2::geom_area(data = df_crit, fill = "red1", alpha = 0.35)
        } else {
          df_crit <- df_curv[df_curv$x <= tc, ]
          p_plot <- p_plot + ggplot2::geom_area(data = df_crit, fill = "red1", alpha = 0.35)
        }

        p_plot <- p_plot +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
          ggplot2::geom_vline(xintercept = as.numeric(t_stat),
                              color = "#1B4F72", linetype = "dashed", linewidth = 1) +
          ggplot2::annotate("text", x = as.numeric(t_stat), y = max(df_curv$y) * 1.05,
                            label = sprintf("t = %s", formatC(as.numeric(t_stat), format="f", digits=3, decimal.mark=",")),
                            hjust = -0.1, color = "#1B4F72", fontface = "bold", size = 4.5) +
          ggplot2::labs(
            title    = sprintf("Distribuição t — Teste para média de %s", var_nome),
            subtitle = sprintf("%s | p %s %s%s", lbs$h1,
                               ifelse(p_val < 0.001, "<", "="), .th_fmt_p(p_val),
                               ifelse(nchar(p_ast) > 0, paste0(" ", p_ast), "")),
            x = "t",
            y = "Densidade"
          ) +
          meu_tema(grade = "dupla")

        suppressMessages(print(p_plot))
      }
    }, error = function(e) {
      message("[Aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico: ", e$message)
    })
  }

  invisible(res)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teste T para Duas Amostras Independentes
#' @description Compara as médias de dois grupos independentes.
#' @param x Vetor numérico com os valores, ou data.frame.
#' @param grupo Vetor com os grupos (deve ter exatamente 2 níveis), ou nome da coluna se \code{x} for data.frame.
#' @param hipotese Direção do teste: \code{"bilateral"} (padrão), \code{"maior"} ou \code{"menor"}.
#' @param variancia_igual Lógico. Se FALSE (padrão), usa o teste de Welch (recomendado).
#' @param confianca Nível de confiança para o IC (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @param grafico Lógico. Se TRUE, exibe boxplot comparativo dos dois grupos.
#' @return Retorna invisivelmente o objeto do teste (\code{htest}).
#' @export
teste_t_duas_amostras <- function(x, grupo, hipotese = c("bilateral", "maior", "menor"),
                                   variancia_igual = FALSE, confianca = 0.95,
                                   decimais = 2, grafico = TRUE) {
  hipotese <- match.arg(hipotese)

  # NSE: aceita data.frame + nome de coluna sem aspas
  grupo_sub  <- substitute(grupo)
  grupo_nome <- sub(".*\\$", "", deparse(grupo_sub))
  x_sub      <- substitute(x)
  x_nome     <- sub(".*\\$", "", deparse(x_sub))

  if (is.data.frame(x)) {
    grupo_vec <- if (grupo_nome %in% names(x)) x[[grupo_nome]] else eval(grupo_sub, parent.frame())
    x_vec     <- NULL  # será extraído por grupo abaixo
  } else {
    grupo_vec <- eval(grupo_sub, parent.frame())
    x_vec     <- x
  }

  grupo_vec <- as.factor(grupo_vec)
  niveis    <- levels(grupo_vec)
  if (length(niveis) != 2) stop("'grupo' deve ter exatamente 2 n\u00edveis.")

  g1 <- niveis[1]; g2 <- niveis[2]

  if (is.data.frame(x)) {
    vals1 <- x[[x_nome]][grupo_vec == g1 & !is.na(grupo_vec)]
    vals2 <- x[[x_nome]][grupo_vec == g2 & !is.na(grupo_vec)]
    var_nome <- x_nome
  } else {
    vals1 <- x_vec[grupo_vec == g1 & !is.na(grupo_vec)]
    vals2 <- x_vec[grupo_vec == g2 & !is.na(grupo_vec)]
    var_nome <- x_nome
  }

  vals1 <- vals1[!is.na(vals1)]
  vals2 <- vals2[!is.na(vals2)]

  alt_r <- switch(hipotese, bilateral = "two.sided", maior = "greater", menor = "less")
  res   <- t.test(vals1, vals2, alternative = alt_r, var.equal = variancia_igual,
                  conf.level = confianca)

  t_stat <- res$statistic
  gl     <- res$parameter
  p_val  <- res$p.value
  ic     <- res$conf.int
  p_ast  <- .th_asterisk(p_val)

  tipo_txt <- if (variancia_igual) "Student (vari\u00e2ncias iguais)" else "Welch (vari\u00e2ncias n\u00e3o assumidas iguais)"
  lbs <- .th_hipotese_labels(hipotese, param = "\u03bc", g1 = g1, g2 = g2)
  w_sep <- 70

  # ── Cabeçalho ──────────────────────────────────────────────────────────────
  .print_titulo("TESTE T PARA DUAS AMOSTRAS INDEPENDENTES")
  cat(sprintf("  Vari\u00e1vel: %s   |   Grupos: %s (%s vs %s)\n", var_nome, grupo_nome, g1, g2))
  cat(sprintf("  Tipo: %s\n\n", tipo_txt))
  cat(sprintf("  %s\n", lbs$h0))
  cat(sprintf("  %s\n\n", lbs$h1))

  # ── Tabela descritiva por grupo ────────────────────────────────────────────
  wg <- c(grp = 16, n = 8, med = 14, dp = 12)
  hdr_g <- paste0(
    .th_pad("Grupo",  wg["grp"], "left"),
    .th_pad("N",      wg["n"],   "center"),
    .th_pad("M\u00e9dia", wg["med"], "center"),
    .th_pad("DP",     wg["dp"],  "center")
  )
  cat(.th_sep(nchar(hdr_g)), "\n")
  cat("  ", hdr_g, "\n", sep = "")
  cat(.th_sep(nchar(hdr_g)), "\n")

  fmt_n <- function(vals, rot) {
    cat("  ",
        .th_pad(rot, wg["grp"], "left"),
        .th_pad(length(vals), wg["n"], "center"),
        .th_pad(formatC(mean(vals), format="f", digits=decimais, decimal.mark=","), wg["med"], "center"),
        .th_pad(formatC(sd(vals),   format="f", digits=decimais, decimal.mark=","), wg["dp"],  "center"),
        "\n", sep = "")
  }
  fmt_n(vals1, g1)
  fmt_n(vals2, g2)
  cat(.th_sep(nchar(hdr_g)), "\n\n")

  # ── Tabela de resultados ───────────────────────────────────────────────────
  wc <- c(t = 13, gl = 12, dif = 17, ic = 22, p = 16)
  ic_str <- sprintf("[%s; %s]",
                    formatC(ic[1], format="f", digits=decimais, decimal.mark=","),
                    formatC(ic[2], format="f", digits=decimais, decimal.mark=","))
  hdr_c <- paste0(
    .th_pad("Estat. t",        wc["t"],   "center"),
    .th_pad("Graus Lib.",      wc["gl"],  "center"),
    .th_pad("Dif. de M\u00e9dias", wc["dif"], "center"),
    .th_pad(sprintf("IC (%.0f%%)", confianca*100), wc["ic"], "center"),
    .th_pad("p-valor",         wc["p"],   "center")
  )
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ",
      .th_pad(formatC(as.numeric(t_stat), format="f", digits=3, decimal.mark=","), wc["t"],   "center"),
      .th_pad(formatC(as.numeric(gl),     format="f", digits=1, decimal.mark=","), wc["gl"],  "center"),
      .th_pad(formatC(mean(vals1)-mean(vals2), format="f", digits=decimais, decimal.mark=","), wc["dif"], "center"),
      .th_pad(ic_str, wc["ic"], "center"),
      .th_pad(paste(.th_fmt_p(p_val), p_ast), wc["p"], "center"),
      "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  .th_nota_sig(p_val)
  cat("\n")

  # ── Interpretação ──────────────────────────────────────────────────────────
  alpha <- 1 - confianca
  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_val < alpha) {
    cat(sprintf("  H\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  A m\u00e9dia de %s difere significativamente entre os grupos\n", var_nome))
    cat(sprintf("  %s (m\u00e9dia = %s) e %s (m\u00e9dia = %s),\n",
                g1, formatC(mean(vals1), format="f", digits=decimais, decimal.mark=","),
                g2, formatC(mean(vals2), format="f", digits=decimais, decimal.mark=",")))
    cat(sprintf("  ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n", alpha * 100))
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  A m\u00e9dia de %s n\u00e3o difere significativamente entre os grupos\n", var_nome))
    cat(sprintf("  %s (m\u00e9dia = %s) e %s (m\u00e9dia = %s),\n",
                g1, formatC(mean(vals1), format="f", digits=decimais, decimal.mark=","),
                g2, formatC(mean(vals2), format="f", digits=decimais, decimal.mark=",")))
    cat(sprintf("  ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n", alpha * 100))
  }
  .print_rodape()

  # ── Gráfico ────────────────────────────────────────────────────────────────
  if (grafico) {
    tryCatch({
      if (exists("meu_tema")) {
        p_sub <- sprintf("%s | p %s %s%s", lbs$h1,
                         ifelse(p_val < 0.001, "<", "="),
                         .th_fmt_p(p_val),
                         ifelse(nchar(p_ast) > 0, paste0(" ", p_ast), ""))

        df_plot <- data.frame(grupo = grupo_vec, y = x_vec)
        df_plot <- df_plot[!is.na(df_plot$grupo) & !is.na(df_plot$y), ]
        max_y <- max(df_plot$y, na.rm = TRUE)
        min_y <- min(df_plot$y, na.rm = TRUE)
        amp <- max_y - min_y
        y_bar <- max_y + amp * 0.08
        y_tick <- max_y + amp * 0.05
        y_text <- max_y + amp * 0.12

        old_w <- getOption("warn"); options(warn = -1)
        med1 <- mean(vals1); med2 <- mean(vals2)
        cor1 <- if (med1 >= med2) "#0072CE" else "#DDDDDD"
        cor2 <- if (med2 > med1) "#0072CE" else "#DDDDDD"
        cores_box <- setNames(c(cor1, cor2), c(g1, g2))

        p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = grupo, y = y, fill = grupo)) +
          ggplot2::stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.6) +
          ggplot2::geom_boxplot(width = 0.5, alpha = 0.85, color = "black", outlier.color = "red3", outlier.shape = 19) +
          ggplot2::stat_summary(fun = mean, geom = "point", shape = 21, size = 3.2, fill = "black", color = "black") +
          ggplot2::annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, color = "black", linewidth = 0.5) +
          ggplot2::annotate("segment", x = 1, xend = 1, y = y_tick, yend = y_bar, color = "black", linewidth = 0.5) +
          ggplot2::annotate("segment", x = 2, xend = 2, y = y_tick, yend = y_bar, color = "black", linewidth = 0.5) +
          ggplot2::annotate("text", x = 1.5, y = y_text, label = p_sub, size = 4.2, fontface = "bold") +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
          ggplot2::scale_fill_manual(values = cores_box) +
          ggplot2::labs(
            title    = sprintf("Comparação: %s por %s", var_nome, grupo_nome),
            x        = grupo_nome,
            y        = var_nome
          ) +
          meu_tema(estilo = 2, grade = "dupla") +
          ggplot2::theme(legend.position = "none")

        suppressMessages(print(p_plot))
        options(warn = old_w)
      }
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o gráfico: ", e$message)
    })
  }

  invisible(res)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teste T Pareado
#' @description Compara as médias de duas medições pareadas (antes/depois, par de medidas).
#' @param x Vetor numérico com as medidas do primeiro momento (ex: antes).
#' @param y Vetor numérico com as medidas do segundo momento (ex: depois).
#' @param hipotese Direção do teste: \code{"bilateral"} (padrão), \code{"maior"} ou \code{"menor"}.
#' @param confianca Nível de confiança para o IC (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @param grafico Lógico. Se TRUE, exibe gráfico de linhas conectando os pares antes/depois.
#' @return Retorna invisivelmente o objeto do teste (\code{htest}).
#' @export
teste_t_pareado <- function(x, y, hipotese = c("bilateral", "maior", "menor"),
                             confianca = 0.95, decimais = 2, grafico = TRUE) {
  hipotese <- match.arg(hipotese)
  if (!is.numeric(x) || !is.numeric(y)) stop("'x' e 'y' devem ser vetores num\u00e9ricos.")
  if (length(x) != length(y)) stop("'x' e 'y' devem ter o mesmo comprimento.")

  x_expr <- deparse(substitute(x)); x_nome <- sub(".*\\$", "", x_expr)
  y_expr <- deparse(substitute(y)); y_nome <- sub(".*\\$", "", y_expr)
  if (identical(x_nome, "x")) x_nome <- "Antes"
  if (identical(y_nome, "y")) y_nome <- "Depois"

  alt_r <- switch(hipotese, bilateral = "two.sided", maior = "greater", menor = "less")
  res   <- t.test(x, y, paired = TRUE, alternative = alt_r, conf.level = confianca)

  t_stat  <- res$statistic
  gl      <- res$parameter
  p_val   <- res$p.value
  ic      <- res$conf.int
  p_ast   <- .th_asterisk(p_val)
  dif     <- x - y
  dif_c   <- dif[!is.na(dif)]
  media_dif <- mean(dif_c)

  lbs   <- .th_hipotese_labels(hipotese, param = "μD", ref = "0")
  w_sep <- 70

  .print_titulo("TESTE T PAREADO")
  cat(sprintf("  Par: %s (antes) vs %s (depois)   |   N pares: %d\n\n",
              x_nome, y_nome, length(dif_c)))

  cat(sprintf("  H0: μD = 0  (a diferen\u00e7a m\u00e9dia \u00e9 zero)\n"))
  cat(sprintf("  H1: μD %s 0  (%s)\n\n",
              switch(hipotese, bilateral = "\u2260", maior = ">", menor = "<"), hipotese))

  # Tabela descritiva
  wg <- c(rot = 12, n = 8, med = 14, dp = 12)
  hdr_g <- paste0(.th_pad("Momento", wg["rot"], "left"),
                  .th_pad("N", wg["n"], "center"),
                  .th_pad("M\u00e9dia", wg["med"], "center"),
                  .th_pad("DP", wg["dp"], "center"))
  cat(.th_sep(nchar(hdr_g)), "\n")
  cat("  ", hdr_g, "\n", sep = "")
  cat(.th_sep(nchar(hdr_g)), "\n")

  x_c <- x[!is.na(x)]; y_c <- y[!is.na(y)]
  for (lst in list(list(x_c, x_nome), list(y_c, y_nome), list(dif_c, "Diferença (x-y)"))) {
    vals <- lst[[1]]; rot <- lst[[2]]
    cat("  ",
        .th_pad(rot, wg["rot"], "left"),
        .th_pad(length(vals), wg["n"], "center"),
        .th_pad(formatC(mean(vals), format="f", digits=decimais, decimal.mark=","), wg["med"], "center"),
        .th_pad(formatC(sd(vals),   format="f", digits=decimais, decimal.mark=","), wg["dp"],  "center"),
        "\n", sep = "")
  }
  cat(.th_sep(nchar(hdr_g)), "\n\n")

  # Tabela de teste
  ic_str <- sprintf("[%s; %s]",
                    formatC(ic[1], format="f", digits=decimais, decimal.mark=","),
                    formatC(ic[2], format="f", digits=decimais, decimal.mark=","))
  wc <- c(dif = 18, t = 14, gl = 12, ic = 20, p = 16)
  hdr_c <- paste0(
    .th_pad("M\u00e9dia das Dif.", wc["dif"], "center"),
    .th_pad("Estat\u00edstica t",  wc["t"],   "center"),
    .th_pad("Graus Lib.",    wc["gl"],  "center"),
    .th_pad(sprintf("IC (%.0f%%)", confianca*100), wc["ic"], "center"),
    .th_pad("p-valor",       wc["p"],   "center")
  )
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ",
      .th_pad(formatC(media_dif,          format="f", digits=decimais, decimal.mark=","), wc["dif"], "center"),
      .th_pad(formatC(as.numeric(t_stat), format="f", digits=3,        decimal.mark=","), wc["t"],   "center"),
      .th_pad(formatC(as.numeric(gl),     format="f", digits=0,        decimal.mark=","), wc["gl"],  "center"),
      .th_pad(ic_str, wc["ic"], "center"),
      .th_pad(paste(.th_fmt_p(p_val), p_ast), wc["p"], "center"),
      "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  .th_nota_sig(p_val)
  cat("\n")

  # Interpretação
  alpha <- 1 - confianca
  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_val < alpha) {
    cat(sprintf("  H\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  A diferen\u00e7a m\u00e9dia entre %s e %s (%s) \u00e9\n",
                x_nome, y_nome, formatC(media_dif, format="f", digits=decimais, decimal.mark=",")))
    cat(sprintf("  estatisticamente diferente de zero ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n",
                alpha * 100))
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  A diferen\u00e7a m\u00e9dia entre %s e %s (%s) n\u00e3o difere\n",
                x_nome, y_nome, formatC(media_dif, format="f", digits=decimais, decimal.mark=",")))
    cat(sprintf("  significativamente de zero ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n",
                alpha * 100))
  }
  .print_rodape()

  # Gráfico
  if (grafico) {
    tryCatch({
      if (exists("meu_tema")) {
        n_pairs <- length(x_c)
        df_plot <- data.frame(
          id     = rep(seq_len(n_pairs), 2),
          momento = factor(c(rep(x_nome, n_pairs), rep(y_nome, n_pairs)), levels = c(x_nome, y_nome)),
          valor  = c(x_c, y_c)
        )

        p_sub <- sprintf("p %s %s%s",
                         ifelse(p_val < 0.001, "<", "="),
                         .th_fmt_p(p_val),
                         ifelse(nchar(p_ast) > 0, paste0(" ", p_ast), ""))

        old_w <- getOption("warn"); options(warn = -1)
        p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = momento, y = valor, group = id)) +
          ggplot2::geom_line(color = "gray70", linewidth = 0.5) +
          ggplot2::geom_point(ggplot2::aes(color = momento), size = 2.5, alpha = 0.8) +
          ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "line",
                                color = "#D90429", linewidth = 1.2, linetype = "dashed") +
          ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "point",
                                color = "#D90429", size = 4, shape = 18) +
          ggplot2::scale_color_manual(values = c("#555555", "#1B4F72")) +
          ggplot2::labs(
            title    = sprintf("Teste Pareado: %s vs %s", x_nome, y_nome),
            subtitle = p_sub,
            x = "Momento", y = "Valor", color = NULL
          ) +
          meu_tema()
        suppressMessages(print(p_plot))
        options(warn = old_w)
      }
    }, error = function(e) {
      message("[Aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico: ", e$message)
    })
  }

  invisible(res)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teste de Proporção (Z)
#' @description Testa se uma proporção amostral é igual a um valor de referência.
#' @param x Número de sucessos observados.
#' @param n Tamanho total da amostra.
#' @param p0 Proporção de referência (entre 0 e 1).
#' @param hipotese Direção do teste: \code{"bilateral"} (padrão), \code{"maior"} ou \code{"menor"}.
#' @param confianca Nível de confiança para o IC (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @param grafico Lógico. Se TRUE, exibe a distribuição normal com a região crítica.
#' @return Retorna invisivelmente o objeto do teste (\code{htest}).
#' @export
teste_proporcao <- function(x, n, p0, hipotese = c("bilateral", "maior", "menor"),
                             confianca = 0.95, decimais = 2, grafico = TRUE) {
  hipotese <- match.arg(hipotese)
  if (x > n) stop("'x' (sucessos) n\u00e3o pode ser maior que 'n' (total).")
  if (p0 <= 0 || p0 >= 1) stop("'p0' deve estar entre 0 e 1.")

  alt_r <- switch(hipotese, bilateral = "two.sided", maior = "greater", menor = "less")
  res   <- prop.test(x, n, p = p0, alternative = alt_r, conf.level = confianca, correct = FALSE)

  p_hat  <- x / n
  z_stat <- (p_hat - p0) / sqrt(p0 * (1 - p0) / n)
  p_val  <- res$p.value
  ic     <- res$conf.int
  p_ast  <- .th_asterisk(p_val)
  w_sep  <- 70

  .print_titulo("TESTE DE PROPOR\u00c7\u00c3O (Z)")
  cat(sprintf("  Sucessos: %d de %d   |   p0: %.1f%%\n\n", x, n, p0 * 100))

  h0_txt <- sprintf("H0: p = %.2f", p0)
  h1_op  <- switch(hipotese, bilateral = "\u2260", maior = ">", menor = "<")
  h1_txt <- sprintf("H1: p %s %.2f (%s)", h1_op, p0, hipotese)
  cat(sprintf("  %s\n  %s\n\n", h0_txt, h1_txt))

  # Tabela de resultados
  ic_str <- sprintf("[%s%%; %s%%]",
                    formatC(ic[1]*100, format="f", digits=1, decimal.mark=","),
                    formatC(ic[2]*100, format="f", digits=1, decimal.mark=","))
  wc <- c(ph = 20, z = 14, ic = 20, p = 16)
  hdr_c <- paste0(
    .th_pad("Prop. amostral", wc["ph"], "center"),
    .th_pad("Estat\u00edstica Z",  wc["z"],  "center"),
    .th_pad(sprintf("IC (%.0f%%)", confianca*100), wc["ic"], "center"),
    .th_pad("p-valor",         wc["p"],  "center")
  )
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ",
      .th_pad(sprintf("%s%%", formatC(p_hat*100, format="f", digits=1, decimal.mark=",")), wc["ph"], "center"),
      .th_pad(formatC(z_stat, format="f", digits=3, decimal.mark=","), wc["z"],  "center"),
      .th_pad(ic_str, wc["ic"], "center"),
      .th_pad(paste(.th_fmt_p(p_val), p_ast), wc["p"], "center"),
      "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  .th_nota_sig(p_val)
  cat("\n")

  # Interpretação
  alpha <- 1 - confianca
  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_val < alpha) {
    cat(sprintf("  H\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  A propor\u00e7\u00e3o observada (%s%%) difere significativamente\n",
                formatC(p_hat*100, format="f", digits=1, decimal.mark=",")))
    cat(sprintf("  de %.1f%% ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n", p0*100, alpha*100))
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  A propor\u00e7\u00e3o observada (%s%%) n\u00e3o difere significativamente\n",
                formatC(p_hat*100, format="f", digits=1, decimal.mark=",")))
    cat(sprintf("  de %.1f%% ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n", p0*100, alpha*100))
  }
  .print_rodape()

  # Gráfico: distribuição normal com região crítica
  if (grafico) {
    tryCatch({
      if (exists("meu_tema")) {
        alpha_n <- 1 - confianca
        zc <- switch(hipotese,
          bilateral = qnorm(1 - alpha_n/2),
          maior     = qnorm(1 - alpha_n),
          menor     = qnorm(alpha_n)
        )
        x_seq   <- seq(-4, 4, length.out = 500)
        df_curv <- data.frame(x = x_seq, y = dnorm(x_seq))

        p_plot <- ggplot2::ggplot(df_curv, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_line(linewidth = 0.9, color = "#333333")

        if (hipotese == "bilateral") {
          p_plot <- p_plot +
            ggplot2::geom_area(data = df_curv[df_curv$x <= -abs(zc), ], fill = "red1", alpha = 0.35) +
            ggplot2::geom_area(data = df_curv[df_curv$x >= abs(zc),  ], fill = "red1", alpha = 0.35)
        } else if (hipotese == "maior") {
          p_plot <- p_plot + ggplot2::geom_area(data = df_curv[df_curv$x >= zc, ], fill = "red1", alpha = 0.35)
        } else {
          p_plot <- p_plot + ggplot2::geom_area(data = df_curv[df_curv$x <= zc, ], fill = "red1", alpha = 0.35)
        }

        p_sub <- sprintf("%s | p %s %s%s", h1_txt,
                         ifelse(p_val < 0.001, "<", "="),
                         .th_fmt_p(p_val),
                         ifelse(nchar(p_ast) > 0, paste0(" ", p_ast), ""))

        p_plot <- p_plot +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
          ggplot2::geom_vline(xintercept = z_stat, color = "#1B4F72",
                              linetype = "dashed", linewidth = 1) +
          ggplot2::annotate("text", x = z_stat, y = max(df_curv$y) * 1.05,
                            label = sprintf("Z = %s", formatC(z_stat, format="f", digits=3, decimal.mark=",")),
                            hjust = -0.1, color = "#1B4F72", fontface = "bold", size = 4.5) +
          ggplot2::labs(
            title    = "Distribuição Normal — Teste de Proporção (Z)",
            subtitle = p_sub,
            x = "Z", y = "Densidade"
          ) +
          meu_tema(grade = "dupla")

        suppressMessages(print(p_plot))
      }
    }, error = function(e) {
      message("[Aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico: ", e$message)
    })
  }

  invisible(res)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teste de Wilcoxon (Uma Amostra / Pareado)
#' @description Alternativa não-paramétrica ao teste t para uma amostra ou para amostras pareadas.
#'   Indicado quando os dados não seguem distribuição normal.
#' @param x Vetor numérico (ou diferenças, no caso pareado).
#' @param y Vetor numérico do segundo momento (opcional). Se informado, o teste é pareado.
#' @param mu Valor de referência para a mediana (padrão \code{0}).
#' @param hipotese Direção do teste: \code{"bilateral"} (padrão), \code{"maior"} ou \code{"menor"}.
#' @param confianca Nível de confiança para o IC (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @return Retorna invisivelmente o objeto do teste (\code{htest}).
#' @export
teste_wilcoxon <- function(x, y = NULL, mu = 0, hipotese = c("bilateral", "maior", "menor"), confianca = 0.95, decimais = 2, grafico = TRUE) {
  hipotese <- match.arg(hipotese)
  alt_r    <- switch(hipotese, bilateral = "two.sided", maior = "greater", menor = "less")

  x_expr <- deparse(substitute(x)); x_nome <- sub(".*\\$", "", x_expr)
  if (identical(x_nome, "x")) x_nome <- "x"

  pareado <- !is.null(y)
  if (pareado) {
    res  <- wilcox.test(x, y, paired = TRUE, alternative = alt_r,
                        conf.int = TRUE, conf.level = confianca, exact = FALSE)
    tipo <- "Teste de Wilcoxon Pareado (Signed-Rank)"
  } else {
    res  <- wilcox.test(x, mu = mu, alternative = alt_r,
                        conf.int = TRUE, conf.level = confianca, exact = FALSE)
    tipo <- "Teste de Wilcoxon para Uma Amostra (Signed-Rank)"
  }

  W      <- res$statistic
  p_val  <- res$p.value
  ic     <- res$conf.int
  p_ast  <- .th_asterisk(p_val)
  w_sep  <- 70

  .print_titulo("TESTE N\u00c3O-PARAM\u00c9TRICO DE WILCOXON")
  cat(sprintf("  Tipo: %s\n", tipo))
  cat(sprintf("  Nota: Alternativa n\u00e3o-param\u00e9trica ao teste t. Testa a mediana.\n\n"))

  h0_txt <- if (pareado) "H0: a mediana das diferenças = 0" else sprintf("H0: mediana = %s", mu)
  h1_op  <- switch(hipotese, bilateral = "\u2260", maior = ">", menor = "<")
  h1_txt <- if (pareado) sprintf("H1: mediana das diferenças %s 0 (%s)", h1_op, hipotese) else
                         sprintf("H1: mediana %s %s (%s)", h1_op, mu, hipotese)
  cat(sprintf("  %s\n  %s\n\n", h0_txt, h1_txt))

  ic_str <- sprintf("[%s; %s]",
                    formatC(ic[1], format="f", digits=decimais, decimal.mark=","),
                    formatC(ic[2], format="f", digits=decimais, decimal.mark=","))
  wc <- c(w = 14, ic = 22, p = 16)
  hdr_c <- paste0(.th_pad("Estat. W", wc["w"], "center"),
                  .th_pad(sprintf("IC (%.0f%%)", confianca*100), wc["ic"], "center"),
                  .th_pad("p-valor", wc["p"], "center"))
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ",
      .th_pad(formatC(as.numeric(W), format="f", digits=0, decimal.mark=","), wc["w"],  "center"),
      .th_pad(ic_str, wc["ic"], "center"),
      .th_pad(paste(.th_fmt_p(p_val), p_ast), wc["p"], "center"),
      "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  .th_nota_sig(p_val)
  cat("\n")

  alpha <- 1 - confianca
  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_val < alpha) {
    cat(sprintf("  H\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    if (pareado) cat("  A mediana das diferenças difere significativamente de zero.\n")
    else cat(sprintf("  A mediana de %s difere significativamente de %s.\n", x_nome, mu))
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    if (pareado) cat("  A mediana das diferenças n\u00e3o difere significativamente de zero.\n")
    else cat(sprintf("  A mediana de %s n\u00e3o difere significativamente de %s.\n", x_nome, mu))
  }
  .print_rodape()
  if (grafico) {
    tryCatch({
      if (exists("meu_tema")) {
        p_sub <- sprintf("%s | p %s %s%s", h1_txt,
                         ifelse(p_val < 0.001, "<", "="),
                         .th_fmt_p(p_val),
                         ifelse(nchar(p_ast) > 0, paste0(" ", p_ast), ""))
        if (!pareado) {
          x_limpo <- x[!is.na(x)]
          df_w <- data.frame(val = x_limpo)
          med_amostra <- median(x_limpo)
          dens_res <- density(x_limpo)
          max_d <- max(dens_res$y, na.rm = TRUE)
          y_txt1 <- max_d * 0.90
          y_txt2 <- max_d * 0.78

          p_plot <- ggplot2::ggplot(df_w, ggplot2::aes(x = val)) +
            ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 15, fill = "#DDDDDD", color = "black", alpha = 0.85) +
            ggplot2::geom_density(color = "#0072CE", linewidth = 1.1) +
            ggplot2::geom_vline(xintercept = med_amostra, color = "#0072CE", linetype = "solid", linewidth = 1) +
            ggplot2::geom_vline(xintercept = mu, color = "red3", linetype = "dashed", linewidth = 1) +
            ggplot2::annotate("text", x = med_amostra, y = y_txt1, label = sprintf("Mediana: %s", formatC(med_amostra, format="f", digits=decimais, decimal.mark=",")), hjust = -0.1, color = "#0072CE", fontface = "bold", size = 3.8) +
            ggplot2::annotate("text", x = mu, y = y_txt2, label = sprintf("Referência: %s", mu), hjust = 1.1, color = "red3", fontface = "bold", size = 3.8) +
            ggplot2::labs(
              title = sprintf("Teste de Wilcoxon: %s", x_nome),
              subtitle = p_sub,
              x = x_nome, y = "Densidade"
            ) +
            meu_tema(estilo = 1, grade = "dupla")
          suppressMessages(print(p_plot))
        } else {
          dif_w <- (x - y)[!is.na(x - y)]
          df_w <- data.frame(dif = dif_w)
          med_dif <- median(dif_w)
          p_plot <- ggplot2::ggplot(df_w, ggplot2::aes(x = dif)) +
            ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 15, fill = "#85b4d4", color = "white", alpha = 0.8) +
            ggplot2::geom_density(color = "#1B4F72", linewidth = 1) +
            ggplot2::geom_vline(xintercept = med_dif, color = "#1B4F72", linetype = "solid", linewidth = 1) +
            ggplot2::geom_vline(xintercept = 0, color = "red3", linetype = "dashed", linewidth = 1) +
            ggplot2::labs(
              title = "Teste de Wilcoxon Pareado (Diferenças)",
              subtitle = p_sub,
              x = "Diferença (x - y)", y = "Densidade"
            ) +
            meu_tema(grade = "dupla")
          suppressMessages(print(p_plot))
        }
      }
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }

  invisible(res)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teste de Mann-Whitney (Wilcoxon para Duas Amostras Independentes)
#' @description Alternativa não-paramétrica ao teste t para dois grupos independentes.
#'   Indicado quando os dados não seguem distribuição normal.
#' @param x Vetor numérico com os valores.
#' @param grupo Vetor com os grupos (deve ter exatamente 2 níveis).
#' @param hipotese Direção do teste: \code{"bilateral"} (padrão), \code{"maior"} ou \code{"menor"}.
#' @param confianca Nível de confiança para o IC (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @return Retorna invisivelmente o objeto do teste (\code{htest}).
#' @export
teste_mann_whitney <- function(x, grupo, hipotese = c("bilateral", "maior", "menor"), confianca = 0.95, decimais = 2, grafico = TRUE) {
  hipotese <- match.arg(hipotese)

  grupo_sub  <- substitute(grupo)
  grupo_nome <- sub(".*\\$", "", deparse(grupo_sub))
  x_expr     <- deparse(substitute(x))
  x_nome     <- sub(".*\\$", "", x_expr)

  grupo_vec <- as.factor(eval(grupo_sub, parent.frame()))
  niveis    <- levels(grupo_vec)
  if (length(niveis) != 2) stop("'grupo' deve ter exatamente 2 n\u00edveis.")

  g1 <- niveis[1]; g2 <- niveis[2]
  vals1 <- x[grupo_vec == g1 & !is.na(grupo_vec)]; vals1 <- vals1[!is.na(vals1)]
  vals2 <- x[grupo_vec == g2 & !is.na(grupo_vec)]; vals2 <- vals2[!is.na(vals2)]

  alt_r <- switch(hipotese, bilateral = "two.sided", maior = "greater", menor = "less")
  res   <- wilcox.test(vals1, vals2, alternative = alt_r,
                       conf.int = TRUE, conf.level = confianca, exact = FALSE)

  W     <- res$statistic
  p_val <- res$p.value
  ic    <- res$conf.int
  p_ast <- .th_asterisk(p_val)
  w_sep <- 70

  .print_titulo("TESTE DE MANN-WHITNEY")
  cat(sprintf("  Vari\u00e1vel: %s   |   Grupos: %s vs %s\n", x_nome, g1, g2))
  cat("  Nota: Alternativa n\u00e3o-param\u00e9trica ao teste t para grupos independentes.\n\n")

  h1_op  <- switch(hipotese, bilateral = "\u2260", maior = ">", menor = "<")
  cat(sprintf("  H0: distribui\u00e7\u00e3o de %s = distribui\u00e7\u00e3o de %s\n", g1, g2))
  cat(sprintf("  H1: distribui\u00e7\u00e3o de %s %s distribui\u00e7\u00e3o de %s (%s)\n\n", g1, h1_op, g2, hipotese))

  # Descritiva
  wg <- c(grp=14, n=8, med=14, mediana=14)
  hdr_g <- paste0(.th_pad("Grupo",   wg["grp"],    "left"),
                  .th_pad("N",       wg["n"],      "center"),
                  .th_pad("M\u00e9dia",  wg["med"],    "center"),
                  .th_pad("Mediana", wg["mediana"], "center"))
  cat(.th_sep(nchar(hdr_g)), "\n")
  cat("  ", hdr_g, "\n", sep = "")
  cat(.th_sep(nchar(hdr_g)), "\n")
  for (lst in list(list(vals1, g1), list(vals2, g2))) {
    vals <- lst[[1]]; rot <- lst[[2]]
    cat("  ",
        .th_pad(rot, wg["grp"], "left"),
        .th_pad(length(vals), wg["n"], "center"),
        .th_pad(formatC(mean(vals),   format="f", digits=decimais, decimal.mark=","), wg["med"],    "center"),
        .th_pad(formatC(median(vals), format="f", digits=decimais, decimal.mark=","), wg["mediana"],"center"),
        "\n", sep = "")
  }
  cat(.th_sep(nchar(hdr_g)), "\n\n")

  ic_str <- sprintf("[%s; %s]",
                    formatC(ic[1], format="f", digits=decimais, decimal.mark=","),
                    formatC(ic[2], format="f", digits=decimais, decimal.mark=","))
  wc <- c(w=14, ic=22, p=16)
  hdr_c <- paste0(.th_pad("Estat. W", wc["w"], "center"),
                  .th_pad(sprintf("IC (%.0f%%)", confianca*100), wc["ic"], "center"),
                  .th_pad("p-valor", wc["p"], "center"))
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ",
      .th_pad(formatC(as.numeric(W), format="f", digits=0, decimal.mark=","), wc["w"],  "center"),
      .th_pad(ic_str, wc["ic"], "center"),
      .th_pad(paste(.th_fmt_p(p_val), p_ast), wc["p"], "center"),
      "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  .th_nota_sig(p_val)
  cat("\n")

  alpha <- 1 - confianca
  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_val < alpha) {
    cat(sprintf("  H\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  As distribui\u00e7\u00f5es de %s e %s diferem significativamente\n", g1, g2))
    cat(sprintf("  ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n", alpha * 100))
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  N\u00e3o h\u00e1 diferen\u00e7a significativa entre as distribui\u00e7\u00f5es\n"))
    cat(sprintf("  de %s e %s ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n", g1, g2, alpha * 100))
  }
  .print_rodape()
  if (grafico) {
    tryCatch({
      if (exists("meu_tema")) {
        p_sub <- sprintf("H1: %s %s %s (%s) | p %s %s%s", g1, h1_op, g2, hipotese,
                         ifelse(p_val < 0.001, "<", "="),
                         .th_fmt_p(p_val),
                         ifelse(nchar(p_ast) > 0, paste0(" ", p_ast), ""))
        df_plot <- data.frame(grupo = grupo_vec, y = x)
        df_plot <- df_plot[!is.na(df_plot$grupo) & !is.na(df_plot$y), ]
        max_y <- max(df_plot$y, na.rm = TRUE)
        min_y <- min(df_plot$y, na.rm = TRUE)
        amp <- max_y - min_y
        y_bar <- max_y + amp * 0.08
        y_tick <- max_y + amp * 0.05
        y_text <- max_y + amp * 0.12

        old_w <- getOption("warn"); options(warn = -1)
        med1 <- median(vals1); med2 <- median(vals2)
        cor1 <- if (med1 >= med2) "#0072CE" else "#DDDDDD"
        cor2 <- if (med2 > med1) "#0072CE" else "#DDDDDD"
        cores_box <- setNames(c(cor1, cor2), c(g1, g2))

        p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = grupo, y = y, fill = grupo)) +
          ggplot2::stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.6) +
          ggplot2::geom_boxplot(width = 0.5, alpha = 0.85, color = "black", outlier.color = "red3", outlier.shape = 19) +
          ggplot2::stat_summary(fun = median, geom = "point", shape = 21, size = 3.2, fill = "black", color = "black") +
          ggplot2::annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, color = "black", linewidth = 0.5) +
          ggplot2::annotate("segment", x = 1, xend = 1, y = y_tick, yend = y_bar, color = "black", linewidth = 0.5) +
          ggplot2::annotate("segment", x = 2, xend = 2, y = y_tick, yend = y_bar, color = "black", linewidth = 0.5) +
          ggplot2::annotate("text", x = 1.5, y = y_text, label = p_sub, size = 4.2, fontface = "bold") +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
          ggplot2::scale_fill_manual(values = c("#3470a3", "#85b4d4")) +
          ggplot2::labs(
            title    = sprintf("Teste de Mann-Whitney: %s por %s", x_nome, grupo_nome),
            x        = grupo_nome,
            y        = x_nome
          ) +
          meu_tema(grade = "dupla") +
          ggplot2::theme(legend.position = "none")

        suppressMessages(print(p_plot))
        options(warn = old_w)
      }
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }

  invisible(res)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teste Qui-Quadrado de Independência
#' @description Testa se duas variáveis categóricas são independentes.
#' @param dados Data frame com os dados.
#' @param var_x Nome da primeira variável categórica (sem aspas).
#' @param var_y Nome da segunda variável categórica (sem aspas).
#' @param confianca Nível de confiança (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @return Retorna invisivelmente o objeto do teste (\code{htest}).
#' @export
teste_qui_quadrado <- function(dados, var_x, var_y, confianca = 0.95, decimais = 2, grafico = TRUE) {
  vx_sub  <- substitute(var_x)
  vy_sub  <- substitute(var_y)
  vx_nome <- as.character(vx_sub)
  vy_nome <- as.character(vy_sub)
  if (length(vx_nome) > 1 || !vx_nome %in% names(dados)) {
    vx_eval <- tryCatch(eval(vx_sub, parent.frame()), error = function(e) NULL)
    if (is.character(vx_eval) && length(vx_eval) == 1 && vx_eval %in% names(dados)) vx_nome <- vx_eval
  }
  if (length(vy_nome) > 1 || !vy_nome %in% names(dados)) {
    vy_eval <- tryCatch(eval(vy_sub, parent.frame()), error = function(e) NULL)
    if (is.character(vy_eval) && length(vy_eval) == 1 && vy_eval %in% names(dados)) vy_nome <- vy_eval
  }

  if (!vx_nome %in% names(dados)) stop(sprintf("Coluna '%s' n\u00e3o encontrada.", vx_nome))
  if (!vy_nome %in% names(dados)) stop(sprintf("Coluna '%s' n\u00e3o encontrada.", vy_nome))

  tab   <- table(dados[[vx_nome]], dados[[vy_nome]])
  res   <- suppressWarnings(chisq.test(tab, correct = FALSE))
  tem_esperado_baixo <- any(res$expected < 5)
  chi2  <- res$statistic
  gl    <- res$parameter
  p_val <- res$p.value
  p_ast <- .th_asterisk(p_val)
  w_sep <- 70

  .print_titulo("TESTE QUI-QUADRADO DE INDEPEND\u00caNCIA")
  cat(sprintf("  Vari\u00e1veis: %s  \u00d7  %s\n\n", vx_nome, vy_nome))
  cat(sprintf("  H0: %s e %s s\u00e3o independentes\n", vx_nome, vy_nome))
  cat(sprintf("  H1: %s e %s n\u00e3o s\u00e3o independentes\n\n", vx_nome, vy_nome))

  # Tabela de contingência
  cat("  Tabela de Conting\u00eancia (frequências observadas):\n\n")
  print(tab)
  cat("\n")

  wc <- c(chi = 16, gl = 12, p = 16)
  hdr_c <- paste0(.th_pad("\u03c7\u00b2", wc["chi"], "center"),
                  .th_pad("Graus Lib.", wc["gl"],  "center"),
                  .th_pad("p-valor",   wc["p"],   "center"))
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ",
      .th_pad(formatC(as.numeric(chi2), format="f", digits=3, decimal.mark=","), wc["chi"], "center"),
      .th_pad(as.integer(gl), wc["gl"], "center"),
      .th_pad(paste(.th_fmt_p(p_val), p_ast), wc["p"], "center"),
      "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  .th_nota_sig(p_val)
  cat("\n")

  alpha <- 1 - confianca
  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_val < alpha) {
    cat(sprintf("  H\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  Existe associa\u00e7\u00e3o estatisticamente significativa entre\n"))
    cat(sprintf("  %s e %s ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n",
                vx_nome, vy_nome, alpha * 100))
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias suficientes para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat(sprintf("  N\u00e3o h\u00e1 associa\u00e7\u00e3o significativa entre %s e %s\n", vx_nome, vy_nome))
    cat(sprintf("  ao n\u00edvel de %.0f%% de signific\u00e2ncia.\n", alpha * 100))
  }
  .print_rodape()
  if (grafico && exists("meu_tema")) {
    tryCatch({
      df_qui <- as.data.frame(tab)
      colnames(df_qui) <- c("X", "Y", "Freq")
      p_qui <- ggplot2::ggplot(df_qui, ggplot2::aes(x = X, y = Freq, fill = Y)) +
        ggplot2::geom_col(position = "fill", color = "black", width = 0.55) +
        ggplot2::scale_y_continuous(labels = scales::percent, expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::labs(
          title = sprintf("Associação: %s x %s", vx_nome, vy_nome),
          subtitle = sprintf("Qui-Quadrado = %.3f | p-valor %s", chi2, if (p_val < 0.001) "< 0,001" else paste0("= ", formatC(p_val, format="f", digits=3, decimal.mark=","))),
          x = vx_nome, y = "Proporção", fill = vy_nome
        ) +
        meu_tema(estilo = 1, grade = "dupla")
      print(p_qui)
    }, error = function(e) NULL)
  }
  invisible(res)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teste de Levene (Homogeneidade de Variâncias)
#' @description Testa se dois grupos possuem variâncias iguais (homogeneidade).
#'   Auxiliar ao teste t para duas amostras.
#' @param x Vetor numérico com os valores.
#' @param grupo Vetor com os grupos (deve ter exatamente 2 níveis).
#' @param confianca Nível de confiança (padrão \code{0.95}).
#' @param decimais Casas decimais para exibição (padrão \code{2}).
#' @return Retorna invisivelmente uma lista com a estatística F e p-valor.
#' @export
teste_levene <- function(x, grupo, confianca = 0.95, decimais = 2, grafico = TRUE) {
  grupo_sub  <- substitute(grupo)
  grupo_nome <- sub(".*\\$", "", deparse(grupo_sub))
  x_expr     <- deparse(substitute(x))
  x_nome     <- sub(".*\\$", "", x_expr)

  grupo_vec <- as.factor(eval(grupo_sub, parent.frame()))
  niveis    <- levels(grupo_vec)
  if (length(niveis) < 2) stop("'grupo' deve ter pelo menos 2 n\u00edveis.")

  # Levene em base R: usa desvios em relação à mediana de cada grupo
  medianas  <- tapply(x, grupo_vec, median, na.rm = TRUE)
  z         <- abs(x - medianas[grupo_vec])
  res_lev   <- summary(aov(z ~ grupo_vec))
  f_stat    <- res_lev[[1]]["grupo_vec", "F value"]
  p_val     <- res_lev[[1]]["grupo_vec", "Pr(>F)"]
  p_ast     <- .th_asterisk(p_val)
  w_sep     <- 70

  .print_titulo("TESTE DE LEVENE (HOMOGENEIDADE DE VARI\u00c2NCIAS)")
  cat(sprintf("  Vari\u00e1vel: %s   |   Grupos: %s\n\n", x_nome, grupo_nome))
  cat("  H0: as vari\u00e2ncias s\u00e3o iguais entre os grupos\n")
  cat("  H1: pelo menos uma vari\u00e2ncia difere\n\n")

  # Variâncias por grupo
  for (niv in niveis) {
    vals_niv <- x[grupo_vec == niv & !is.na(grupo_vec)]
    cat(sprintf("  Vari\u00e2ncia %s: %s\n", niv,
                formatC(var(vals_niv, na.rm=TRUE), format="f", digits=decimais, decimal.mark=",")))
  }
  cat("\n")

  wc <- c(f=14, p=16)
  hdr_c <- paste0(.th_pad("Estat. F", wc["f"], "center"),
                  .th_pad("p-valor",  wc["p"],  "center"))
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  cat("  ",
      .th_pad(formatC(f_stat, format="f", digits=3, decimal.mark=","), wc["f"], "center"),
      .th_pad(paste(.th_fmt_p(p_val), p_ast), wc["p"], "center"),
      "\n", sep = "")
  cat(.th_sep(nchar(hdr_c)), "\n")
  .th_nota_sig(p_val)
  cat("\n")

  alpha <- 1 - confianca
  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_val < alpha) {
    cat(sprintf("  H\u00e1 evid\u00eancias para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat("  As vari\u00e2ncias dos grupos s\u00e3o significativamente diferentes.\n")
    cat("  Recomenda-se usar o teste t de Welch (variancia_igual = FALSE).\n")
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias para rejeitar H0 (%s).\n",
                .th_fmt_p_interp(p_val)))
    cat("  As vari\u00e2ncias dos grupos s\u00e3o homog\u00eaneas.\n")
    cat("  O teste t de Student (variancia_igual = TRUE) pode ser utilizado.\n")
  }
  .print_rodape()
  if (grafico && exists("meu_tema")) {
    tryCatch({
      df_lev <- data.frame(grupo = grupo_vec, z = z)
      df_lev <- df_lev[!is.na(df_lev$grupo) & !is.na(df_lev$z), ]
      max_z_grp <- tapply(df_lev$z, df_lev$grupo, median, na.rm = TRUE)
      maior_grp <- names(which.max(max_z_grp))
      cores_lev <- setNames(ifelse(levels(grupo_vec) == maior_grp, "#0072CE", "#DDDDDD"), levels(grupo_vec))

      p_lev <- ggplot2::ggplot(df_lev, ggplot2::aes(x = grupo, y = z, fill = grupo)) +
        ggplot2::stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.6) +
        ggplot2::geom_boxplot(width = 0.5, alpha = 0.85, color = "black", outlier.color = "red3") +
        ggplot2::stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "black", color = "black") +
        ggplot2::scale_fill_manual(values = cores_lev) +
        ggplot2::labs(
          title = "Teste de Levene: Desvios da Mediana por Grupo",
          subtitle = sprintf("F = %.3f | p-valor %s", f_stat, if (p_val < 0.001) "< 0,001" else paste0("= ", formatC(p_val, format="f", digits=3, decimal.mark=","))),
          x = grupo_nome, y = "|y - Mediana|"
        ) +
        meu_tema(estilo = 2, grade = "dupla") +
        ggplot2::theme(legend.position = "none")
      print(p_lev)
    }, error = function(e) NULL)
  }
  invisible(list(F = f_stat, p.value = p_val))
}
