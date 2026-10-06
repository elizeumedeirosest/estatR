# ── Modulo: distribuicao() ────────────────────────────────────────────────────
# Visao geral pedagogica de uma distribuicao: parametros, momentos,
# quantis, regra empirica (Normal), grafico pedagogico ou painel de simulacao.

#' @title Distribuição Estatística
#' @description Exibe uma visão completa e pedagógica de uma distribuição de probabilidade:
#' parâmetros, momentos teóricos, quantis-chave e gráfico ilustrativo.
#' Suporta: "normal", "binomial", "poisson", "exponencial".
#'
#' @param tipo Tipo da distribuição: "normal", "binomial", "poisson" ou "exponencial".
#' @param media Parâmetro da Normal: média (μ). Padrão: 0.
#' @param dp Parâmetro da Normal: desvio-padrão (σ). Padrão: 1.
#' @param n Parâmetro da Binomial: número de ensaios. Padrão: 10.
#' @param p Parâmetro da Binomial: probabilidade de sucesso. Padrão: 0.5.
#' @param lambda Parâmetro da Poisson: taxa de ocorrência (λ). Padrão: 1.
#' @param taxa Parâmetro da Exponencial: taxa (λ). Padrão: 1.
#' @param grafico Lógico. Se TRUE, exibe o gráfico pedagógico da distribuição.
#' @param simulacao Lógico. Se TRUE, exibe painel com histogramas simulados
#'   para n = 10, 50, 100, 250, 500 e 1000 observações.
#' @param semente Semente para reprodutibilidade da simulação.
#' @export
distribuicao <- function(tipo = c("normal", "binomial", "poisson", "exponencial"),
                         media = 0, dp = 1, n = 10, p = 0.5, lambda = 1, taxa = 1,
                         grafico = TRUE, simulacao = FALSE, semente = 42) {

  tipo <- match.arg(tipo)

  # ── Helpers internos ──────────────────────────────────────────────────────
  .pad <- function(s, w, align = "left") {
    if (length(s) > 1) return(sapply(s, .pad, w = w, align = align))
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left")  return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .f2 <- function(x) formatC(x, format = "f", digits = 2, decimal.mark = ",")
  .f3 <- function(x) formatC(x, format = "f", digits = 4, decimal.mark = ",")

  # ── Calculos por tipo ─────────────────────────────────────────────────────
  if (tipo == "normal") {
    if (dp <= 0) stop("O desvio-padrão (dp) deve ser > 0.")
    teo_media  <- media
    teo_var    <- dp^2
    teo_dp     <- dp
    teo_assim  <- 0
    teo_curt   <- 0
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O NORMAL  \u2014  \u03bc = %g  |  \u03c3 = %g", media, dp)
    params_txt <- sprintf("  M\u00e9dia (\u03bc): %g    Desvio Padr\u00e3o (\u03c3): %g\n", media, dp)
    quantis <- qnorm(c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
                     mean = media, sd = dp)
    regra68  <- c(media - dp,  media + dp)
    regra95  <- c(media - 2*dp, media + 2*dp)
    regra997 <- c(media - 3*dp, media + 3*dp)
    tem_regra <- TRUE
    pdf_formula <- sprintf("f(x) = (1 / (\u03c3\u221a2\u03c0)) \u00d7 exp(-\u00bd((x-\u03bc)/\u03c3)\u00b2)")

  } else if (tipo == "binomial") {
    if (n <= 0 || n != floor(n)) stop("'n' deve ser inteiro positivo.")
    if (p < 0 || p > 1)          stop("'p' deve estar entre 0 e 1.")
    teo_media  <- n * p
    teo_var    <- n * p * (1 - p)
    teo_dp     <- sqrt(teo_var)
    teo_assim  <- (1 - 2*p) / sqrt(n * p * (1 - p))
    teo_curt   <- (1 - 6*p*(1-p)) / (n * p * (1-p))
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O BINOMIAL  \u2014  n = %d  |  p = %g", n, p)
    params_txt <- sprintf("  Ensaios (n): %d    Probabilidade (p): %g\n", n, p)
    quantis <- qbinom(c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99), n, p)
    tem_regra <- FALSE
    pdf_formula <- sprintf("P(X=k) = C(n,k) \u00d7 p\u1d4f \u00d7 (1-p)^(n-k)")

  } else if (tipo == "poisson") {
    if (lambda <= 0) stop("'lambda' deve ser > 0.")
    teo_media  <- lambda
    teo_var    <- lambda
    teo_dp     <- sqrt(lambda)
    teo_assim  <- 1 / sqrt(lambda)
    teo_curt   <- 1 / lambda
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O POISSON  \u2014  \u03bb = %g", lambda)
    params_txt <- sprintf("  Taxa (\u03bb): %g\n", lambda)
    quantis <- qpois(c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99), lambda)
    tem_regra <- FALSE
    pdf_formula <- sprintf("P(X=k) = (e^-\u03bb \u00d7 \u03bb\u1d4f) / k!")

  } else if (tipo == "exponencial") {
    if (taxa <= 0) stop("'taxa' deve ser > 0.")
    teo_media  <- 1 / taxa
    teo_var    <- 1 / taxa^2
    teo_dp     <- 1 / taxa
    teo_assim  <- 2
    teo_curt   <- 6
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O EXPONENCIAL  \u2014  \u03bb = %g", taxa)
    params_txt <- sprintf("  Taxa (\u03bb): %g    M\u00e9dia (1/\u03bb): %g\n", taxa, 1/taxa)
    quantis <- qexp(c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99), rate = taxa)
    tem_regra <- FALSE
    pdf_formula <- sprintf("f(x) = \u03bb \u00d7 e^(-\u03bb\u00d7x),  x \u2265 0")
  }

  # ── Saida Console ─────────────────────────────────────────────────────────
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo(titulo_dist)
  } else {
    cat(sprintf("\n\u2500\u2500 %s \u2500\u2500\n", titulo_dist))
  }

  if (exists(".print_topico", mode = "function")) .print_topico("PAR\u00c2METROS")
  cat(params_txt)
  cat(sprintf("  F\u00f3rmula:  %s\n\n", pdf_formula))

  if (exists(".print_topico", mode = "function")) .print_topico("MOMENTOS TE\u00d3RICOS")
  wm <- 22; wv <- 12
  sep_m <- paste0("  ", strrep("\u2500", wm + wv + 1))
  cat(sep_m, "\n")
  cat("  ", .pad("Medida", wm, "left"), .pad("Valor", wv, "center"), "\n", sep = "")
  cat(sep_m, "\n")
  cat("  ", .pad("M\u00e9dia", wm, "left"),      .pad(.f3(teo_media), wv, "center"), "\n", sep = "")
  cat("  ", .pad("Vari\u00e2ncia", wm, "left"),  .pad(.f3(teo_var),   wv, "center"), "\n", sep = "")
  cat("  ", .pad("Desvio Padr\u00e3o", wm, "left"), .pad(.f3(teo_dp), wv, "center"), "\n", sep = "")
  cat("  ", .pad("Assimetria", wm, "left"),    .pad(.f3(teo_assim), wv, "center"), "\n", sep = "")
  cat("  ", .pad("Curtose (excesso)", wm, "left"), .pad(.f3(teo_curt), wv, "center"), "\n", sep = "")
  cat(sep_m, "\n\n")

  if (tem_regra) {
    if (exists(".print_topico", mode = "function")) .print_topico("REGRA EMP\u00cdRICA (68 \u2014 95 \u2014 99,7%)")
    cat(sprintf("  \u03bc \u00b1 1\u03c3  \u2192  [%s ; %s]  \u2192  68,27%% dos dados\n", .f2(regra68[1]),  .f2(regra68[2])))
    cat(sprintf("  \u03bc \u00b1 2\u03c3  \u2192  [%s ; %s]  \u2192  95,45%% dos dados\n", .f2(regra95[1]),  .f2(regra95[2])))
    cat(sprintf("  \u03bc \u00b1 3\u03c3  \u2192  [%s ; %s]  \u2192  99,73%% dos dados\n\n", .f2(regra997[1]), .f2(regra997[2])))
  }

  if (exists(".print_topico", mode = "function")) .print_topico("QUANTIS")
  pcts <- c("P01","P05","P10","P25","P50","P75","P90","P95","P99")
  wq <- 8
  cat("  ", paste(.pad(pcts, wq, "center"), collapse = ""), "\n", sep = "")
  cat("  ", strrep("\u2500", wq * length(pcts)), "\n", sep = "")
  cat("  ", paste(.pad(.f2(quantis), wq, "center"), collapse = ""), "\n\n", sep = "")

  if (exists(".print_rodape", mode = "function")) .print_rodape()

  # ── Funcao de amostra por tipo ────────────────────────────────────────────
  .gera <- function(nn) {
    set.seed(semente + nn)
    switch(tipo,
      normal      = rnorm(nn, mean = media, sd = dp),
      binomial    = rbinom(nn, size = n, prob = p),
      poisson     = rpois(nn, lambda = lambda),
      exponencial = rexp(nn, rate = taxa)
    )
  }
  .curva_teo <- function(gg, nn) {
    switch(tipo,
      normal = gg +
        ggplot2::stat_function(fun = dnorm, args = list(mean = media, sd = dp),
                               color = "#D90429", linewidth = 1),
      exponencial = gg +
        ggplot2::stat_function(fun = dexp, args = list(rate = taxa),
                               color = "#D90429", linewidth = 1),
      binomial = {
        am <- .gera(nn)
        xs <- min(am):max(am)
        df_teo <- data.frame(x = xs, y = dbinom(xs, size = n, prob = p))
        gg + ggplot2::geom_segment(data = df_teo,
               ggplot2::aes(x = x, xend = x, y = 0, yend = y),
               color = "#D90429", linewidth = 1) +
             ggplot2::geom_point(data = df_teo, ggplot2::aes(x = x, y = y),
               color = "#D90429", size = 2)
      },
      poisson = {
        am <- .gera(nn)
        xs <- min(am):max(am)
        df_teo <- data.frame(x = xs, y = dpois(xs, lambda = lambda))
        gg + ggplot2::geom_segment(data = df_teo,
               ggplot2::aes(x = x, xend = x, y = 0, yend = y),
               color = "#D90429", linewidth = 1) +
             ggplot2::geom_point(data = df_teo, ggplot2::aes(x = x, y = y),
               color = "#D90429", size = 2)
      }
    )
  }
  
  .tema <- function(gg) {
    if (exists("tema_estatR", mode = "function")) {
      gg + tema_estatR(estilo = 2)
    } else if (requireNamespace("estatR", quietly = TRUE)) {
      fn <- get("tema_estatR", envir = asNamespace("estatR"))
      gg + fn(estilo = 2)
    } else {
      gg + ggplot2::theme_minimal()
    }
  }

  # ── Modo Simulacao: painel 2x3 ────────────────────────────────────────────
  if (simulacao && grafico) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      message("[Aviso] Instale o pacote 'patchwork' para o painel de simulacao.")
      return(invisible(NULL))
    }

    ns <- c(10, 50, 100, 250, 500, 1000)
    plots <- lapply(ns, function(nn) {
      am <- .gera(nn)
      df <- data.frame(X = am)

      eh_continua <- tipo %in% c("normal", "exponencial")

      if (eh_continua) {
        gg <- ggplot2::ggplot(df, ggplot2::aes(x = X)) +
          ggplot2::geom_histogram(
            ggplot2::aes(y = ggplot2::after_stat(density)),
            bins    = max(8, min(30, round(sqrt(nn)))),
            fill    = "#AAAAAA", color = "#555555", alpha = 0.7
          )
      } else {
        gg <- ggplot2::ggplot(df, ggplot2::aes(x = X)) +
          ggplot2::geom_bar(
            ggplot2::aes(y = ggplot2::after_stat(prop)),
            fill = "#AAAAAA", color = "#555555", alpha = 0.7
          )
      }

      gg <- .curva_teo(gg, nn)
      gg <- gg + ggplot2::labs(
        title = sprintf("n = %d", nn),
        x = "x",
        y = if (eh_continua) "Densidade" else "Proporção"
      )
      .tema(gg)
    })

    painel <- (plots[[1]] | plots[[2]] | plots[[3]]) /
              (plots[[4]] | plots[[5]] | plots[[6]]) +
      patchwork::plot_annotation(
        title = sprintf("Simulação \u2014 %s", tools::toTitleCase(tipo)),
        subtitle = "Convergência da distribuição amostral para a distribuição teórica (curva vermelha)"
      )

    print(painel)
    return(invisible(NULL))
  }

  # ── Modo Grafico Normal (68-95-99.7 ou curva padrao) ──────────────────────
  if (grafico && !simulacao) {

    if (tipo == "normal") {
      # Grafico estilo da imagem: 3 regioes sombreadas em tons ambar
      x_seq  <- seq(media - 4*dp, media + 4*dp, length.out = 800)
      df_c   <- data.frame(x = x_seq, y = dnorm(x_seq, mean = media, sd = dp))

      .poly_region <- function(a, b) {
        sub <- df_c[df_c$x >= a & df_c$x <= b, ]
        rbind(data.frame(x = sub$x[1], y = 0), sub,
              data.frame(x = sub$x[nrow(sub)], y = 0))
      }

      df_3s <- .poly_region(media - 3*dp, media + 3*dp)
      df_2s <- .poly_region(media - 2*dp, media + 2*dp)
      df_1s <- .poly_region(media - dp,   media + dp)

      breaks_x <- c(media - 3*dp, media - 2*dp, media - dp,
                    media,
                    media + dp,   media + 2*dp,  media + 3*dp)
      labs_x <- c(
        expression(mu - 3*sigma), expression(mu - 2*sigma), expression(mu - sigma),
        expression(mu),
        expression(mu + sigma),  expression(mu + 2*sigma),  expression(mu + 3*sigma)
      )

      p <- ggplot2::ggplot(df_c, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_polygon(data = df_3s, ggplot2::aes(x = x, y = y), fill = "#F5DEB3", alpha = 1) +
        ggplot2::geom_polygon(data = df_2s, ggplot2::aes(x = x, y = y), fill = "#D4A843", alpha = 1) +
        ggplot2::geom_polygon(data = df_1s, ggplot2::aes(x = x, y = y), fill = "#B8860B", alpha = 1) +
        ggplot2::geom_line(color = "#8B6914", linewidth = 1.4) +
        ggplot2::scale_x_continuous(breaks = breaks_x, labels = labs_x) +
        ggplot2::annotate("text", x = media - 0.5*dp, y = max(df_c$y) * 0.35, label = "34,13%", size = 5.5, color = "white", fontface = "bold") +
        ggplot2::annotate("text", x = media + 0.5*dp, y = max(df_c$y) * 0.35, label = "34,13%", size = 5.5, color = "white", fontface = "bold") +
        ggplot2::annotate("text", x = media - 1.5*dp, y = max(df_c$y) * 0.12, label = "13,59%", size = 4.5, color = "#4A3000", fontface = "bold") +
        ggplot2::annotate("text", x = media + 1.5*dp, y = max(df_c$y) * 0.12, label = "13,59%", size = 4.5, color = "#4A3000", fontface = "bold") +
        ggplot2::annotate("text", x = media - 2.5*dp, y = max(df_c$y) * 0.03, label = "2,14%", size = 4, color = "#4A3000") +
        ggplot2::annotate("text", x = media + 2.5*dp, y = max(df_c$y) * 0.03, label = "2,14%", size = 4, color = "#4A3000") +
        ggplot2::labs(title = sprintf("Distribuição Normal  \u2014  \u03bc = %g, \u03c3 = %g", media, dp),
                      subtitle = "Regra Empírica: 68% \u2014 95% \u2014 99,7%", x = "x", y = "f(x)")
      p <- .tema(p)

    } else if (tipo == "exponencial") {
      x_max  <- qexp(0.999, rate = taxa)
      x_seq  <- seq(0, x_max * 1.05, length.out = 500)
      df_c   <- data.frame(x = x_seq, y = dexp(x_seq, rate = taxa))
      med_x  <- 1 / taxa

      p <- ggplot2::ggplot(df_c, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_area(fill = "#AACCFF", alpha = 0.5) +
        ggplot2::geom_line(color = "#1A5DB5", linewidth = 1.2) +
        ggplot2::geom_vline(xintercept = med_x, color = "#D90429",
                            linetype = "dashed", linewidth = 0.9) +
        ggplot2::annotate("text", x = med_x + x_max * 0.03,
                          y = dexp(0, taxa) * 0.7,
                          label = sprintf("\u03bc = %g", med_x),
                          color = "#D90429", size = 6, hjust = 0, fontface = "bold") +
        ggplot2::labs(
          title    = sprintf("Distribuição Exponencial  \u2014  \u03bb = %g", taxa),
          subtitle = sprintf("Média (1/\u03bb) = %g", med_x),
          x = "x", y = "f(x)"
        )
      p <- .tema(p)

    } else if (tipo == "binomial") {
      xs  <- 0:n
      df_b <- data.frame(x = xs, y = dbinom(xs, n, p))
      med_b <- n * p

      p <- ggplot2::ggplot(df_b, ggplot2::aes(x = factor(x), y = y)) +
        ggplot2::geom_col(
          fill = ifelse(df_b$x == round(med_b), "#1A5DB5", "#AACCFF"),
          color = "#555555", alpha = 0.85
        ) +
        ggplot2::geom_point(ggplot2::aes(x = as.numeric(factor(x))),
                            color = "#D90429", size = 2.5) +
        ggplot2::geom_line(ggplot2::aes(x = as.numeric(factor(x))),
                           color = "#D90429", linewidth = 1, group = 1) +
        ggplot2::labs(
          title    = sprintf("Distribuição Binomial  \u2014  n = %d, p = %g", n, p),
          subtitle = sprintf("Média = %g    Variância = %g", med_b, n*p*(1-p)),
          x = "Número de Sucessos (k)", y = "P(X = k)"
        )
      p <- .tema(p)

    } else if (tipo == "poisson") {
      x_max_p <- qpois(0.999, lambda)
      xs_p    <- 0:x_max_p
      df_p    <- data.frame(x = xs_p, y = dpois(xs_p, lambda))

      p <- ggplot2::ggplot(df_p, ggplot2::aes(x = factor(x), y = y)) +
        ggplot2::geom_col(
          fill = ifelse(df_p$x == round(lambda), "#1A5DB5", "#AACCFF"),
          color = "#555555", alpha = 0.85
        ) +
        ggplot2::geom_point(ggplot2::aes(x = as.numeric(factor(x))),
                            color = "#D90429", size = 2.5) +
        ggplot2::geom_line(ggplot2::aes(x = as.numeric(factor(x))),
                           color = "#D90429", linewidth = 1, group = 1) +
        ggplot2::labs(
          title    = sprintf("Distribuição Poisson  \u2014  \u03bb = %g", lambda),
          subtitle = sprintf("Média = \u03bb = %g    Variância = \u03bb = %g", lambda, lambda),
          x = "Número de Ocorrências (k)", y = "P(X = k)"
        )
      p <- .tema(p)
    }

    suppressMessages(print(p))
  }

  invisible(list(media = teo_media, variancia = teo_var, dp = teo_dp,
                 quantis = stats::setNames(quantis,
                   c("P01","P05","P10","P25","P50","P75","P90","P95","P99"))))
}
