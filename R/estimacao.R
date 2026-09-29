# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: ESTIMAÇÃO
# Funções de Estimação por Máxima Verossimilhança (EMV) e Método dos Momentos
# ─────────────────────────────────────────────────────────────────────────────

# ── Helpers internos ──────────────────────────────────────────────────────────
.est_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
  paste0(s, strrep(" ", sp))
}

.est_sep <- function(w) paste0("  ", strrep("\u2500", w))
.est_fmt <- function(x, d = 4) formatC(x, format = "f", digits = d, decimal.mark = ",")


# ─────────────────────────────────────────────────────────────────────────────
# Funções internas de log-verossimilhança
# ─────────────────────────────────────────────────────────────────────────────
.ll_normal <- function(params, x) {
  mu <- params[1]; sigma <- params[2]
  if (sigma <= 0) return(-Inf)
  sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

.ll_poisson <- function(lambda, x) {
  if (lambda <= 0) return(-Inf)
  sum(dpois(x, lambda = lambda, log = TRUE))
}

.ll_exponencial <- function(taxa, x) {
  if (taxa <= 0) return(-Inf)
  sum(dexp(x, rate = taxa, log = TRUE))
}

.ll_gamma <- function(params, x) {
  shape <- params[1]; rate <- params[2]
  if (shape <= 0 || rate <= 0) return(-Inf)
  sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Estimação por Máxima Verossimilhança (EMV)
#' @description Estima os parâmetros de uma distribuição teórica via Máxima Verossimilhança.
#'   Plota a curva da função de log-verossimilhança com o EMV destacado.
#' @param x Vetor numérico com os dados observados.
#' @param distribuicao Distribuição teórica: \code{"normal"}, \code{"poisson"},
#'   \code{"exponencial"} ou \code{"gamma"}.
#' @param decimais Casas decimais para exibição (padrão 4).
#' @param grafico Lógico. Se TRUE, plota a curva de log-verossimilhança.
#' @return Retorna invisivelmente uma lista com os parâmetros estimados.
#' @export
estimar_verossimilhanca <- function(x, distribuicao = c("normal", "poisson", "exponencial", "gamma"),
                                     decimais = 4, grafico = TRUE) {
  distribuicao <- match.arg(distribuicao)

  x_expr   <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x")) var_nome <- "x"

  x_c <- x[!is.na(x)]
  n   <- length(x_c)
  if (n < 5) stop("S\u00e3o necess\u00e1rias pelo menos 5 observa\u00e7\u00f5es.")

  w_sep <- 67

  .print_titulo("ESTIMAÇÃO POR M\u00c1XIMA VEROSSIMILHAN\u00c7A (EMV) %s")
  cat(sprintf("  Vari\u00e1vel: %s   |   N: %d   |   Distribui\u00e7\u00e3o: %s\n\n",
              var_nome, n, tools::toTitleCase(distribuicao)))

  # ── Estimação e anotações específicas por distribuição ──────────────────────
  params_emv <- list()
  theta_nome  <- character(0)
  theta_valor <- numeric(0)
  ll_max <- NULL
  xlab_plot <- NULL

  if (distribuicao == "normal") {
    mu_hat    <- mean(x_c)
    sigma_hat <- sqrt(sum((x_c - mu_hat)^2) / n)   # EMV usa /n (não /n-1)
    params_emv <- list(mu = mu_hat, sigma = sigma_hat)
    theta_nome  <- c("\u03bc (m\u00e9dia)", "\u03c3 (desvio padr\u00e3o)", "\u03c3\u00b2 (vari\u00e2ncia)")
    theta_valor <- c(mu_hat, sigma_hat, sigma_hat^2)
    ll_max      <- .ll_normal(c(mu_hat, sigma_hat), x_c)
    xlab_plot   <- "\u03bc"

  } else if (distribuicao == "poisson") {
    if (any(x_c != floor(x_c)) || any(x_c < 0)) stop("Poisson requer dados inteiros n\u00e3o-negativos.")
    lambda_hat <- mean(x_c)
    params_emv <- list(lambda = lambda_hat)
    theta_nome  <- c("\u03bb (m\u00e9dia / taxa)")
    theta_valor <- c(lambda_hat)
    ll_max      <- .ll_poisson(lambda_hat, x_c)
    xlab_plot   <- "\u03bb"

  } else if (distribuicao == "exponencial") {
    if (any(x_c <= 0)) stop("Exponencial requer dados estritamente positivos.")
    taxa_hat <- 1 / mean(x_c)
    params_emv <- list(taxa = taxa_hat)
    theta_nome  <- c("\u03bb (taxa)", "\u03bc = 1/\u03bb (m\u00e9dia)")
    theta_valor <- c(taxa_hat, 1 / taxa_hat)
    ll_max      <- .ll_exponencial(taxa_hat, x_c)
    xlab_plot   <- "\u03bb"

  } else if (distribuicao == "gamma") {
    if (any(x_c <= 0)) stop("Gamma requer dados estritamente positivos.")
    m1 <- mean(x_c)
    m2 <- mean(x_c^2)
    shape0 <- m1^2 / (m2 - m1^2)
    rate0  <- m1   / (m2 - m1^2)
    opt <- optim(c(shape = shape0, rate = rate0), function(p) -.ll_gamma(p, x_c),
                 method = "L-BFGS-B", lower = c(1e-6, 1e-6))
    shape_hat <- opt$par[1]
    rate_hat  <- opt$par[2]
    params_emv <- list(shape = shape_hat, rate = rate_hat)
    theta_nome  <- c("\u03b1 (forma / shape)", "\u03b2 (taxa / rate)", "\u03bc = \u03b1/\u03b2 (m\u00e9dia)")
    theta_valor <- c(shape_hat, rate_hat, shape_hat / rate_hat)
    ll_max      <- .ll_gamma(c(shape_hat, rate_hat), x_c)
    xlab_plot   <- "\u03b1"
  }

  # ── Tabela de resultados ─────────────────────────────────────────────────────
  cat("  \u25b6 PAR\u00c2METROS ESTIMADOS (EMV)\n")
  wc <- c(par = 28, val = 20)
  hdr <- paste0(.est_pad("Par\u00e2metro", wc["par"], "left"), .est_pad("Estimativa", wc["val"], "center"))
  cat(.est_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.est_sep(nchar(hdr)), "\n")
  for (i in seq_along(theta_nome)) {
    cat("  ",
        .est_pad(theta_nome[i],                wc["par"], "left"),
        .est_pad(.est_fmt(theta_valor[i], decimais), wc["val"], "center"),
        "\n", sep = "")
  }
  cat(.est_sep(nchar(hdr)), "\n")
  cat(sprintf("  Log-verossimilhan\u00e7a m\u00e1xima: %s\n\n", .est_fmt(ll_max, 3)))

  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  cat(sprintf("  Os par\u00e2metros da distribui\u00e7\u00e3o %s foram estimados\n", tools::toTitleCase(distribuicao)))
  cat(sprintf("  maximizando a fun\u00e7\u00e3o de verossimilhan\u00e7a com base nas\n"))
  cat(sprintf("  %d observa\u00e7\u00f5es fornecidas. Os valores apresentados s\u00e3o\n", n))
  cat(sprintf("  os que tornam os dados observados 'mais prov\u00e1veis' de ocorrer.\n"))
  .print_rodape()

  # ── Gráfico: curva de log-verossimilhança ─────────────────────────────────────
  if (grafico && exists("meu_tema")) {
    tryCatch({

      if (distribuicao == "normal") {
        sigma_fix <- params_emv$sigma
        mu_hat    <- params_emv$mu
        theta_seq <- seq(mu_hat - 4 * sigma_fix / sqrt(n), mu_hat + 4 * sigma_fix / sqrt(n), length.out = 300)
        ll_seq    <- sapply(theta_seq, function(m) .ll_normal(c(m, sigma_fix), x_c))
        emv_val   <- mu_hat
        emv_ll    <- ll_max
        titulo    <- sprintf("Curva de Log-Verossimilhan\u00e7a — %s", tools::toTitleCase(distribuicao))
        sub_txt   <- sprintf("EMV: \u03bc\u0302 = %s  (com \u03c3 = %s fixado)", .est_fmt(mu_hat, 3), .est_fmt(sigma_fix, 3))

      } else if (distribuicao == "poisson") {
        lam_hat   <- params_emv$lambda
        lam_min   <- max(0.001, lam_hat * 0.3)
        lam_max   <- lam_hat * 2
        theta_seq <- seq(lam_min, lam_max, length.out = 300)
        ll_seq    <- sapply(theta_seq, .ll_poisson, x = x_c)
        emv_val   <- lam_hat
        emv_ll    <- ll_max
        titulo    <- "Curva de Log-Verossimilhan\u00e7a — Poisson"
        sub_txt   <- sprintf("EMV: \u03bb\u0302 = %s", .est_fmt(lam_hat, 3))

      } else if (distribuicao == "exponencial") {
        taxa_hat  <- params_emv$taxa
        taxa_min  <- taxa_hat * 0.3
        taxa_max  <- taxa_hat * 2.5
        theta_seq <- seq(taxa_min, taxa_max, length.out = 300)
        ll_seq    <- sapply(theta_seq, .ll_exponencial, x = x_c)
        emv_val   <- taxa_hat
        emv_ll    <- ll_max
        titulo    <- "Curva de Log-Verossimilhan\u00e7a — Exponencial"
        sub_txt   <- sprintf("EMV: \u03bb\u0302 = %s  (\u03bc\u0302 = %s)", .est_fmt(taxa_hat, 3), .est_fmt(1/taxa_hat, 3))

      } else { # gamma — perfil em relação ao shape
        shape_hat <- params_emv$shape
        rate_hat  <- params_emv$rate
        sh_min    <- max(0.01, shape_hat * 0.4)
        sh_max    <- shape_hat * 2
        theta_seq <- seq(sh_min, sh_max, length.out = 300)
        ll_seq    <- sapply(theta_seq, function(s) .ll_gamma(c(s, rate_hat), x_c))
        emv_val   <- shape_hat
        emv_ll    <- ll_max
        titulo    <- "Curva de Log-Verossimilhan\u00e7a — Gamma"
        sub_txt   <- sprintf("EMV: \u03b1\u0302 = %s  (com \u03b2 = %s fixado)", .est_fmt(shape_hat, 3), .est_fmt(rate_hat, 3))
        xlab_plot <- "\u03b1"
      }

      df_curv <- data.frame(theta = theta_seq, ll = ll_seq)
      df_curv <- df_curv[is.finite(df_curv$ll), ]

      old_w <- getOption("warn"); options(warn = -1)
      p_plot <- ggplot2::ggplot(df_curv, ggplot2::aes(x = theta, y = ll)) +
        ggplot2::geom_line(color = "#333333", linewidth = 1.1) +
        ggplot2::geom_vline(xintercept = emv_val, linetype = "dashed",
                            color = "#1B4F72", linewidth = 0.8) +
        ggplot2::geom_point(data = data.frame(x = emv_val, y = emv_ll),
                            ggplot2::aes(x = x, y = y),
                            color = "#D90429", size = 5, shape = 18) +
        ggplot2::annotate("text", x = emv_val, y = emv_ll,
                          label = sprintf(" \u03b8\u0302 = %s", .est_fmt(emv_val, 3)),
                          hjust = -0.1, vjust = 0.5, size = 4.5,
                          fontface = "bold", color = "#D90429") +
        ggplot2::labs(title = titulo, subtitle = sub_txt,
                      x = xlab_plot, y = "log L(\u03b8)", caption = "estatR") +
        meu_tema()

      suppressMessages(print(p_plot))
      options(warn = old_w)
    }, error = function(e) {
      message("[Aviso] Falha ao plotar a curva de verossimilhan\u00e7a: ", e$message)
    })
  }

  invisible(params_emv)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Estimação pelo Método dos Momentos
#' @description Estima os parâmetros de uma distribuição igualando os momentos
#'   amostrais aos momentos populacionais teóricos.
#' @param x Vetor numérico com os dados observados.
#' @param distribuicao Distribuição teórica: \code{"normal"}, \code{"poisson"},
#'   \code{"exponencial"} ou \code{"gamma"}.
#' @param decimais Casas decimais para exibição (padrão 4).
#' @param grafico Lógico. Se TRUE, plota histograma com a curva teórica ajustada.
#' @return Retorna invisivelmente uma lista com os parâmetros estimados.
#' @export
estimar_momentos <- function(x, distribuicao = c("normal", "poisson", "exponencial", "gamma"),
                              decimais = 4, grafico = TRUE) {
  distribuicao <- match.arg(distribuicao)

  x_expr   <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x")) var_nome <- "x"

  x_c <- x[!is.na(x)]
  n   <- length(x_c)
  if (n < 5) stop("S\u00e3o necess\u00e1rias pelo menos 5 observa\u00e7\u00f5es.")

  m1 <- mean(x_c)
  m2 <- mean(x_c^2)
  v  <- var(x_c)
  dp <- sd(x_c)

  w_sep <- 67

  .print_titulo("ESTIMAÇÃO PELO M\u00c9TODO DOS MOMENTOS %s")
  cat(sprintf("  Vari\u00e1vel: %s   |   N: %d   |   Distribui\u00e7\u00e3o: %s\n\n",
              var_nome, n, tools::toTitleCase(distribuicao)))

  # ── Tabela de Momentos Amostrais ─────────────────────────────────────────────
  cat("  \u25b6 MOMENTOS AMOSTRAIS UTILIZADOS\n")
  wc <- c(mom = 28, val = 20)
  hdr_m <- paste0(.est_pad("Momento", wc["mom"], "left"), .est_pad("Valor", wc["val"], "center"))
  cat(.est_sep(nchar(hdr_m)), "\n")
  cat("  ", hdr_m, "\n", sep = "")
  cat(.est_sep(nchar(hdr_m)), "\n")
  cat("  ", .est_pad("1\u00ba Momento Amostral (m\u00e9dia)", wc["mom"], "left"), .est_pad(.est_fmt(m1, decimais), wc["val"], "center"), "\n", sep = "")
  cat("  ", .est_pad("2\u00ba Momento Amostral (variância)", wc["mom"], "left"), .est_pad(.est_fmt(v, decimais), wc["val"], "center"), "\n", sep = "")
  cat(.est_sep(nchar(hdr_m)), "\n\n")

  # ── Estimação por distribuição ────────────────────────────────────────────────
  params_mm <- list()
  theta_nome  <- character(0)
  theta_valor <- numeric(0)
  formula_txt <- character(0)

  if (distribuicao == "normal") {
    mu_hat    <- m1
    sigma_hat <- dp      # MM para normal: μ = m1, σ² = s² (variância amostral)
    params_mm  <- list(mu = mu_hat, sigma = sigma_hat)
    theta_nome  <- c("\u03bc\u0302 (m\u00e9dia)", "\u03c3\u0302 (desvio padr\u00e3o)", "\u03c3\u0302\u00b2 (vari\u00e2ncia)")
    theta_valor <- c(mu_hat, sigma_hat, sigma_hat^2)
    formula_txt <- c("\u03bc\u0302 = m\u2081 = \u0078\u0305", "\u03c3\u0302 = \u221a(m\u2082 - m\u2081\u00b2) = s")

  } else if (distribuicao == "poisson") {
    if (any(x_c != floor(x_c)) || any(x_c < 0)) stop("Poisson requer dados inteiros n\u00e3o-negativos.")
    lambda_hat <- m1    # Para Poisson: E[X] = λ = m1
    params_mm  <- list(lambda = lambda_hat)
    theta_nome  <- c("\u03bb\u0302 (m\u00e9dia / taxa)")
    theta_valor <- c(lambda_hat)
    formula_txt <- c("\u03bb\u0302 = m\u2081 = \u0078\u0305")

  } else if (distribuicao == "exponencial") {
    if (any(x_c <= 0)) stop("Exponencial requer dados estritamente positivos.")
    taxa_hat <- 1 / m1   # Para Exp: E[X] = 1/λ, logo λ = 1/m1
    params_mm  <- list(taxa = taxa_hat)
    theta_nome  <- c("\u03bb\u0302 (taxa)", "\u03bc\u0302 = 1/\u03bb\u0302 (m\u00e9dia)")
    theta_valor <- c(taxa_hat, m1)
    formula_txt <- c("\u03bb\u0302 = 1/m\u2081 = 1/\u0078\u0305")

  } else { # gamma
    if (any(x_c <= 0)) stop("Gamma requer dados estritamente positivos.")
    # Para Gamma: E[X] = α/β, Var[X] = α/β²  →  α = m1²/s², β = m1/s²
    shape_hat <- m1^2 / v
    rate_hat  <- m1   / v
    params_mm  <- list(shape = shape_hat, rate = rate_hat)
    theta_nome  <- c("\u03b1\u0302 (forma / shape)", "\u03b2\u0302 (taxa / rate)", "\u03bc\u0302 = \u03b1\u0302/\u03b2\u0302 (m\u00e9dia)")
    theta_valor <- c(shape_hat, rate_hat, shape_hat / rate_hat)
    formula_txt <- c("\u03b1\u0302 = m\u2081\u00b2 / s\u00b2", "\u03b2\u0302 = m\u2081 / s\u00b2")
  }

  # ── Fórmulas utilizadas ───────────────────────────────────────────────────────
  cat("  \u25b6 EQUA\u00c7\u00f5ES DOS MOMENTOS\n")
  for (f in formula_txt) cat(sprintf("  %s\n", f))
  cat("\n")

  # ── Tabela de parâmetros estimados ────────────────────────────────────────────
  cat("  \u25b6 PAR\u00c2METROS ESTIMADOS (Momentos)\n")
  hdr_p <- paste0(.est_pad("Par\u00e2metro", wc["mom"], "left"), .est_pad("Estimativa", wc["val"], "center"))
  cat(.est_sep(nchar(hdr_p)), "\n")
  cat("  ", hdr_p, "\n", sep = "")
  cat(.est_sep(nchar(hdr_p)), "\n")
  for (i in seq_along(theta_nome)) {
    cat("  ",
        .est_pad(theta_nome[i],                     wc["mom"], "left"),
        .est_pad(.est_fmt(theta_valor[i], decimais), wc["val"], "center"),
        "\n", sep = "")
  }
  cat(.est_sep(nchar(hdr_p)), "\n\n")

  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  cat(sprintf("  Os par\u00e2metros foram obtidos igualando os momentos\n"))
  cat(sprintf("  amostrais (m\u00e9dia e vari\u00e2ncia) aos momentos te\u00f3ricos da\n"))
  cat(sprintf("  distribui\u00e7\u00e3o %s e resolvendo o sistema de equa\u00e7\u00f5es.\n", tools::toTitleCase(distribuicao)))
  .print_rodape()

  # ── Gráfico: histograma + curva teórica ajustada ─────────────────────────────
  if (grafico && exists("meu_tema")) {
    tryCatch({
      df_obs <- data.frame(x = x_c)

      x_min <- min(x_c); x_max <- max(x_c)
      x_seq <- seq(x_min * 0.9, x_max * 1.1, length.out = 300)

      if (distribuicao == "normal") {
        df_teo <- data.frame(x = x_seq, y = dnorm(x_seq, mean = params_mm$mu, sd = params_mm$sigma))
        sub_txt <- sprintf("\u03bc\u0302 = %s,  \u03c3\u0302 = %s", .est_fmt(params_mm$mu, 3), .est_fmt(params_mm$sigma, 3))
      } else if (distribuicao == "poisson") {
        x_int   <- seq(floor(x_min), ceiling(x_max))
        df_teo  <- data.frame(x = x_int, y = dpois(x_int, lambda = params_mm$lambda))
        sub_txt <- sprintf("\u03bb\u0302 = %s", .est_fmt(params_mm$lambda, 3))
      } else if (distribuicao == "exponencial") {
        x_seq_e <- seq(0, x_max * 1.2, length.out = 300)
        df_teo  <- data.frame(x = x_seq_e, y = dexp(x_seq_e, rate = params_mm$taxa))
        sub_txt <- sprintf("\u03bb\u0302 = %s  (\u03bc\u0302 = %s)", .est_fmt(params_mm$taxa, 3), .est_fmt(1/params_mm$taxa, 3))
      } else { # gamma
        x_seq_g <- seq(0.001, x_max * 1.2, length.out = 300)
        df_teo  <- data.frame(x = x_seq_g, y = dgamma(x_seq_g, shape = params_mm$shape, rate = params_mm$rate))
        sub_txt <- sprintf("\u03b1\u0302 = %s,  \u03b2\u0302 = %s", .est_fmt(params_mm$shape, 3), .est_fmt(params_mm$rate, 3))
      }

      old_w <- getOption("warn"); options(warn = -1)

      if (distribuicao == "poisson") {
        # Para discretas: barras do histograma + pontos/linhas da distribuição teórica
        df_freq <- as.data.frame(table(x_c))
        df_freq$x_c <- as.numeric(as.character(df_freq$x_c))
        df_freq$prop <- df_freq$Freq / sum(df_freq$Freq)

        p_plot <- ggplot2::ggplot() +
          ggplot2::geom_col(data = df_freq, ggplot2::aes(x = x_c, y = prop),
                            fill = "#888888", color = "white", alpha = 0.7, width = 0.6) +
          ggplot2::geom_point(data = df_teo, ggplot2::aes(x = x, y = y),
                              color = "#D90429", size = 3) +
          ggplot2::geom_line(data = df_teo, ggplot2::aes(x = x, y = y),
                             color = "#D90429", linewidth = 0.9) +
          ggplot2::labs(
            title    = sprintf("Ajuste Poisson — M\u00e9todo dos Momentos"),
            subtitle = sub_txt,
            x = var_nome, y = "Propor\u00e7\u00e3o", caption = "estatR"
          ) +
          meu_tema()
      } else {
        p_plot <- ggplot2::ggplot(df_obs, ggplot2::aes(x = x)) +
          ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                  bins = min(30, ceiling(sqrt(n))),
                                  fill = "#888888", color = "white", alpha = 0.7) +
          ggplot2::geom_line(data = df_teo, ggplot2::aes(x = x, y = y),
                             color = "#D90429", linewidth = 1.2) +
          ggplot2::labs(
            title    = sprintf("Ajuste %s — M\u00e9todo dos Momentos", tools::toTitleCase(distribuicao)),
            subtitle = sub_txt,
            x = var_nome, y = "Densidade", caption = "estatR"
          ) +
          meu_tema()
      }

      suppressMessages(print(p_plot))
      options(warn = old_w)
    }, error = function(e) {
      message("[Aviso] Falha ao plotar o ajuste: ", e$message)
    })
  }

  invisible(params_mm)
}
