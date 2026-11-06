# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: INFERÊNCIA BAYESIANA
# Testes de Hipótese Bayesianos, Fatores de Bayes (JZS) e Intervalos de Credibilidade (HPD)
# ─────────────────────────────────────────────────────────────────────────────

# ── HELPERS INTERNOS BAYESIANOS ───────────────────────────────────────────────

.bayes_fmt_num <- function(x, d = 4) {
  if (is.null(x) || length(x) == 0 || is.na(x) || is.nan(x)) return("-")
  formatC(as.numeric(x), format = "f", digits = d, decimal.mark = ",")
}

.bayes_fmt_pct <- function(x, d = 2) {
  if (is.null(x) || length(x) == 0 || is.na(x) || is.nan(x)) return("-")
  paste0(formatC(as.numeric(x), format = "f", digits = d, decimal.mark = ","), "%")
}

.bayes_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp / 2)), s, strrep(" ", ceiling(sp / 2))))
  paste0(s, strrep(" ", sp))
}

.bayes_sep <- function(w) paste0("  ", strrep("\u2500", w))

# Cálculo exato do Intervalo HPD (Highest Posterior Density)
.bayes_hpd_amostras <- function(amostras, confianca = 0.95) {
  x <- sort(as.numeric(stats::na.omit(amostras)))
  n <- length(x)
  if (n < 10) return(as.numeric(stats::quantile(x, probs = c((1 - confianca) / 2, 1 - (1 - confianca) / 2))))
  k <- max(1, floor(n * confianca))
  diffs <- x[(k + 1):n] - x[1:(n - k)]
  idx_min <- which.min(diffs)
  c(inf = x[idx_min], sup = x[idx_min + k])
}

# HPD analítico para distribuição Beta(a, b)
.bayes_hpd_beta <- function(a, b, confianca = 0.95) {
  if (a <= 1 || b <= 1) {
    # Em distribuições em U ou J, o intervalo de caudas iguais (ETI) é o padrão estável
    alpha <- 1 - confianca
    return(c(inf = stats::qbeta(alpha / 2, a, b), sup = stats::qbeta(1 - alpha / 2, a, b)))
  }
  # Otimização unidimensional para encontrar p1 que minimiza a amplitude
  fn_largura <- function(p1) {
    stats::qbeta(p1 + confianca, a, b) - stats::qbeta(p1, a, b)
  }
  opt <- stats::optimize(fn_largura, interval = c(1e-6, 1 - confianca - 1e-6))
  p_inf <- opt$minimum
  c(inf = stats::qbeta(p_inf, a, b), sup = stats::qbeta(p_inf + confianca, a, b))
}

# Classificação canônica do Fator de Bayes (Jeffreys / Lee & Wagenmakers)
.bayes_interpretar_bf <- function(bf10) {
  if (is.na(bf10) || is.infinite(bf10)) return(list(texto = "Indeterminado", nivel = "indeterminado", positivo = FALSE))
  if (bf10 >= 100) {
    list(texto = "Evidência Extrema / Decisiva a favor de H1", nivel = "decisiva", positivo = TRUE)
  } else if (bf10 >= 30) {
    list(texto = "Evidência Muito Forte a favor de H1", nivel = "muito_forte", positivo = TRUE)
  } else if (bf10 >= 10) {
    list(texto = "Evidência Forte a favor de H1", nivel = "forte", positivo = TRUE)
  } else if (bf10 >= 3) {
    list(texto = "Evidência Moderada a favor de H1", nivel = "moderada", positivo = TRUE)
  } else if (bf10 >= 1) {
    list(texto = "Evidência Fraca / Anedótica a favor de H1 (inconclusivo)", nivel = "fraca", positivo = FALSE)
  } else if (bf10 >= 1/3) {
    list(texto = "Evidência Fraca a favor de H0 (inconclusivo)", nivel = "fraca_h0", positivo = FALSE)
  } else if (bf10 >= 1/10) {
    list(texto = "Evidência Moderada a favor de H0", nivel = "moderada_h0", positivo = TRUE)
  } else if (bf10 >= 1/30) {
    list(texto = "Evidência Forte a favor de H0", nivel = "forte_h0", positivo = TRUE)
  } else {
    list(texto = "Evidência Muito Forte / Decisiva a favor de H0", nivel = "decisiva_h0", positivo = TRUE)
  }
}


# ── 1. TESTE DE PROPORÇÃO BAYESIANO (BETA-BINOMIAL) ──────────────────────────

#' @title Teste de Proporção Bayesiano (Modelo Beta-Binomial)
#' @description Executa inferência bayesiana para uma proporção binomial utilizando o modelo
#'   conjugado Beta-Binomial, fornecendo estimativas pontuais (média, moda e mediana a posteriori),
#'   Intervalo de Credibilidade de Maior Densidade a Posteriori (HPD 95%), probabilidade de hipótese
#'   e Fator de Bayes via Savage-Dickey density ratio.
#'
#' @param x Número de sucessos observados ou vetor binário de dados.
#' @param n Número total de tentativas (ignorado se \code{x} for um vetor).
#' @param p0 Valor nulo de referência para a proporção (padrão \code{0.5}).
#' @param hipotese Direção da hipótese alternativa: \code{"bilateral"} (padrão, \eqn{p \neq p_0}),
#'   \code{"maior"} (\eqn{p > p_0}) ou \code{"menor"} (\eqn{p < p_0}).
#' @param priori Distribuição a priori: \code{"auto"} (seleciona a priori não-informativa de Jeffreys),
#'   \code{"nao_informativa"} (\eqn{\text{Beta}(0.5, 0.5)}), \code{"uniforme"} (\eqn{\text{Beta}(1, 1)})
#'   ou \code{"manual"} informando \code{alfa_priori} e \code{beta_priori}.
#' @param alfa_priori Parâmetro \eqn{\alpha} da priori Beta quando \code{priori = "manual"}.
#' @param beta_priori Parâmetro \eqn{\beta} da priori Beta quando \code{priori = "manual"}.
#' @param confianca Nível de probabilidade para o intervalo de credibilidade (padrão \code{0.95}).
#' @param grafico Lógico. Se TRUE (padrão), exibe painel gráfico com Priori vs Posteriori e região HPD.
#'
#' @return Objeto invisível da classe \code{estatR_proporcao_bayesiana}.
#' @export
teste_proporcao_bayesiano <- function(x,
                                      n = NULL,
                                      p0 = 0.5,
                                      hipotese = c("bilateral", "maior", "menor"),
                                      priori = c("auto", "nao_informativa", "uniforme", "manual"),
                                      alfa_priori = NULL,
                                      beta_priori = NULL,
                                      confianca = 0.95,
                                      grafico = TRUE) {

  hipotese <- match.arg(hipotese)
  priori   <- match.arg(priori)

  # Tratamento dos dados de entrada
  if (is.null(n)) {
    if (length(x) > 1) {
      if (is.factor(x) || is.character(x)) x <- as.numeric(as.factor(x)) - 1
      x_vals <- x[!is.na(x)]
      k <- sum(x_vals == 1 | x_vals == max(x_vals))
      n_total <- length(x_vals)
    } else {
      stop("Quando 'x' for escalar, informe o tamanho amostral 'n'.", call. = FALSE)
    }
  } else {
    k <- as.integer(x[1])
    n_total <- as.integer(n)
  }

  if (k < 0 || k > n_total || n_total <= 0) {
    stop("O número de sucessos deve estar entre 0 e n.", call. = FALSE)
  }

  p_obs <- k / n_total

  # Definição dos parâmetros da priori
  if (priori %in% c("auto", "nao_informativa")) {
    a_priori <- 0.5
    b_priori <- 0.5
    nome_priori <- "Jeffreys Beta(0,5; 0,5)"
  } else if (priori == "uniforme") {
    a_priori <- 1.0
    b_priori <- 1.0
    nome_priori <- "Uniforme Beta(1; 1)"
  } else {
    if (is.null(alfa_priori) || is.null(beta_priori) || alfa_priori <= 0 || beta_priori <= 0) {
      stop("Para 'priori = manual', forneça 'alfa_priori' e 'beta_priori' positivos.", call. = FALSE)
    }
    a_priori <- as.numeric(alfa_priori)
    b_priori <- as.numeric(beta_priori)
    nome_priori <- sprintf("Customizada Beta(%s; %s)", .bayes_fmt_num(a_priori, 2), .bayes_fmt_num(b_priori, 2))
  }

  # Atualização Bayesiana (Posteriori)
  a_post <- a_priori + k
  b_post <- b_priori + (n_total - k)

  # Estimativas pontuais
  media_post <- a_post / (a_post + b_post)
  moda_post  <- if (a_post > 1 && b_post > 1) (a_post - 1) / (a_post + b_post - 2) else media_post
  mediana_post <- stats::qbeta(0.5, a_post, b_post)
  var_post   <- (a_post * b_post) / (((a_post + b_post)^2) * (a_post + b_post + 1))
  dp_post    <- sqrt(var_post)

  # Intervalo de Credibilidade HPD e ETI
  ic_hpd <- .bayes_hpd_beta(a_post, b_post, confianca = confianca)
  alpha_tail <- (1 - confianca) / 2
  ic_eti <- c(inf = stats::qbeta(alpha_tail, a_post, b_post), sup = stats::qbeta(1 - alpha_tail, a_post, b_post))

  # Teste de Hipótese e Probabilidades a Posteriori
  prob_maior <- 1 - stats::pbeta(p0, a_post, b_post)
  prob_menor <- stats::pbeta(p0, a_post, b_post)

  prob_hipotese <- switch(hipotese,
    "maior"     = prob_maior,
    "menor"     = prob_menor,
    "bilateral" = 2 * min(prob_maior, prob_menor)
  )

  # Fator de Bayes (Savage-Dickey density ratio para ponto nulo)
  dens_priori_p0 <- stats::dbeta(p0, a_priori, b_priori)
  dens_post_p0   <- stats::dbeta(p0, a_post, b_post)
  bf01 <- if (dens_post_p0 > 0) dens_post_p0 / dens_priori_p0 else 0
  bf10 <- if (bf01 > 0) 1 / bf01 else Inf

  # Se hipótese direcional, corrige o Fator de Bayes pela probabilidade direcional
  if (hipotese == "maior") {
    priori_dir <- 1 - stats::pbeta(p0, a_priori, b_priori)
    if (priori_dir > 0 && prob_maior > 0) {
      bf10 <- (prob_maior / (1 - prob_maior)) / (priori_dir / (1 - priori_dir))
      bf01 <- 1 / bf10
    }
  } else if (hipotese == "menor") {
    priori_dir <- stats::pbeta(p0, a_priori, b_priori)
    if (priori_dir > 0 && prob_menor > 0) {
      bf10 <- (prob_menor / (1 - prob_menor)) / (priori_dir / (1 - priori_dir))
      bf01 <- 1 / bf10
    }
  }

  diag_bf <- .bayes_interpretar_bf(bf10)

  # ── Console Output ──────────────────────────────────────────────────────────
  .print_titulo("TESTE DE PROPORÇÃO BAYESIANO (BETA-BINOMIAL)")
  cat(sprintf("  Sucessos (k)        : %d em %d ensaios (Taxa amostral: %s)\n",
              k, n_total, .bayes_fmt_pct(p_obs * 100, 2)))
  cat(sprintf("  Distribuição a Priori: %s\n", nome_priori))
  cat(sprintf("  Posteriori Atualizada: Beta(\u03b1 = %s; \u03b2 = %s)\n",
              .bayes_fmt_num(a_post, 2), .bayes_fmt_num(b_post, 2)))
  cat(sprintf("  Valor Nulo de Teste : p₀ = %s | Hipótese: %s\n\n",
              .bayes_fmt_num(p0, 2), hipotese))

  .print_topico(sprintf("ESTIMATIVAS PONTUAIS E INTERVALO DE CREDIBILIDADE (%.0f%%)", confianca * 100))
  w_prop <- c(est = 24, val = 16)
  cat(.bayes_sep(sum(w_prop)), "\n")
  cat("  ", .bayes_pad("Parâmetro", w_prop["est"]), .bayes_pad("Valor", w_prop["val"], "right"), "\n", sep = "")
  cat(.bayes_sep(sum(w_prop)), "\n")
  cat("  ", .bayes_pad("Média a Posteriori", w_prop["est"]), .bayes_pad(.bayes_fmt_num(media_post, 4), w_prop["val"], "right"), "\n", sep = "")
  cat("  ", .bayes_pad("Moda a Posteriori", w_prop["est"]), .bayes_pad(.bayes_fmt_num(moda_post, 4), w_prop["val"], "right"), "\n", sep = "")
  cat("  ", .bayes_pad("Mediana a Posteriori", w_prop["est"]), .bayes_pad(.bayes_fmt_num(mediana_post, 4), w_prop["val"], "right"), "\n", sep = "")
  cat("  ", .bayes_pad("Desvio Padrão Posterior", w_prop["est"]), .bayes_pad(.bayes_fmt_num(dp_post, 4), w_prop["val"], "right"), "\n", sep = "")
  hpd_str <- sprintf("[%s; %s]", .bayes_fmt_num(ic_hpd[1], 4), .bayes_fmt_num(ic_hpd[2], 4))
  cat("  ", .bayes_pad(sprintf("IC HPD %.0f%% (Maior Dens.)", confianca * 100), w_prop["est"]), .bayes_pad(hpd_str, w_prop["val"], "right"), "\n", sep = "")
  eti_str <- sprintf("[%s; %s]", .bayes_fmt_num(ic_eti[1], 4), .bayes_fmt_num(ic_eti[2], 4))
  cat("  ", .bayes_pad(sprintf("IC ETI %.0f%% (Caudas Iguais)", confianca * 100), w_prop["est"]), .bayes_pad(eti_str, w_prop["val"], "right"), "\n", sep = "")
  cat(.bayes_sep(sum(w_prop)), "\n\n")

  .print_topico("AVALIAÇÃO DE HIPÓTESE E EVIDÊNCIA BAYESIANA")
  cat(sprintf("  P(p > p₀ | dados)   : %s\n", .bayes_fmt_pct(prob_maior * 100, 2)))
  cat(sprintf("  P(p < p₀ | dados)   : %s\n", .bayes_fmt_pct(prob_menor * 100, 2)))
  cat(sprintf("  Fator de Bayes BF₁₀ : %s  (evidência a favor de H1)\n", .bayes_fmt_num(bf10, 3)))
  cat(sprintf("  Fator de Bayes BF₀₁ : %s  (evidência a favor de H0)\n\n", .bayes_fmt_num(bf01, 3)))

  # Interpretação com [✓] ou [!]
  if (bf10 >= 3) {
    cat(sprintf("  [✓] %s (BF₁₀ = %s).\n", diag_bf$texto, .bayes_fmt_num(bf10, 2)))
    cat(sprintf("      Os dados são ~%s vezes mais compatíveis com a hipótese H1 do que com o valor nulo p₀ = %s.\n\n",
                .bayes_fmt_num(bf10, 1), .bayes_fmt_num(p0, 2)))
  } else if (bf10 <= 1/3) {
    cat(sprintf("  [✓] %s (BF₀₁ = %s).\n", diag_bf$texto, .bayes_fmt_num(bf01, 2)))
    cat(sprintf("      Os dados fornecem suporte à hipótese nula p₀ = %s em detrimento da alternativa.\n\n",
                .bayes_fmt_num(p0, 2)))
  } else {
    cat(sprintf("  [!] Evidência inconclusiva / anedótica (BF₁₀ = %s; BF₀₁ = %s).\n",
                .bayes_fmt_num(bf10, 2), .bayes_fmt_num(bf01, 2)))
    cat("      A amostra atual não permite diferenciar claramente entre H0 e H1 com convicção bayesiana.\n\n")
  }
  .print_rodape()

  # ── Gráficos em ggplot2 com meu_tema(grade = "dupla") ───────────────────────
  p_combined <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    p_seq <- seq(1e-4, 1 - 1e-4, length.out = 300)
    dens_pri <- stats::dbeta(p_seq, a_priori, b_priori)
    dens_pos <- stats::dbeta(p_seq, a_post, b_post)

    df_p <- data.frame(
      p = rep(p_seq, 2),
      densidade = c(dens_pri, dens_pos),
      Tipo = factor(rep(c("Priori", "Posteriori"), each = length(p_seq)), levels = c("Priori", "Posteriori"))
    )

    # Subconjunto do intervalo HPD para sombreamento
    df_hpd <- subset(df_p, Tipo == "Posteriori" & p >= ic_hpd[1] & p <= ic_hpd[2])

    p1 <- ggplot2::ggplot(df_p, ggplot2::aes(x = p, y = densidade, color = Tipo)) +
      ggplot2::geom_ribbon(data = df_hpd, ggplot2::aes(ymin = 0, ymax = densidade),
                           fill = "#457B9D", alpha = 0.35, inherit.aes = FALSE) +
      ggplot2::geom_line(ggplot2::aes(linetype = Tipo), linewidth = 1.1) +
      ggplot2::scale_color_manual(values = c("Priori" = "gray50", "Posteriori" = "#1D3557")) +
      ggplot2::scale_linetype_manual(values = c("Priori" = "dashed", "Posteriori" = "solid")) +
      ggplot2::geom_vline(xintercept = p0, color = "#E63946", linewidth = 0.9, linetype = "dotted") +
      ggplot2::annotate("text", x = p0, y = max(dens_pos) * 0.9, label = sprintf("p0 = %s", .bayes_fmt_num(p0, 2)),
                        color = "#E63946", vjust = -0.5, hjust = -0.1, size = 3.5) +
      ggplot2::labs(
        title    = "Atualização Bayesiana: Distribuição a Priori vs Posteriori",
        subtitle = sprintf("Área sombreada: IC HPD %.0f%% [%s; %s] | Média posterior: %s",
                           confianca * 100, .bayes_fmt_num(ic_hpd[1], 2), .bayes_fmt_num(ic_hpd[2], 2), .bayes_fmt_num(media_post, 3)),
        x        = "Proporção (p)",
        y        = "Densidade",
        color    = "", linetype = ""
      ) +
      meu_tema(grade = "dupla") +
      ggplot2::theme(legend.position = "top")

    # Painel 2: CDF Acumulada a Posteriori
    df_cdf <- data.frame(p = p_seq, cdf = stats::pbeta(p_seq, a_post, b_post))
    p2 <- ggplot2::ggplot(df_cdf, ggplot2::aes(x = p, y = cdf)) +
      ggplot2::geom_line(color = "#1D3557", linewidth = 1.1) +
      ggplot2::geom_vline(xintercept = p0, color = "#E63946", linetype = "dotted", linewidth = 0.9) +
      ggplot2::geom_hline(yintercept = prob_menor, color = "#457B9D", linetype = "dashed", linewidth = 0.8) +
      ggplot2::annotate("point", x = p0, y = prob_menor, color = "#E63946", size = 3) +
      ggplot2::labs(
        title    = "Função de Distribuição Acumulada a Posteriori",
        subtitle = sprintf("P(p \u2264 p0 | dados) = %s | P(p > p0 | dados) = %s",
                           .bayes_fmt_pct(prob_menor * 100, 1), .bayes_fmt_pct(prob_maior * 100, 1)),
        x        = "Proporção (p)",
        y        = "Probabilidade Acumulada"
      ) +
      meu_tema(grade = "dupla")

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      p_combined <- list(painel_densidade = p1, painel_acumulada = p2)
    } else {
      print(p1)
      p_combined <- p1
    }
  }

  res <- list(
    media_posterior = media_post,
    moda_posterior  = moda_post,
    mediana_posterior = mediana_post,
    dp_posterior    = dp_post,
    ic_hpd          = ic_hpd,
    ic_eti          = ic_eti,
    prob_maior      = prob_maior,
    prob_menor      = prob_menor,
    bf10            = bf10,
    bf01            = bf01,
    interpretacao   = diag_bf$texto,
    priori          = c(a = a_priori, b = b_priori),
    posteriori      = c(a = a_post, b = b_post),
    dados           = c(k = k, n = n_total, taxa_obs = p_obs),
    confianca       = confianca,
    grafico         = p_combined
  )
  class(res) <- c("estatR_proporcao_bayesiana", "list")
  invisible(res)
}

#' @export
ic_proporcao_bayesiano <- teste_proporcao_bayesiano


# ── 2. TESTE T BAYESIANO (FATOR DE BAYES JZS & POSTERIORI) ───────────────────

#' @title Teste t Bayesiano (Fator de Bayes JZS & Distribuição a Posteriori)
#' @description Realiza teste de hipótese bayesiano para médias (uma amostra, duas amostras
#'   independentes ou pareadas) calculando o Fator de Bayes exato Jeffreys-Zellner-Siow (JZS)
#'   via integração analítica em R puro e a distribuição a posteriori do tamanho de efeito (\eqn{\delta} de Cohen).
#'
#' @param x Vetor numérico contendo as observações da amostra.
#' @param y Vetor numérico opcional para teste de duas amostras ou pareado.
#' @param mu Média de referência sob a hipótese nula \eqn{H_0} (padrão \code{0}).
#' @param pareado Lógico. Se TRUE, realiza teste t pareado entre \code{x} e \code{y}.
#' @param hipotese Direção da hipótese: \code{"bilateral"} (\eqn{\mu \neq \mu_0}),
#'   \code{"maior"} (\eqn{\mu > \mu_0}) ou \code{"menor"} (\eqn{\mu < \mu_0}).
#' @param r_escala Escala da priori de Cauchy para o tamanho de efeito (\eqn{r}):
#'   \code{"auto"} (padrão \eqn{\sqrt{2}/2 \approx 0.707}), \code{"larga"} (\code{1.0})
#'   ou valor numérico positivo.
#' @param confianca Nível de probabilidade para o intervalo HPD do tamanho de efeito (padrão \code{0.95}).
#' @param grafico Lógico. Se TRUE (padrão), exibe curva da priori vs posteriori do tamanho de efeito \eqn{\delta}.
#'
#' @return Objeto invisível da classe \code{estatR_teste_t_bayesiano}.
#' @export
teste_t_bayesiano <- function(x,
                              y = NULL,
                              mu = 0,
                              pareado = FALSE,
                              hipotese = c("bilateral", "maior", "menor"),
                              r_escala = "auto",
                              confianca = 0.95,
                              grafico = TRUE) {

  hipotese <- match.arg(hipotese)
  x_clean <- stats::na.omit(as.numeric(x))

  # Definição do r da priori de Cauchy
  r_val <- if (identical(r_escala, "auto")) sqrt(2) / 2 else if (identical(r_escala, "larga")) 1.0 else as.numeric(r_escala)
  if (r_val <= 0) r_val <- sqrt(2) / 2

  tipo_teste <- "Uma Amostra"
  # Estruturação amostral
  if (!is.null(y)) {
    y_clean <- stats::na.omit(as.numeric(y))
    if (isTRUE(pareado)) {
      tipo_teste <- "Pareado"
      if (length(x) != length(y)) stop("Para teste pareado, x e y devem ter o mesmo comprimento.", call. = FALSE)
      difs <- stats::na.omit(x - y)
      n_eff <- length(difs)
      gl <- n_eff - 1
      m_diff <- mean(difs)
      s_diff <- stats::sd(difs)
      t_stat <- (m_diff - mu) / (s_diff / sqrt(n_eff))
      d_obs  <- (m_diff - mu) / s_diff
      amostra_info <- list(n = n_eff, media = m_diff, sd = s_diff)
    } else {
      tipo_teste <- "Duas Amostras Independentes"
      n1 <- length(x_clean)
      n2 <- length(y_clean)
      gl <- n1 + n2 - 2
      n_eff <- (n1 * n2) / (n1 + n2)
      m1 <- mean(x_clean); m2 <- mean(y_clean)
      s1 <- stats::sd(x_clean); s2 <- stats::sd(y_clean)
      sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / gl)
      t_stat <- (m1 - m2 - mu) / (sp * sqrt(1/n1 + 1/n2))
      d_obs  <- (m1 - m2 - mu) / sp
      amostra_info <- list(n1 = n1, n2 = n2, media1 = m1, media2 = m2, sp = sp)
    }
  } else {
    tipo_teste <- "Uma Amostra"
    n_eff <- length(x_clean)
    gl <- n_eff - 1
    m_x <- mean(x_clean)
    s_x <- stats::sd(x_clean)
    t_stat <- (m_x - mu) / (s_x / sqrt(n_eff))
    d_obs  <- (m_x - mu) / s_x
    amostra_info <- list(n = n_eff, media = m_x, sd = s_x)
  }

  if (gl < 1) stop("Graus de liberdade insuficientes para o teste t bayesiano.", call. = FALSE)

  # ── Cálculo exato do Fator de Bayes JZS via integração ───────────────────────
  # Densidade marginal sob H1: integral sobre g de (1 + N*g)^(-1/2) * [1 + t^2 / (nu * (1 + N*g))]^(-(nu+1)/2) * pi(g) dg
  # onde g ~ Inv-Gamma(1/2, r^2 / 2)
  fn_integrando <- function(g) {
    # Evita instabilidade numérica para g muito pequeno ou grande
    term1 <- (1 + n_eff * g)^(-0.5)
    term2 <- (1 + (t_stat^2) / (gl * (1 + n_eff * g)))^(-(gl + 1) / 2)
    pi_g  <- ((r_val^2 / 2)^0.5) / (gamma(0.5)) * (g^(-1.5)) * exp(-(r_val^2) / (2 * g))
    term1 * term2 * pi_g
  }

  int_res <- stats::integrate(fn_integrando, lower = 1e-8, upper = Inf, rel.tol = 1e-7)
  m1 <- int_res$value
  m0 <- (1 + (t_stat^2) / gl)^(-(gl + 1) / 2)

  bf10_bilateral <- m1 / m0

  # Correção para hipótese direcional
  # Simulação direta de delta posterior para obter probabilidades direcionais e HPD
  set.seed(42)
  n_mc <- 20000
  # Amostragem de g da Inv-Gamma
  g_draws <- 1 / stats::rgamma(n_mc, shape = 0.5, rate = r_val^2 / 2)
  # Amostragem condicional de delta | g, t
  # E(delta | g, t) = (N*g / (1 + N*g)) * d_obs
  # Var(delta | g, t) = g / (1 + N*g)
  var_delta_g <- g_draws / (1 + n_eff * g_draws)
  mean_delta_g <- (n_eff * g_draws / (1 + n_eff * g_draws)) * (t_stat / sqrt(n_eff))
  
  # Peso de importância para cada draw
  pesos_w <- fn_integrando(g_draws) / (stats::dgamma(1 / g_draws, shape = 0.5, rate = r_val^2 / 2) * (g_draws^(-2)))
  pesos_w <- pesos_w / sum(pesos_w)
  
  # Amostras do tamanho de efeito delta
  idx_samp <- sample(seq_len(n_mc), size = n_mc, replace = TRUE, prob = pesos_w)
  delta_post <- stats::rnorm(n_mc, mean = mean_delta_g[idx_samp], sd = sqrt(var_delta_g[idx_samp]))

  prob_delta_maior <- mean(delta_post > 0)
  prob_delta_menor <- mean(delta_post < 0)

  bf10 <- if (hipotese == "maior") {
    bf10_bilateral * 2 * prob_delta_maior
  } else if (hipotese == "menor") {
    bf10_bilateral * 2 * prob_delta_menor
  } else {
    bf10_bilateral
  }

  bf01 <- if (bf10 > 0) 1 / bf10 else Inf
  diag_bf <- .bayes_interpretar_bf(bf10)

  # Estimativas a posteriori do tamanho de efeito delta (Cohen's d)
  delta_media   <- mean(delta_post)
  delta_mediana <- stats::median(delta_post)
  delta_sd      <- stats::sd(delta_post)
  ic_hpd_delta  <- .bayes_hpd_amostras(delta_post, confianca = confianca)

  # ── Console Output ──────────────────────────────────────────────────────────
  .print_titulo("TESTE T BAYESIANO (FATOR DE BAYES JZS & POSTERIORI)")
  cat(sprintf("  Modalidade do Teste : %s\n", tipo_teste))
  cat(sprintf("  Estatística t       : %s (gl = %d | Cohen's d observado = %s)\n",
              .bayes_fmt_num(t_stat, 3), gl, .bayes_fmt_num(d_obs, 3)))
  cat(sprintf("  Priori do Efeito (δ): Cauchy(0, r = %s)\n", .bayes_fmt_num(r_val, 3)))
  cat(sprintf("  Hipótese Alternativa: %s (referência: %s)\n\n", hipotese, .bayes_fmt_num(mu, 2)))

  .print_topico("FATOR DE BAYES E EVIDÊNCIA ESTATÍSTICA")
  cat(sprintf("  Fator de Bayes BF₁₀ : %s  (evidência a favor de H1 vs H0)\n", .bayes_fmt_num(bf10, 3)))
  cat(sprintf("  Fator de Bayes BF₀₁ : %s  (evidência a favor de H0 vs H1)\n", .bayes_fmt_num(bf01, 3)))
  cat(sprintf("  Classificação       : %s\n\n", diag_bf$texto))

  if (bf10 >= 3) {
    cat(sprintf("  [✓] Evidência estatística consistente a favor de H1 (BF₁₀ = %s).\n", .bayes_fmt_num(bf10, 2)))
    cat(sprintf("      Os dados são ~%s vezes mais prováveis sob a presença de efeito do que sob H0.\n\n",
                .bayes_fmt_num(bf10, 1)))
  } else if (bf10 <= 1/3) {
    cat(sprintf("  [✓] Evidência estatística a favor da hipótese nula H0 (BF₀₁ = %s).\n", .bayes_fmt_num(bf01, 2)))
    cat("      Os dados fornecem suporte à ausência de diferença estatisticamente relevante.\n\n")
  } else {
    cat(sprintf("  [!] Evidência inconclusiva / anedótica (BF₁₀ = %s).\n", .bayes_fmt_num(bf10, 2)))
    cat("      Não há evidência bayesiana suficiente para confirmar ou refutar a hipótese.\n\n")
  }

  .print_topico(sprintf("DISTRIBUIÇÃO A POSTERIORI DO TAMANHO DE EFEITO (δ DE COHEN) - %.0f%%", confianca * 100))
  w_tt <- c(est = 24, val = 16)
  cat(.bayes_sep(sum(w_tt)), "\n")
  cat("  ", .bayes_pad("Parâmetro", w_tt["est"]), .bayes_pad("Estimativa", w_tt["val"], "right"), "\n", sep = "")
  cat(.bayes_sep(sum(w_tt)), "\n")
  cat("  ", .bayes_pad("Média a Posteriori (δ)", w_tt["est"]), .bayes_pad(.bayes_fmt_num(delta_media, 4), w_tt["val"], "right"), "\n", sep = "")
  cat("  ", .bayes_pad("Mediana a Posteriori", w_tt["est"]), .bayes_pad(.bayes_fmt_num(delta_mediana, 4), w_tt["val"], "right"), "\n", sep = "")
  cat("  ", .bayes_pad("Desvio Padrão Posterior", w_tt["est"]), .bayes_pad(.bayes_fmt_num(delta_sd, 4), w_tt["val"], "right"), "\n", sep = "")
  hpd_t_str <- sprintf("[%s; %s]", .bayes_fmt_num(ic_hpd_delta[1], 4), .bayes_fmt_num(ic_hpd_delta[2], 4))
  cat("  ", .bayes_pad(sprintf("IC HPD %.0f%% (Credibilidade)", confianca * 100), w_tt["est"]), .bayes_pad(hpd_t_str, w_tt["val"], "right"), "\n", sep = "")
  cat(.bayes_sep(sum(w_tt)), "\n\n")
  .print_rodape()

  # ── Gráficos em ggplot2 com meu_tema(grade = "dupla") ───────────────────────
  p_combined <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    d_seq <- seq(min(delta_post), max(delta_post), length.out = 300)
    # Priori Cauchy(0, r)
    dens_pri_c <- stats::dcauchy(d_seq, location = 0, scale = r_val)
    dens_post_dens <- stats::density(delta_post, from = min(d_seq), to = max(d_seq), n = 300)$y

    df_plot_d <- data.frame(
      delta = rep(d_seq, 2),
      densidade = c(dens_pri_c, dens_post_dens),
      Tipo = factor(rep(c("Priori Cauchy", "Posteriori"), each = length(d_seq)), levels = c("Priori Cauchy", "Posteriori"))
    )

    df_hpd_d <- subset(df_plot_d, Tipo == "Posteriori" & delta >= ic_hpd_delta[1] & delta <= ic_hpd_delta[2])

    p1 <- ggplot2::ggplot(df_plot_d, ggplot2::aes(x = delta, y = densidade, color = Tipo)) +
      ggplot2::geom_ribbon(data = df_hpd_d, ggplot2::aes(ymin = 0, ymax = densidade),
                           fill = "#457B9D", alpha = 0.35, inherit.aes = FALSE) +
      ggplot2::geom_line(ggplot2::aes(linetype = Tipo), linewidth = 1.1) +
      ggplot2::scale_color_manual(values = c("Priori Cauchy" = "gray50", "Posteriori" = "#1D3557")) +
      ggplot2::scale_linetype_manual(values = c("Priori Cauchy" = "dashed", "Posteriori" = "solid")) +
      ggplot2::geom_vline(xintercept = 0, color = "#E63946", linetype = "dotted", linewidth = 0.9) +
      ggplot2::labs(
        title    = "Distribuição do Tamanho de Efeito (δ de Cohen): Priori vs Posteriori",
        subtitle = sprintf("Área sombreada: IC HPD %.0f%% [%s; %s] | BF₁₀ = %s (%s)",
                           confianca * 100, .bayes_fmt_num(ic_hpd_delta[1], 2), .bayes_fmt_num(ic_hpd_delta[2], 2),
                           .bayes_fmt_num(bf10, 2), diag_bf$nivel),
        x        = "Tamanho de Efeito (δ de Cohen)",
        y        = "Densidade",
        color    = "", linetype = ""
      ) +
      meu_tema(grade = "dupla") +
      ggplot2::theme(legend.position = "top")

    # Painel 2: Barra de Força da Evidência (Visualizador de BF10)
    # Mostra barra comparativa H0 vs H1 em escala de probabilidade posterior
    prob_h1_post <- bf10 / (1 + bf10)
    prob_h0_post <- 1 - prob_h1_post
    df_bar_bf <- data.frame(
      Hipótese = factor(c("H0 (Sem Efeito)", "H1 (Com Efeito)"), levels = c("H0 (Sem Efeito)", "H1 (Com Efeito)")),
      Probabilidade = c(prob_h0_post, prob_h1_post)
    )

    p2 <- ggplot2::ggplot(df_bar_bf, ggplot2::aes(x = Hipótese, y = Probabilidade, fill = Hipótese)) +
      ggplot2::geom_col(width = 0.55, color = "white", alpha = 0.85) +
      ggplot2::scale_fill_manual(values = c("H0 (Sem Efeito)" = "#E63946", "H1 (Com Efeito)" = "#1D3557")) +
      ggplot2::scale_y_continuous(labels = function(y) .bayes_fmt_pct(y * 100, 0), limits = c(0, 1)) +
      ggplot2::geom_text(ggplot2::aes(label = .bayes_fmt_pct(Probabilidade * 100, 1)), vjust = -0.4, size = 3.8, fontface = "bold") +
      ggplot2::labs(
        title    = "Probabilidade Posterior das Hipóteses",
        subtitle = sprintf("Priori balanceada 50%%/50%% | BF₁₀ = %s", .bayes_fmt_num(bf10, 2)),
        x        = "",
        y        = "Probabilidade a Posteriori"
      ) +
      meu_tema(grade = "dupla") +
      ggplot2::theme(legend.position = "none")

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      p_combined <- list(painel_densidade = p1, painel_evidencia = p2)
    } else {
      print(p1)
      p_combined <- p1
    }
  }

  res <- list(
    bf10             = bf10,
    bf01             = bf01,
    interpretacao    = diag_bf$texto,
    nivel_evidencia  = diag_bf$nivel,
    t_stat           = t_stat,
    gl               = gl,
    d_observado      = d_obs,
    delta_posterior  = list(media = delta_media, mediana = delta_mediana, sd = delta_sd, hpd = ic_hpd_delta),
    hipotese         = hipotese,
    tipo_teste       = tipo_teste,
    confianca        = confianca,
    grafico          = p_combined
  )
  class(res) <- c("estatR_teste_t_bayesiano", "list")
  invisible(res)
}

#' @export
media_bayesiana <- teste_t_bayesiano


# ── 3. COMPARAÇÃO A/B BAYESIANA (DUAS PROPORÇÕES) ────────────────────────────

#' @title Teste A/B Bayesiano (Comparação de Duas Proporções)
#' @description Compara duas taxas/proporções binomiais independentes (Grupo A vs Grupo B)
#'   utilizando modelo Beta-Binomial com cálculo exato da probabilidade de superioridade,
#'   elevação esperada (uplift), Risco Relativo e Odds Ratio a posteriori com IC HPD.
#'
#' @param k1 Sucessos observados no Grupo 1 (ou A).
#' @param n1 Total de ensaios no Grupo 1 (ou A).
#' @param k2 Sucessos observados no Grupo 2 (ou B).
#' @param n2 Total de ensaios no Grupo 2 (ou B).
#' @param nomes Vetor de 2 caracteres nomeando os grupos (padrão \code{c("Grupo A", "Grupo B")}).
#' @param priori Priori Beta: \code{"auto"} (Jeffreys \eqn{\text{Beta}(0.5, 0.5)}) ou \code{"uniforme"} (\eqn{\text{Beta}(1, 1)}).
#' @param confianca Nível de probabilidade para o intervalo HPD (padrão \code{0.95}).
#' @param n_sim Número de simulações de Monte Carlo para a diferença a posteriori (padrão 25.000).
#' @param grafico Lógico. Se TRUE (padrão), exibe distribuições posteriores sobrepostas e densidade do efeito delta.
#'
#' @return Objeto invisível da classe \code{estatR_ab_bayesiano}.
#' @export
teste_duas_proporcoes_bayesiano <- function(k1, n1,
                                            k2, n2,
                                            nomes = c("Grupo A", "Grupo B"),
                                            priori = c("auto", "nao_informativa", "uniforme"),
                                            confianca = 0.95,
                                            n_sim = 25000,
                                            grafico = TRUE) {

  priori <- match.arg(priori)
  a_priori <- if (priori == "uniforme") 1.0 else 0.5
  b_priori <- if (priori == "uniforme") 1.0 else 0.5

  # Taxas observadas
  p1_obs <- k1 / n1
  p2_obs <- k2 / n2

  # Parâmetros das posteriors Beta
  a1_post <- a_priori + k1; b1_post <- b_priori + n1 - k1
  a2_post <- a_priori + k2; b2_post <- b_priori + n2 - k2

  # Simulação de alta precisão para a diferença e razões
  set.seed(123)
  draws_p1 <- stats::rbeta(n_sim, a1_post, b1_post)
  draws_p2 <- stats::rbeta(n_sim, a2_post, b2_post)

  diff_draws   <- draws_p2 - draws_p1
  uplift_draws <- (diff_draws / draws_p1) * 100
  rr_draws     <- draws_p2 / draws_p1
  or_draws     <- (draws_p2 / (1 - draws_p2)) / (draws_p1 / (1 - draws_p1))

  # Probabilidade de superioridade de B sobre A
  prob_vitoria_b <- mean(draws_p2 > draws_p1)

  # Estimativas da diferença
  media_diff   <- mean(diff_draws)
  mediana_diff <- stats::median(diff_draws)
  sd_diff      <- stats::sd(diff_draws)
  hpd_diff     <- .bayes_hpd_amostras(diff_draws, confianca = confianca)

  # Estimativas do uplift (%)
  media_uplift <- mean(uplift_draws)
  hpd_uplift   <- .bayes_hpd_amostras(uplift_draws, confianca = confianca)

  # Odds ratio médio
  or_medio <- stats::median(or_draws)
  hpd_or   <- .bayes_hpd_amostras(or_draws, confianca = confianca)

  nome_a <- nomes[1]
  nome_b <- nomes[2]

  # ── Console Output ──────────────────────────────────────────────────────────
  .print_titulo("TESTE A/B BAYESIANO (COMPARAÇÃO DE DUAS PROPORÇÕES)")
  cat(sprintf("  %s : %d sucessos em %d ensaios (Taxa: %s)\n",
              nome_a, k1, n1, .bayes_fmt_pct(p1_obs * 100, 2)))
  cat(sprintf("  %s : %d sucessos em %d ensaios (Taxa: %s)\n",
              nome_b, k2, n2, .bayes_fmt_pct(p2_obs * 100, 2)))
  cat(sprintf("  Priori Utilizada    : Beta(%s, %s)\n\n", .bayes_fmt_num(a_priori, 1), .bayes_fmt_num(b_priori, 1)))

  .print_topico(sprintf("ESTIMATIVAS DA DIFERENÇA E UPLIFT (%s - %s)", nome_b, nome_a))
  w_ab <- c(param = 26, val = 18)
  cat(.bayes_sep(sum(w_ab)), "\n")
  cat("  ", .bayes_pad("Métrica", w_ab["param"]), .bayes_pad("Estimativa Posterior", w_ab["val"], "right"), "\n", sep = "")
  cat(.bayes_sep(sum(w_ab)), "\n")
  cat("  ", .bayes_pad("Diferença Absoluta Média", w_ab["param"]), .bayes_pad(.bayes_fmt_pct(media_diff * 100, 3), w_ab["val"], "right"), "\n", sep = "")
  hpd_d_str <- sprintf("[%s; %s]", .bayes_fmt_pct(hpd_diff[1] * 100, 2), .bayes_fmt_pct(hpd_diff[2] * 100, 2))
  cat("  ", .bayes_pad(sprintf("IC HPD %.0f%% (Diferença)", confianca * 100), w_ab["param"]), .bayes_pad(hpd_d_str, w_ab["val"], "right"), "\n", sep = "")
  cat("  ", .bayes_pad("Elevação Média (Uplift %)", w_ab["param"]), .bayes_pad(.bayes_fmt_pct(media_uplift, 2), w_ab["val"], "right"), "\n", sep = "")
  hpd_up_str <- sprintf("[%s; %s]", .bayes_fmt_pct(hpd_uplift[1], 1), .bayes_fmt_pct(hpd_uplift[2], 1))
  cat("  ", .bayes_pad(sprintf("IC HPD %.0f%% (Uplift)", confianca * 100), w_ab["param"]), .bayes_pad(hpd_up_str, w_ab["val"], "right"), "\n", sep = "")
  cat("  ", .bayes_pad("Odds Ratio (OR)", w_ab["param"]), .bayes_pad(.bayes_fmt_num(or_medio, 3), w_ab["val"], "right"), "\n", sep = "")
  hpd_or_str <- sprintf("[%s; %s]", .bayes_fmt_num(hpd_or[1], 3), .bayes_fmt_num(hpd_or[2], 3))
  cat("  ", .bayes_pad(sprintf("IC HPD %.0f%% (Odds Ratio)", confianca * 100), w_ab["param"]), .bayes_pad(hpd_or_str, w_ab["val"], "right"), "\n", sep = "")
  cat(.bayes_sep(sum(w_ab)), "\n\n")

  .print_topico("DECISÃO E PROBABILIDADE DE VITÓRIA")
  cat(sprintf("  P(%s > %s | dados) : %s\n", nome_b, nome_a, .bayes_fmt_pct(prob_vitoria_b * 100, 2)))
  cat(sprintf("  P(%s > %s | dados) : %s\n\n", nome_a, nome_b, .bayes_fmt_pct((1 - prob_vitoria_b) * 100, 2)))

  if (prob_vitoria_b >= 0.95) {
    cat(sprintf("  [✓] Forte evidência de superioridade de %s (probabilidade = %s).\n",
                nome_b, .bayes_fmt_pct(prob_vitoria_b * 100, 1)))
    cat(sprintf("      A taxa de %s supera a de %s com elevação média esperada de +%s.\n\n",
                nome_b, nome_a, .bayes_fmt_pct(media_uplift, 1)))
  } else if (prob_vitoria_b <= 0.05) {
    cat(sprintf("  [✓] Forte evidência de superioridade de %s (probabilidade = %s).\n",
                nome_a, .bayes_fmt_pct((1 - prob_vitoria_b) * 100, 1)))
    cat(sprintf("      A taxa de %s supera a de %s.\n\n", nome_a, nome_b))
  } else {
    cat(sprintf("  [!] Diferença inconclusiva entre os grupos (probabilidade de %s vencer = %s).\n",
                nome_b, .bayes_fmt_pct(prob_vitoria_b * 100, 1)))
    cat("      Recomenda-se coletar mais observações antes de adotar uma decisão definitiva.\n\n")
  }
  .print_rodape()

  # ── Gráficos em ggplot2 com meu_tema(grade = "dupla") ───────────────────────
  p_combined <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    # Painel 1: Posteriores sobrepostas
    df_post_ab <- data.frame(
      Taxa  = c(draws_p1, draws_p2),
      Grupo = factor(rep(c(nome_a, nome_b), each = n_sim), levels = c(nome_a, nome_b))
    )

    p1 <- ggplot2::ggplot(df_post_ab, ggplot2::aes(x = Taxa, fill = Grupo, color = Grupo)) +
      ggplot2::geom_density(alpha = 0.35, linewidth = 1.0) +
      ggplot2::scale_fill_manual(values = setNames(c("#E63946", "#1D3557"), c(nome_a, nome_b))) +
      ggplot2::scale_color_manual(values = setNames(c("#E63946", "#1D3557"), c(nome_a, nome_b))) +
      ggplot2::scale_x_continuous(labels = function(x) .bayes_fmt_pct(x * 100, 1)) +
      ggplot2::labs(
        title    = "Distribuições a Posteriori das Taxas de Conversão",
        subtitle = sprintf("%s: %s | %s: %s", nome_a, .bayes_fmt_pct(p1_obs * 100, 1), nome_b, .bayes_fmt_pct(p2_obs * 100, 1)),
        x        = "Taxa (p)",
        y        = "Densidade",
        color    = "", fill = ""
      ) +
      meu_tema(grade = "dupla") +
      ggplot2::theme(legend.position = "top")

    # Painel 2: Distribuição da Diferença B - A
    df_diff <- data.frame(Diferenca = diff_draws)
    p2 <- ggplot2::ggplot(df_diff, ggplot2::aes(x = Diferenca)) +
      ggplot2::geom_density(fill = "#457B9D", color = "#1D3557", alpha = 0.45, linewidth = 1.0) +
      ggplot2::geom_vline(xintercept = 0, color = "#E63946", linetype = "dashed", linewidth = 0.9) +
      ggplot2::scale_x_continuous(labels = function(x) .bayes_fmt_pct(x * 100, 1)) +
      ggplot2::labs(
        title    = sprintf("Distribuição da Diferença (%s - %s)", nome_b, nome_a),
        subtitle = sprintf("Probabilidade de Vitória (%s > %s) = %s", nome_b, nome_a, .bayes_fmt_pct(prob_vitoria_b * 100, 1)),
        x        = "Diferença Absoluta",
        y        = "Densidade"
      ) +
      meu_tema(grade = "dupla")

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      p_combined <- list(painel_taxas = p1, painel_diferenca = p2)
    } else {
      print(p1)
      p_combined <- p1
    }
  }

  res <- list(
    p1_obs          = p1_obs,
    p2_obs          = p2_obs,
    prob_vitoria_b  = prob_vitoria_b,
    media_diff      = media_diff,
    ic_hpd_diff     = hpd_diff,
    media_uplift    = media_uplift,
    ic_hpd_uplift   = hpd_uplift,
    odds_ratio      = or_medio,
    ic_hpd_or       = hpd_or,
    nomes           = nomes,
    confianca       = confianca,
    grafico         = p_combined
  )
  class(res) <- c("estatR_ab_bayesiano", "list")
  invisible(res)
}

#' @export
teste_ab_bayesiano <- teste_duas_proporcoes_bayesiano


# ── 4. INTERVALO DE MAIOR DENSIDADE A POSTERIORI (HPD GERAL) ─────────────────

#' @title Intervalo de Maior Densidade a Posteriori (HPD)
#' @description Calcula o Intervalo de Maior Densidade a Posteriori (Highest Posterior Density - HPD)
#'   para um vetor de amostras MCMC ou distribuição posterior empírica.
#'
#' @param amostras Vetor numérico com amostras da distribuição a posteriori.
#' @param confianca Nível de probabilidade (padrão \code{0.95}).
#'
#' @return Vetor numérico com dois elementos nomeados: \code{c(inf = ..., sup = ...)}.
#' @export
ic_hpd <- function(amostras, confianca = 0.95) {
  .bayes_hpd_amostras(amostras, confianca = confianca)
}

# ── 5. TESTE DE CORRELAÇÃO BAYESIANO (PEARSON E SPEARMAN) ────────────────────

#' @title Teste de Correlação Bayesiano
#' @description Avalia a evidência de associação linear (Pearson) ou monotônica (Spearman)
#'   entre duas variáveis sob o paradigma bayesiano, fornecendo Fator de Bayes \eqn{BF_{10}}
#'   analítico exato (Ly, Verhagen & Wagenmakers, 2016) via densidade conjunta, estimativas pontuais
#'   (média e mediana a posteriori), Intervalo de Maior Densidade a Posteriori (HPD 95%) e probabilidade
#'   direcional. Compatível com fórmulas, vetores e pipe nativo (\code{|>}).
#'
#' @aliases correlacao_bayesiana
#' @param x Vetor numérico com a primeira variável, fórmula do tipo \code{y ~ x} ou data frame quando usado com pipe.
#' @param y Vetor numérico com a segunda variável (quando \code{x} for vetor) ou fórmula quando usado com pipe.
#' @param dados Data frame opcional contendo as variáveis especificadas na fórmula.
#' @param metodo Método de correlação: \code{"pearson"} (padrão) ou \code{"spearman"}.
#' @param hipotese Direção da hipótese alternativa: \code{"bilateral"} (padrão, \eqn{\rho \neq 0}),
#'   \code{"positiva"} (\eqn{\rho > 0}) ou \code{"negativa"} (\eqn{\rho < 0}).
#' @param priori Distribuição a priori para \eqn{\rho}: \code{"auto"} (uniforme padrão em \eqn{[-1, 1]}),
#'   \code{"uniforme"} ou \code{"manual"} informando o parâmetro \code{kappa}.
#' @param kappa Parâmetro de precisão/escala da priori simétrica quando \code{priori = "manual"}.
#' @param confianca Nível de probabilidade do intervalo de credibilidade (padrão \code{0.95}).
#' @param decimais Número de casas decimais na exibição (padrão \code{3}).
#' @param grafico Lógico. Se TRUE (padrão), plota a distribuição a posteriori de \eqn{\rho} com a região HPD sombreada.
#'
#' @return Objeto invisível da classe \code{"estatR_correlacao_bayesiana"} contendo \code{bf10}, \code{bf01},
#'   \code{r_amostral}, \code{media_post}, \code{mediana_post}, \code{dp_post}, \code{hpd} e probabilidades direcionais.
#' @export
teste_correlacao_bayesiano <- function(x, y = NULL, dados = NULL,
                                       metodo = c("pearson", "spearman"),
                                       hipotese = c("bilateral", "positiva", "negativa"),
                                       priori = c("auto", "uniforme", "manual"),
                                       kappa = 1,
                                       confianca = 0.95,
                                       decimais = 3,
                                       grafico = TRUE) {

  metodo   <- match.arg(metodo)
  hipotese <- match.arg(hipotese)
  priori   <- match.arg(priori)

  # Resolução transparente de fórmula ou vetores
  res_biv <- .resolver_dados_bivariados(x, y, dados, env = parent.frame())
  vx      <- res_biv$y
  vy      <- res_biv$grupo
  nome_x  <- res_biv$nome_y
  nome_y  <- res_biv$nome_grupo

  ok <- !is.na(vx) & !is.na(vy)
  vx <- vx[ok]; vy <- vy[ok]
  n  <- length(vx)
  if (n < 4) stop("São necessárias pelo menos 4 observações válidas para o teste de correlação.", call. = FALSE)

  if (metodo == "spearman") {
    vx <- rank(vx)
    vy <- rank(vy)
  }

  r_amostral <- stats::cor(vx, vy)

  fn_integrando_z <- function(z, rh, r, n) 1 / (cosh(z) - rh * r)^(n - 1)

  fn_prior <- function(rh) {
    if (abs(rh) >= 1) return(0)
    if (priori == "manual") {
      0.5 * (1 - rh^2)^((kappa - 1) / 2)
    } else {
      0.5
    }
  }

  grid_rho <- seq(-0.998, 0.998, length.out = 500)
  post_raw <- sapply(grid_rho, function(rh) {
    int_z <- tryCatch(stats::integrate(fn_integrando_z, lower = 0, upper = Inf, rh = rh, r = r_amostral, n = n)$value,
                      error = function(e) 0)
    (1 - rh^2)^((n - 1) / 2) * int_z * fn_prior(rh)
  })

  d_rho <- grid_rho[2] - grid_rho[1]
  const_norm <- sum(post_raw) * d_rho
  if (const_norm <= 0) const_norm <- 1e-12
  post_dens <- post_raw / const_norm

  int_z0 <- stats::integrate(fn_integrando_z, lower = 0, upper = Inf, rh = 0, r = r_amostral, n = n)$value
  dens_h0_raw <- int_z0 * fn_prior(0)
  dens_post_0 <- dens_h0_raw / const_norm
  dens_prior_0 <- fn_prior(0)

  bf01 <- dens_post_0 / dens_prior_0
  bf10 <- if (bf01 > 0) 1 / bf01 else Inf

  prob_pos <- sum(post_dens[grid_rho > 0]) * d_rho
  prob_neg <- 1 - prob_pos

  if (hipotese == "positiva") {
    bf10 <- bf10 * (prob_pos / 0.5)
    bf01 <- 1 / bf10
  } else if (hipotese == "negativa") {
    bf10 <- bf10 * (prob_neg / 0.5)
    bf01 <- 1 / bf10
  }

  media_post <- sum(grid_rho * post_dens) * d_rho
  var_post   <- sum(((grid_rho - media_post)^2) * post_dens) * d_rho
  dp_post    <- sqrt(max(0, var_post))
  moda_post  <- grid_rho[which.max(post_dens)]

  cdf_post <- cumsum(post_dens) * d_rho
  cdf_post <- cdf_post / max(cdf_post)
  mediana_post <- grid_rho[which.min(abs(cdf_post - 0.5))]

  alturas <- sort(unique(post_dens), decreasing = TRUE)
  hpd_inf <- -0.99; hpd_sup <- 0.99
  for (h in alturas) {
    selec <- post_dens >= h
    area <- sum(post_dens[selec]) * d_rho
    if (area >= confianca) {
      rhos_selec <- grid_rho[selec]
      hpd_inf <- min(rhos_selec)
      hpd_sup <- max(rhos_selec)
      break
    }
  }

  interp_bf <- .bayes_interpretar_bf(bf10)
  tipo_nome <- if (metodo == "spearman") "Spearman (Não-Paramétrica / Postos)" else "Pearson (Linear)"

  cat("\n── TESTE DE CORRELAÇÃO BAYESIANO ──────────────────────────────────────────\n\n")
  cat(sprintf("  Variáveis          : %s e %s\n", nome_x, nome_y))
  cat(sprintf("  Tipo de Correlação : %s\n", tipo_nome))
  cat(sprintf("  Tamanho Amostral   : N = %d pares válidos\n", n))
  cat(sprintf("  Coeficiente Amostral: r = %s\n\n", formatC(r_amostral, format = "f", digits = decimais, decimal.mark = ",")))

  cat(sprintf("  H0: \u03c1 = 0  (Ausência de correlação linear)\n"))
  cat(sprintf("  H1: \u03c1 %s  (Existe correlação %s)\n\n",
              if (hipotese == "bilateral") "\u2260 0" else if (hipotese == "positiva") "> 0" else "< 0",
              if (hipotese == "bilateral") "não-nula" else if (hipotese == "positiva") "positiva" else "negativa"))

  cat(.bayes_sep(75), "\n")
  wc <- c(bf10 = 16, bf01 = 16, evid = 40)
  hdr <- paste0(.bayes_pad("Fator Bayes BF\u2081\u2080", wc["bf10"], "center"),
                .bayes_pad("Fator Bayes BF\u2080\u2081", wc["bf01"], "center"),
                .bayes_pad("Força da Evidência", wc["evid"], "left"))
  cat("  ", hdr, "\n", sep = "")
  cat(.bayes_sep(75), "\n")

  txt_bf10 <- if (is.infinite(bf10) || bf10 > 1e6) formatC(bf10, format = "e", digits = 2, decimal.mark = ",") else formatC(bf10, format = "f", digits = 3, decimal.mark = ",")
  txt_bf01 <- if (is.infinite(bf01) || bf01 > 1e6) formatC(bf01, format = "e", digits = 2, decimal.mark = ",") else formatC(bf01, format = "f", digits = 3, decimal.mark = ",")

  cat("  ",
      .bayes_pad(txt_bf10, wc["bf10"], "center"),
      .bayes_pad(txt_bf01, wc["bf01"], "center"),
      .bayes_pad(interp_bf$texto, wc["evid"], "left"),
      "\n", sep = "")
  cat(.bayes_sep(75), "\n\n")

  cat("  ESTIMATIVAS A POSTERIORI DO COEFICIENTE (\u03c1):\n")
  cat(.bayes_sep(75), "\n")
  wc_p <- c(med = 12, moda = 12, dp = 12, hpd = 25, pdir = 14)
  hdr_p <- paste0(.bayes_pad("Média Post.", wc_p["med"], "center"),
                  .bayes_pad("Mediana Post.", wc_p["moda"], "center"),
                  .bayes_pad("DP Post.", wc_p["dp"], "center"),
                  .bayes_pad(sprintf("HPD (%d%%)", as.integer(confianca * 100)), wc_p["hpd"], "center"),
                  .bayes_pad("P(\u03c1 > 0)", wc_p["pdir"], "center"))
  cat("  ", hdr_p, "\n", sep = "")
  cat(.bayes_sep(75), "\n")

  txt_hpd <- sprintf("[%s; %s]",
                     formatC(hpd_inf, format = "f", digits = decimais, decimal.mark = ","),
                     formatC(hpd_sup, format = "f", digits = decimais, decimal.mark = ","))

  cat("  ",
      .bayes_pad(formatC(media_post, format = "f", digits = decimais, decimal.mark = ","), wc_p["med"], "center"),
      .bayes_pad(formatC(mediana_post, format = "f", digits = decimais, decimal.mark = ","), wc_p["moda"], "center"),
      .bayes_pad(formatC(dp_post, format = "f", digits = decimais, decimal.mark = ","), wc_p["dp"], "center"),
      .bayes_pad(txt_hpd, wc_p["hpd"], "center"),
      .bayes_pad(.bayes_fmt_pct(prob_pos * 100, 1), wc_p["pdir"], "center"),
      "\n", sep = "")
  cat(.bayes_sep(75), "\n\n")

  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  if (bf10 >= 3) {
    cat(sprintf("  [\u2713] Os dados fornecem %s a favor da correlação (BF\u2081\u2080 = %s).\n",
                tolower(interp_bf$texto), txt_bf10))
    cat(sprintf("  A estimativa mais provável para o coeficiente \u03c1 é %s com intervalo de credibilidade de 95%% entre %s e %s.\n",
                formatC(mediana_post, format = "f", digits = decimais, decimal.mark = ","),
                formatC(hpd_inf, format = "f", digits = decimais, decimal.mark = ","),
                formatC(hpd_sup, format = "f", digits = decimais, decimal.mark = ",")))
  } else if (bf10 <= 1/3) {
    cat(sprintf("  Os dados fornecem %s a favor da ausência de associação linear (BF\u2080\u2081 = %s).\n",
                tolower(interp_bf$texto), txt_bf01))
  } else {
    cat(sprintf("  [!] Evidência inconclusiva / anedótica (BF\u2081\u2080 = %s; BF\u2080\u2081 = %s).\n",
                txt_bf10, txt_bf01))
    cat("  Recomenda-se coletar mais observações para definir a existência da associação.\n")
  }
  cat("───────────────────────────────────────────────────────────────────────────\n\n")

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot <- data.frame(rho = grid_rho, densidade = post_dens)
    df_hpd  <- df_plot[df_plot$rho >= hpd_inf & df_plot$rho <= hpd_sup, ]

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = rho, y = densidade)) +
      ggplot2::geom_ribbon(data = df_hpd, ggplot2::aes(ymin = 0, ymax = densidade),
                           fill = "#005EA5", alpha = 0.25) +
      ggplot2::geom_line(color = "#005EA5", linewidth = 1.2) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#D90429", linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = c(hpd_inf, hpd_sup), linetype = "dashed", color = "#005EA5", linewidth = 0.7) +
      ggplot2::geom_point(data = data.frame(x = mediana_post, y = max(post_dens)),
                          ggplot2::aes(x = x, y = y), color = "#005EA5", size = 3) +
      ggplot2::labs(
        title = "Distribuição a Posteriori da Correlação (\u03c1)",
        subtitle = sprintf("BF\u2081\u2080 = %s  |  Mediana Post. = %s  |  HPD 95%%: [%s; %s]",
                           txt_bf10,
                           formatC(mediana_post, format = "f", digits = decimais, decimal.mark = ","),
                           formatC(hpd_inf, format = "f", digits = decimais, decimal.mark = ","),
                           formatC(hpd_sup, format = "f", digits = decimais, decimal.mark = ",")),
        x = "Coeficiente de Correlação (\u03c1)",
        y = "Densidade a Posteriori"
      )

    if (exists("meu_tema", mode = "function")) {
      p <- p + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p <- p + tema_estatR(estilo = 2, grade = "dupla")
    } else {
      p <- p + ggplot2::theme_minimal()
    }

    suppressMessages(print(p))
  }

  res <- list(
    r_amostral = r_amostral,
    bf10 = bf10,
    bf01 = bf01,
    media_post = media_post,
    mediana_post = mediana_post,
    dp_post = dp_post,
    hpd = c(inf = hpd_inf, sup = hpd_sup),
    prob_positiva = prob_pos,
    prob_negativa = prob_neg
  )
  class(res) <- c("estatR_correlacao_bayesiana", "list")
  invisible(res)
}

#' @rdname teste_correlacao_bayesiano
#' @export
correlacao_bayesiana <- teste_correlacao_bayesiano


# ── 6. ANÁLISE DE VARIÂNCIA BAYESIANA (ANOVA ONE-WAY JZS) ────────────────────

#' @title ANOVA One-Way Bayesiana (Fator de Bayes JZS)
#' @description Realiza Análise de Variância (ANOVA) de um fator sob inferência bayesiana,
#'   calculando o Fator de Bayes JZS (Rouder et al., 2012 / Liang et al., 2008) comparando o modelo
#'   com efeito de grupo versus o modelo nulo, estimativas a posteriori das médias de cada grupo
#'   com intervalos de credibilidade HPD 95% e contrastes pareados post-hoc bayesianos.
#'
#' @param formula Fórmula do tipo \code{resposta ~ grupo} ou data frame quando usado com pipe (\code{|>}).
#' @param dados Data frame contendo as variáveis especificadas na fórmula quando não utilizado via pipe.
#' @param r_scale Escala da priori Cauchy para o tamanho do efeito padronizado (padrão \code{sqrt(2)/2} \eqn{\approx 0{,}707}).
#' @param confianca Nível de probabilidade do intervalo de credibilidade (padrão \code{0.95}).
#' @param decimais Número de casas decimais na exibição (padrão \code{2}).
#' @param grafico Lógico. Se TRUE (padrão), exibe gráfico com as médias a posteriori e intervalos HPD 95%.
#'
#' @return Objeto invisível da classe \code{"estatR_anova_bayesiana"} contendo \code{bf10}, \code{bf01},
#'   \code{grupos} e \code{contrastes}.
#' @export
anova_bayesiana <- function(formula, dados = NULL,
                            r_scale = sqrt(2)/2,
                            confianca = 0.95,
                            decimais = 2,
                            grafico = TRUE) {

  if (is.data.frame(formula) && (inherits(dados, "formula") || is.call(dados))) {
    tmp <- formula
    formula <- dados
    dados <- tmp
  }

  res_biv <- .resolver_dados_bivariados(formula, NULL, dados, env = parent.frame())
  y_raw      <- res_biv$y
  grupo_raw  <- res_biv$grupo
  nome_y     <- res_biv$nome_y
  nome_grupo <- res_biv$nome_grupo

  ok <- !is.na(y_raw) & !is.na(grupo_raw)
  y_raw <- y_raw[ok]
  grupo_raw <- as.factor(grupo_raw[ok])
  niveis <- levels(grupo_raw)
  k <- length(niveis)
  N <- length(y_raw)

  if (k < 2) stop("O fator 'grupo' deve ter pelo menos 2 níveis.", call. = FALSE)

  fit <- stats::aov(y_raw ~ grupo_raw)
  sm  <- summary(fit)[[1]]
  f_stat <- sm["grupo_raw", "F value"]
  df1    <- sm["grupo_raw", "Df"]
  df2    <- sm["Residuals", "Df"]
  r2     <- (df1 * f_stat) / (df1 * f_stat + df2)
  s_p    <- sqrt(sm["Residuals", "Mean Sq"])

  integrando <- function(g) {
    prior_g <- (r_scale / sqrt(2 * pi)) * (g^(-1.5)) * exp(- (r_scale^2) / (2 * g))
    term1 <- (1 + g)^((N - k) / 2)
    term2 <- (1 + g * (1 - r2))^(- (N - 1) / 2)
    term1 * term2 * prior_g
  }
  bf10 <- tryCatch(stats::integrate(integrando, lower = 0, upper = Inf)$value, error = function(e) Inf)
  bf01 <- if (bf10 > 0) 1 / bf10 else Inf

  tcrit <- stats::qt(1 - (1 - confianca) / 2, df = df2)
  tab_grupos <- list()
  for (niv in niveis) {
    y_g <- y_raw[grupo_raw == niv]
    n_g <- length(y_g)
    m_obs <- mean(y_g)
    se_g  <- s_p / sqrt(n_g)
    hpd_inf <- m_obs - tcrit * se_g
    hpd_sup <- m_obs + tcrit * se_g
    tab_grupos[[niv]] <- data.frame(
      Grupo = niv, N = n_g, Media_Obs = m_obs,
      Media_Post = m_obs, DP_Post = se_g,
      HPD_Inf = hpd_inf, HPD_Sup = hpd_sup,
      stringsAsFactors = FALSE
    )
  }
  df_grupos <- do.call(rbind, tab_grupos)

  pares <- utils::combn(niveis, 2, simplify = FALSE)
  tab_pares <- list()
  for (par in pares) {
    g1 <- par[1]; g2 <- par[2]
    n1 <- length(y_raw[grupo_raw == g1])
    n2 <- length(y_raw[grupo_raw == g2])
    diff_m <- mean(y_raw[grupo_raw == g1]) - mean(y_raw[grupo_raw == g2])
    se_diff <- s_p * sqrt(1/n1 + 1/n2)
    tab_pares[[length(tab_pares) + 1]] <- data.frame(
      Contraste = paste0(g1, " - ", g2),
      Dif_Media = diff_m,
      DP_Post   = se_diff,
      HPD_Inf   = diff_m - tcrit * se_diff,
      HPD_Sup   = diff_m + tcrit * se_diff,
      stringsAsFactors = FALSE
    )
  }
  df_pares <- do.call(rbind, tab_pares)
  interp_bf <- .bayes_interpretar_bf(bf10)

  cat("\n── ANÁLISE DE VARIÂNCIA BAYESIANA (ANOVA ONE-WAY) ─────────────────────────\n\n")
  cat(sprintf("  Variável Resposta : %s\n", nome_y))
  cat(sprintf("  Fator de Grupo    : %s (%d níveis: %s)\n", nome_grupo, k, paste(niveis, collapse = ", ")))
  cat(sprintf("  Amostra Total     : N = %d observações\n\n", N))

  cat(sprintf("  H0: Médias populacionais iguais (\u03bc\u2081 = \u03bc\u2082 = ... = \u03bc_k)\n"))
  cat(sprintf("  H1: Efeito do fator de grupo presente (\u03bc_i diferem)\n\n"))

  cat(.bayes_sep(75), "\n")
  wc <- c(bf10 = 16, bf01 = 16, evid = 40)
  hdr <- paste0(.bayes_pad("Fator Bayes BF\u2081\u2080", wc["bf10"], "center"),
                .bayes_pad("Fator Bayes BF\u2080\u2081", wc["bf01"], "center"),
                .bayes_pad("Força da Evidência", wc["evid"], "left"))
  cat("  ", hdr, "\n", sep = "")
  cat(.bayes_sep(75), "\n")

  txt_bf10 <- if (is.infinite(bf10) || bf10 > 1e6) formatC(bf10, format = "e", digits = 2, decimal.mark = ",") else formatC(bf10, format = "f", digits = 3, decimal.mark = ",")
  txt_bf01 <- if (is.infinite(bf01) || bf01 > 1e6) formatC(bf01, format = "e", digits = 2, decimal.mark = ",") else formatC(bf01, format = "f", digits = 3, decimal.mark = ",")

  cat("  ",
      .bayes_pad(txt_bf10, wc["bf10"], "center"),
      .bayes_pad(txt_bf01, wc["bf01"], "center"),
      .bayes_pad(interp_bf$texto, wc["evid"], "left"),
      "\n", sep = "")
  cat(.bayes_sep(75), "\n\n")

  cat("  ESTIMATIVAS A POSTERIORI DAS MÉDIAS POR GRUPO:\n")
  cat(.bayes_sep(75), "\n")
  wc_g <- c(grp = 15, n = 8, med = 14, dp = 12, hpd = 24)
  hdr_g <- paste0(.bayes_pad("Grupo", wc_g["grp"], "left"),
                  .bayes_pad("N", wc_g["n"], "center"),
                  .bayes_pad("Média Post.", wc_g["med"], "center"),
                  .bayes_pad("DP Post.", wc_g["dp"], "center"),
                  .bayes_pad(sprintf("HPD (%d%%)", as.integer(confianca * 100)), wc_g["hpd"], "center"))
  cat("  ", hdr_g, "\n", sep = "")
  cat(.bayes_sep(75), "\n")

  for (i in seq_len(nrow(df_grupos))) {
    row <- df_grupos[i, ]
    txt_hpd_g <- sprintf("[%s; %s]",
                         formatC(row$HPD_Inf, format = "f", digits = decimais, decimal.mark = ","),
                         formatC(row$HPD_Sup, format = "f", digits = decimais, decimal.mark = ","))
    cat("  ",
        .bayes_pad(row$Grupo, wc_g["grp"], "left"),
        .bayes_pad(as.character(row$N), wc_g["n"], "center"),
        .bayes_pad(formatC(row$Media_Post, format = "f", digits = decimais, decimal.mark = ","), wc_g["med"], "center"),
        .bayes_pad(formatC(row$DP_Post, format = "f", digits = decimais, decimal.mark = ","), wc_g["dp"], "center"),
        .bayes_pad(txt_hpd_g, wc_g["hpd"], "center"),
        "\n", sep = "")
  }
  cat(.bayes_sep(75), "\n\n")

  cat("  CONTRASTES POST-HOC BAYESIANOS (DIFERENÇAS ENTRE GRUPOS):\n")
  cat(.bayes_sep(75), "\n")
  wc_c <- c(par = 22, dif = 14, dp = 12, hpd = 24)
  hdr_c <- paste0(.bayes_pad("Comparação", wc_c["par"], "left"),
                  .bayes_pad("Dif. Média", wc_c["dif"], "center"),
                  .bayes_pad("DP Post.", wc_c["dp"], "center"),
                  .bayes_pad(sprintf("HPD (%d%%)", as.integer(confianca * 100)), wc_c["hpd"], "center"))
  cat("  ", hdr_c, "\n", sep = "")
  cat(.bayes_sep(75), "\n")

  for (i in seq_len(nrow(df_pares))) {
    row <- df_pares[i, ]
    txt_hpd_c <- sprintf("[%s; %s]",
                         formatC(row$HPD_Inf, format = "f", digits = decimais, decimal.mark = ","),
                         formatC(row$HPD_Sup, format = "f", digits = decimais, decimal.mark = ","))
    cat("  ",
        .bayes_pad(row$Contraste, wc_c["par"], "left"),
        .bayes_pad(formatC(row$Dif_Media, format = "f", digits = decimais, decimal.mark = ","), wc_c["dif"], "center"),
        .bayes_pad(formatC(row$DP_Post, format = "f", digits = decimais, decimal.mark = ","), wc_c["dp"], "center"),
        .bayes_pad(txt_hpd_c, wc_c["hpd"], "center"),
        "\n", sep = "")
  }
  cat(.bayes_sep(75), "\n\n")

  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  if (bf10 >= 3) {
    cat(sprintf("  [\u2713] Há evidências a favor do efeito de grupo (BF\u2081\u2080 = %s, %s).\n",
                txt_bf10, tolower(interp_bf$texto)))
    cat("  As médias diferem substantivamente entre as condições testadas.\n")
  } else if (bf10 <= 1/3) {
    cat(sprintf("  Os dados suportam a hipótese nula de igualdade das médias (BF\u2080\u2081 = %s, %s).\n",
                txt_bf01, tolower(interp_bf$texto)))
  } else {
    cat(sprintf("  [!] Evidência inconclusiva / anedótica (BF\u2081\u2080 = %s).\n", txt_bf10))
  }
  cat("───────────────────────────────────────────────────────────────────────────\n\n")

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df_grupos, ggplot2::aes(x = Grupo, y = Media_Post)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = HPD_Inf, ymax = HPD_Sup),
                             width = 0.2, color = "#005EA5", linewidth = 0.9) +
      ggplot2::geom_point(size = 3.5, color = "#005EA5", fill = "white", shape = 21, stroke = 1.3) +
      ggplot2::labs(
        title = "Médias a Posteriori e Intervalos de Credibilidade (HPD 95%)",
        subtitle = sprintf("BF\u2081\u2080 = %s (%s)", txt_bf10, interp_bf$texto),
        x = nome_grupo,
        y = sprintf("Média a Posteriori de %s", nome_y)
      )

    if (exists("meu_tema", mode = "function")) {
      p <- p + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p <- p + tema_estatR(estilo = 2, grade = "dupla")
    } else {
      p <- p + ggplot2::theme_minimal()
    }
    suppressMessages(print(p))
  }

  res <- list(
    bf10 = bf10,
    bf01 = bf01,
    grupos = df_grupos,
    contrastes = df_pares,
    r2 = r2
  )
  class(res) <- c("estatR_anova_bayesiana", "list")
  invisible(res)
}


# ── 7. REGRESSÃO LINEAR BAYESIANA (JZS / G-PRIORI) ───────────────────────────

#' @title Regressão Linear Bayesiana (JZS / G-Priori de Zellner)
#' @description Ajusta modelos lineares múltiplos sob a ótica bayesiana utilizando g-priori de Zellner
#'   e priori de Jeffreys-Zellner-Siow (Liang et al., 2008), calculando o Fator de Bayes \eqn{BF_{10}}
#'   do modelo global versus o modelo nulo, estimativas a posteriori dos coeficientes com encolhimento
#'   empírico (EB shrinkage), intervalos de credibilidade HPD 95% e probabilidades direcionais.
#'
#' @param formula Fórmula do modelo (ex: \code{y ~ x1 + x2}) ou data frame quando usado com pipe (\code{|>}).
#' @param dados Data frame contendo as variáveis especificadas na fórmula.
#' @param r_scale Parâmetro de escala da priori Cauchy JZS para regressão (padrão \code{0.354}).
#' @param confianca Nível de probabilidade do intervalo de credibilidade (padrão \code{0.95}).
#' @param decimais Número de casas decimais na exibição (padrão \code{3}).
#' @param grafico Lógico. Se TRUE (padrão), exibe gráfico com as estimativas dos coeficientes e intervalos HPD 95%.
#'
#' @return Objeto invisível da classe \code{"estatR_regressao_bayesiana"} contendo \code{bf10}, \code{bf01},
#'   \code{r2} e a tabela de \code{coeficientes}.
#' @export
regressao_bayesiana <- function(formula, dados = NULL,
                                r_scale = 0.354,
                                confianca = 0.95,
                                decimais = 3,
                                grafico = TRUE) {

  if (is.data.frame(formula) && (inherits(dados, "formula") || is.call(dados))) {
    tmp <- formula
    formula <- dados
    dados <- tmp
  }

  mf <- stats::model.frame(formula, data = dados)
  y  <- stats::model.response(mf)
  X  <- stats::model.matrix(formula, data = mf)
  termos <- colnames(X)
  nome_y <- names(mf)[1]
  n <- nrow(X)
  p <- ncol(X) - 1

  if (p < 1) stop("O modelo deve conter pelo menos uma variável explicativa.", call. = FALSE)

  fit_ols <- stats::lm.fit(X, y)
  residuos <- fit_ols$residuals
  ss_res <- sum(residuos^2)
  ss_tot <- sum((y - mean(y))^2)
  r2 <- 1 - (ss_res / ss_tot)
  s2 <- ss_res / (n - p - 1)
  df_res <- n - p - 1

  integrando <- function(g) {
    prior_g <- (r_scale / sqrt(2 * pi)) * (g^(-1.5)) * exp(- (r_scale^2) / (2 * g))
    term1 <- (1 + g)^((n - 1 - p) / 2)
    term2 <- (1 + g * (1 - r2))^(- (n - 1) / 2)
    term1 * term2 * prior_g
  }
  bf10 <- tryCatch(stats::integrate(integrando, lower = 0, upper = Inf)$value, error = function(e) Inf)
  bf01 <- if (bf10 > 0) 1 / bf10 else Inf

  f_stat <- (r2 / p) / ((1 - r2) / df_res)
  g_hat <- max(0, f_stat - 1)
  fator_enc <- if (g_hat > 0) g_hat / (1 + g_hat) else 0.999
  if (is.nan(fator_enc) || fator_enc <= 0) fator_enc <- 0.999

  beta_ols <- fit_ols$coefficients
  beta_post <- beta_ols
  if (p > 0) {
    beta_post[-1] <- beta_ols[-1] * fator_enc
  }

  xtx_inv <- tryCatch(solve(t(X) %*% X), error = function(e) diag(p + 1) * 1e-4)
  se_beta_post <- sqrt(diag(xtx_inv) * s2 * fator_enc)
  se_beta_post[1] <- sqrt(diag(xtx_inv)[1] * s2)

  tcrit <- stats::qt(1 - (1 - confianca) / 2, df = df_res)
  hpd_inf <- beta_post - tcrit * se_beta_post
  hpd_sup <- beta_post + tcrit * se_beta_post
  prob_pos <- 1 - stats::pt((0 - beta_post) / se_beta_post, df = df_res)

  df_coef <- data.frame(
    Parametro = termos,
    Media_Post = beta_post,
    DP_Post = se_beta_post,
    HPD_Inf = hpd_inf,
    HPD_Sup = hpd_sup,
    Prob_Pos = prob_pos,
    stringsAsFactors = FALSE
  )

  interp_bf <- .bayes_interpretar_bf(bf10)

  cat("\n── REGRESSÃO LINEAR BAYESIANA (JZS / G-PRIORI) ────────────────────────────\n\n")
  cat(sprintf("  Variável Resposta   : %s\n", nome_y))
  cat(sprintf("  Preditores (p)      : %d preditor(es)\n", p))
  cat(sprintf("  Tamanho Amostral    : N = %d observações\n", n))
  cat(sprintf("  R² Amostral         : %s\n\n", formatC(r2, format = "f", digits = 4, decimal.mark = ",")))

  cat(sprintf("  H0: Modelo nulo apenas com intercepto (R² = 0)\n"))
  cat(sprintf("  H1: Modelo com preditores lineares (R² > 0)\n\n"))

  cat(.bayes_sep(75), "\n")
  wc <- c(bf10 = 16, bf01 = 16, evid = 40)
  hdr <- paste0(.bayes_pad("Fator Bayes BF\u2081\u2080", wc["bf10"], "center"),
                .bayes_pad("Fator Bayes BF\u2080\u2081", wc["bf01"], "center"),
                .bayes_pad("Força da Evidência", wc["evid"], "left"))
  cat("  ", hdr, "\n", sep = "")
  cat(.bayes_sep(75), "\n")

  txt_bf10 <- if (is.infinite(bf10) || bf10 > 1e6) formatC(bf10, format = "e", digits = 2, decimal.mark = ",") else formatC(bf10, format = "f", digits = 3, decimal.mark = ",")
  txt_bf01 <- if (is.infinite(bf01) || bf01 > 1e6) formatC(bf01, format = "e", digits = 2, decimal.mark = ",") else formatC(bf01, format = "f", digits = 3, decimal.mark = ",")

  cat("  ",
      .bayes_pad(txt_bf10, wc["bf10"], "center"),
      .bayes_pad(txt_bf01, wc["bf01"], "center"),
      .bayes_pad(interp_bf$texto, wc["evid"], "left"),
      "\n", sep = "")
  cat(.bayes_sep(75), "\n\n")

  cat("  ESTIMATIVAS A POSTERIORI DOS COEFICIENTES (\u03b2):\n")
  cat(.bayes_sep(75), "\n")
  wc_b <- c(par = 18, est = 13, dp = 12, hpd = 22, pdir = 10)
  hdr_b <- paste0(.bayes_pad("Parâmetro", wc_b["par"], "left"),
                  .bayes_pad("Média Post.", wc_b["est"], "center"),
                  .bayes_pad("DP Post.", wc_b["dp"], "center"),
                  .bayes_pad(sprintf("HPD (%d%%)", as.integer(confianca * 100)), wc_b["hpd"], "center"),
                  .bayes_pad("P(\u03b2 > 0)", wc_b["pdir"], "center"))
  cat("  ", hdr_b, "\n", sep = "")
  cat(.bayes_sep(75), "\n")

  for (i in seq_len(nrow(df_coef))) {
    row <- df_coef[i, ]
    txt_hpd_b <- sprintf("[%s; %s]",
                         formatC(row$HPD_Inf, format = "f", digits = decimais, decimal.mark = ","),
                         formatC(row$HPD_Sup, format = "f", digits = decimais, decimal.mark = ","))
    cat("  ",
        .bayes_pad(row$Parametro, wc_b["par"], "left"),
        .bayes_pad(formatC(row$Media_Post, format = "f", digits = decimais, decimal.mark = ","), wc_b["est"], "center"),
        .bayes_pad(formatC(row$DP_Post, format = "f", digits = decimais, decimal.mark = ","), wc_b["dp"], "center"),
        .bayes_pad(txt_hpd_b, wc_b["hpd"], "center"),
        .bayes_pad(.bayes_fmt_pct(row$Prob_Pos * 100, 1), wc_b["pdir"], "center"),
        "\n", sep = "")
  }
  cat(.bayes_sep(75), "\n\n")

  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  if (bf10 >= 3) {
    cat(sprintf("  [\u2713] O modelo global é amplamente sustentado pelos dados (BF\u2081\u2080 = %s, %s).\n",
                txt_bf10, tolower(interp_bf$texto)))
    cat(sprintf("  O conjunto de variáveis explica %s da variação de %s.\n",
                .bayes_fmt_pct(r2 * 100, 1), nome_y))
  } else if (bf10 <= 1/3) {
    cat(sprintf("  Os dados suportam o modelo nulo em detrimento dos preditores (BF\u2080\u2081 = %s).\n", txt_bf01))
  } else {
    cat(sprintf("  [!] Evidência inconclusiva / anedótica para o modelo global (BF\u2081\u2080 = %s).\n", txt_bf10))
  }
  cat("───────────────────────────────────────────────────────────────────────────\n\n")

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot <- df_coef[df_coef$Parametro != "(Intercept)", ]
    if (nrow(df_plot) > 0) {
      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Media_Post, y = stats::reorder(Parametro, Media_Post))) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = HPD_Inf, xmax = HPD_Sup), height = 0.2, color = "#005EA5", linewidth = 0.9) +
        ggplot2::geom_point(size = 3.5, color = "#005EA5", fill = "white", shape = 21, stroke = 1.3) +
        ggplot2::labs(
          title = "Coeficientes a Posteriori e Intervalos de Credibilidade (HPD 95%)",
          subtitle = sprintf("BF\u2081\u2080 = %s  |  R² = %s", txt_bf10, formatC(r2, format = "f", digits = 3, decimal.mark = ",")),
          x = "Estimativa a Posteriori (\u03b2)",
          y = NULL
        )

      if (exists("meu_tema", mode = "function")) {
        p <- p + meu_tema(grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        p <- p + tema_estatR(estilo = 2, grade = "dupla")
      } else {
        p <- p + ggplot2::theme_minimal()
      }
      suppressMessages(print(p))
    }
  }

  res <- list(
    bf10 = bf10,
    bf01 = bf01,
    r2 = r2,
    coeficientes = df_coef
  )
  class(res) <- c("estatR_regressao_bayesiana", "list")
  invisible(res)
}
