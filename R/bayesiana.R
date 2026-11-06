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
