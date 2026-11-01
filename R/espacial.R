# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: ESTATÍSTICA ESPACIAL
# Autocorrelação Espacial (Moran Global e LISA) e Modelos de Regressão Espacial
# ─────────────────────────────────────────────────────────────────────────────

# ── HELPER: CONSTRUÇÃO E PADRONIZAÇÃO DA MATRIZ DE PESOS ESPACIAIS ───────────

#' @noRd
.construir_matriz_pesos <- function(coords, k = 5, padronizacao = "linha") {
  coords <- as.matrix(coords)
  n <- nrow(coords)

  # Distância Euclidiana pareada
  dists <- as.matrix(stats::dist(coords))
  diag(dists) <- Inf

  W <- matrix(0, nrow = n, ncol = n)
  k_efetivo <- min(k, n - 1)

  for (i in seq_len(n)) {
    vizinhos <- order(dists[i, ])[seq_len(k_efetivo)]
    W[i, vizinhos] <- 1
  }

  # Padronização
  if (padronizacao == "linha") {
    soma_linhas <- rowSums(W)
    soma_linhas[soma_linhas == 0] <- 1
    W <- W / soma_linhas
  } else if (padronizacao == "global") {
    soma_total <- sum(W)
    if (soma_total > 0) W <- W / soma_total
  } else if (padronizacao == "binaria") {
    W <- (W > 0) * 1
  }
  # "nenhuma" mantém como está

  W
}


# ── 1. AUTOCORRELAÇÃO ESPACIAL (MORAN GLOBAL E LISA) ─────────────────────────

#' @title Análise de Autocorrelação Espacial Global e Local (Moran e LISA)
#' @description Calcula o Índice de Moran Global e o Índice de Moran Local (LISA)
#'   para identificar dependência espacial e padrões de agrupamento (clusters Alto-Alto,
#'   Baixo-Baixo e outliers espaciais).
#'
#' @param dados Um \code{data.frame} com as variáveis de interesse.
#' @param variavel Nome da variável numérica a analisar (símbolo ou string).
#' @param coordenadas Matriz ou data frame com coordenadas espaciais (X e Y, ou Long e Lat).
#'   Se omitido, tenta detectar colunas com nomes como \code{x, y, long, lat}.
#' @param matriz_pesos Matriz de vizinhança espacial \eqn{W} opcional (quadrada \eqn{n \times n}).
#' @param k_vizinhos Número de vizinhos mais próximos para construir a matriz caso \code{matriz_pesos} não seja informada.
#' @param n_permutacoes Número de permutações Monte Carlo para o pseudo-p-valor (padrão: 999).
#' @param padronizacao Método de padronização da matriz de pesos: \code{"linha"} (padrão),
#'   \code{"binaria"}, \code{"global"} ou \code{"nenhuma"}.
#' @param grafico Lógico. Se TRUE (padrão), gera o Diagrama de Dispersão de Moran.
#'
#' @return Objeto invisível da classe \code{estatR_autocorrelacao_espacial}.
#' @export
analise_autocorrelacao <- function(dados, variavel,
                                   coordenadas = NULL,
                                   matriz_pesos = NULL,
                                   k_vizinhos = 5,
                                   n_permutacoes = 999,
                                   padronizacao = c("linha", "binaria", "global", "nenhuma"),
                                   grafico = TRUE) {

  padronizacao <- match.arg(padronizacao)

  if (!is.data.frame(dados)) {
    stop("O argumento 'dados' deve ser um data.frame.", call. = FALSE)
  }

  var_sub  <- substitute(variavel)
  var_nome <- if (is.symbol(var_sub)) as.character(var_sub) else if (is.character(variavel)) variavel else NULL
  if (is.null(var_nome) || !var_nome %in% names(dados)) {
    stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada nos dados.", var_nome), call. = FALSE)
  }

  y <- as.numeric(dados[[var_nome]])
  idx_val <- !is.na(y)
  y <- y[idx_val]
  n <- length(y)

  if (n < 5) {
    stop("S\u00e3o necess\u00e1rias pelo menos 5 observa\u00e7\u00f5es v\u00e1lidas.", call. = FALSE)
  }

  # ── Matriz de Pesos W ───────────────────────────────────────────────────────
  if (!is.null(matriz_pesos)) {
    W <- as.matrix(matriz_pesos)[idx_val, idx_val, drop = FALSE]
  } else {
    coords <- NULL
    if (!is.null(coordenadas)) {
      coords <- as.matrix(coordenadas)[idx_val, , drop = FALSE]
    } else {
      # Detecção automática de coordenadas
      cand_cols <- grep("^(x|y|long|lat|lon|latitude|longitude)$", names(dados), ignore.case = TRUE, value = TRUE)
      if (length(cand_cols) >= 2) {
        coords <- as.matrix(dados[idx_val, cand_cols[1:2], drop = FALSE])
      } else {
        # Se não há coordenadas, cria grade sequencial 1D como fallback
        coords <- matrix(c(seq_len(n), rep(1, n)), ncol = 2)
      }
    }
    W <- .construir_matriz_pesos(coords, k = k_vizinhos, padronizacao = padronizacao)
  }

  # ── Cálculo do Índice de Moran Global ───────────────────────────────────────
  z <- y - mean(y)
  s2 <- sum(z^2) / n
  S0 <- sum(W)

  # Defasagem espacial (Spatial Lag)
  Wz <- as.numeric(W %*% z)

  moran_i <- (n / S0) * (sum(z * Wz) / sum(z^2))
  esperanca_i <- -1 / (n - 1)

  # Permutação Monte Carlo
  perm_i <- numeric(n_permutacoes)
  for (p in seq_len(n_permutacoes)) {
    z_perm <- sample(z)
    Wz_perm <- as.numeric(W %*% z_perm)
    perm_i[p] <- (n / S0) * (sum(z_perm * Wz_perm) / sum(z_perm^2))
  }

  var_perm <- stats::var(perm_i)
  ep_i     <- sqrt(var_perm)
  z_score  <- (moran_i - esperanca_i) / ep_i
  p_pseudo <- (sum(abs(perm_i) >= abs(moran_i)) + 1) / (n_permutacoes + 1)

  # ── Cálculo do Moran Local (LISA) e Quadrantes ─────────────────────────────
  lisa_i <- numeric(n)
  m2 <- sum(z^2) / n
  for (i in seq_len(n)) {
    lisa_i[i] <- (z[i] / m2) * sum(W[i, ] * z)
  }

  # Classificação dos Quadrantes
  z_pad  <- (y - mean(y)) / stats::sd(y)
  Wz_pad <- as.numeric(W %*% z_pad)

  quadrante <- character(n)
  for (i in seq_len(n)) {
    if (z_pad[i] >= 0 && Wz_pad[i] >= 0) {
      quadrante[i] <- "Alto-Alto (Q1)"
    } else if (z_pad[i] < 0 && Wz_pad[i] < 0) {
      quadrante[i] <- "Baixo-Baixo (Q2)"
    } else if (z_pad[i] >= 0 && Wz_pad[i] < 0) {
      quadrante[i] <- "Alto-Baixo (Q3)"
    } else {
      quadrante[i] <- "Baixo-Alto (Q4)"
    }
  }

  tab_quad <- table(quadrante)
  cores_quad <- c(
    "Alto-Alto (Q1)"   = "#D90429", # Vermelho quente
    "Baixo-Baixo (Q2)" = "#2B5C8F", # Azul frio
    "Alto-Baixo (Q3)"  = "#E5989B", # Rosa transicao
    "Baixo-Alto (Q4)"  = "#85B4D4"  # Azul claro transicao
  )

  # ── Gráfico: Diagrama de Dispersão de Moran ─────────────────────────────────
  p_moran <- NULL
  if (isTRUE(grafico)) {
    df_plot <- data.frame(
      Z_Padronizado = z_pad,
      Spatial_Lag   = Wz_pad,
      Quadrante     = factor(quadrante, levels = c("Alto-Alto (Q1)", "Baixo-Baixo (Q2)",
                                                   "Alto-Baixo (Q3)", "Baixo-Alto (Q4)")),
      stringsAsFactors = FALSE
    )

    cores_quad <- c(
      "Alto-Alto (Q1)"   = "#D90429", # Vermelho quente
      "Baixo-Baixo (Q2)" = "#2B5C8F", # Azul frio
      "Alto-Baixo (Q3)"  = "#E5989B", # Rosa transição
      "Baixo-Alto (Q4)"  = "#85B4D4"  # Azul claro transição
    )

    p_moran <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Z_Padronizado, y = Spatial_Lag)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      ggplot2::geom_point(ggplot2::aes(color = Quadrante), size = 2.8, alpha = 0.8) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#111111", fill = "gray70", alpha = 0.25, linewidth = 0.9) +
      ggplot2::scale_color_manual(values = cores_quad) +
      ggplot2::labs(
        title    = sprintf("Diagrama de Dispers\u00e3o de Moran: %s", var_nome),
        subtitle = sprintf("I de Moran: %.3f (z = %.2f, p = %s)", moran_i, z_score, .fmt_p(p_pseudo)),
        x        = sprintf("Valor Padronizado (%s)", var_nome),
        y        = "Defasagem Espacial (W \u00d7 z)"
      )

    if (exists("meu_tema", mode = "function")) {
      p_moran <- p_moran + meu_tema()
    } else {
      p_moran <- p_moran + ggplot2::theme_minimal()
    }
    print(p_moran)
  }

  # ── Saída no Console ───────────────────────────────────────────────────────
  pad <- function(s, w, align = "left") {
    s <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right") return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }
  sep_l <- function(w) paste0("  ", strrep("\u2500", w))

  .print_titulo("AUTOCORRELA\u00c7\u00c3O ESPACIAL (MORAN GLOBAL E LISA)")

  cat(sprintf("  Vari\u00e1vel Analisada : %s\n", var_nome))
  cat(sprintf("  N\u00ba Observa\u00e7\u00f5es      : %d\n", n))
  cat(sprintf("  Matriz de Pesos     : Padroniza\u00e7\u00e3o por %s\n", padronizacao))
  cat(sprintf("  Permuta\u00e7\u00f5es Monte C.: %d replica\u00e7\u00f5es\n\n", n_permutacoes))

  .print_topico("\u00cdNDICE DE MORAN GLOBAL")
  w_m <- c(ind = 18, esp = 14, ep = 12, z = 10, p = 14)
  w_tot <- sum(w_m) + length(w_m) - 1

  hdr_m <- paste0(
    pad("I de Moran",  w_m["ind"], "center"),
    pad("Esperan\u00e7a",  w_m["esp"], "center"),
    pad("Erro Pad.",  w_m["ep"],  "center"),
    pad("Estat. Z",   w_m["z"],   "center"),
    pad("p-valor",    w_m["p"],   "center")
  )
  cat(sep_l(w_tot), "\n")
  cat("  ", hdr_m, "\n", sep = "")
  cat(sep_l(w_tot), "\n")
  cat("  ",
      pad(.fmt_num(moran_i, 4),     w_m["ind"], "center"),
      pad(.fmt_num(esperanca_i, 4), w_m["esp"], "center"),
      pad(.fmt_num(ep_i, 4),        w_m["ep"],  "center"),
      pad(.fmt_num(z_score, 2),     w_m["z"],   "center"),
      pad(.fmt_p(p_pseudo),         w_m["p"],   "center"),
      "\n", sep = "")
  cat(sep_l(w_tot), "\n\n")

  cat("  Notas:\n")
  cat("  Signific\u00e2ncia: *** p < 0,001   ** p < 0,01   * p < 0,05\n")
  cat("  I > E(I) indica autocorrela\u00e7\u00e3o positiva (agrupamento/clusters).\n")
  cat("  I < E(I) indica autocorrela\u00e7\u00e3o negativa (dispers\u00e3o espacial).\n\n")

  .print_topico("DISTRIBUI\u00c7\u00c3O DOS QUADRANTES DE MORAN LOCAL (LISA)")
  w_q <- c(quad = 22, n = 12, pct = 14)
  w_tot_q <- sum(w_q) + length(w_q) - 1

  hdr_q <- paste0(
    pad("Quadrante",        w_q["quad"], "left"),
    pad("Quantidade",       w_q["n"],    "center"),
    pad("Propor\u00e7\u00e3o (%)",  w_q["pct"],  "center")
  )
  cat(sep_l(w_tot_q), "\n")
  cat("  ", hdr_q, "\n", sep = "")
  cat(sep_l(w_tot_q), "\n")

  for (qd in names(cores_quad)) {
    cnt <- if (qd %in% names(tab_quad)) as.numeric(tab_quad[qd]) else 0
    pct <- (cnt / n) * 100
    cat("  ",
        pad(qd,                   w_q["quad"], "left"),
        pad(cnt,                  w_q["n"],    "center"),
        pad(.fmt_pct(pct / 100),  w_q["pct"],  "center"),
        "\n", sep = "")
  }
  cat(sep_l(w_tot_q), "\n\n")

  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_pseudo < 0.05) {
    if (moran_i > esperanca_i) {
      cat(sprintf("  Foi detectada AUTOCORRELA\u00c7\u00c3O ESPACIAL POSITIVA estatisticamente significativa\n  (I = %.4f; z = %.2f; p = %s).\n",
                  moran_i, z_score, .fmt_p(p_pseudo)))
      cat("  Regiões com valores elevados tendem a ser vizinhas de regi\u00f5es com valores elevados,\n  e áreas com baixos valores agrupam-se entre si (clusters espaciais definidos).\n\n")
    } else {
      cat(sprintf("  Foi detectada AUTOCORRELA\u00c7\u00c3O ESPACIAL NEGATIVA estatisticamente significativa\n  (I = %.4f; z = %.2f; p = %s).\n",
                  moran_i, z_score, .fmt_p(p_pseudo)))
      cat("  O padr\u00e3o espacial revela dispers\u00e3o sistem\u00e1tica ou altern\u00e2ncia em formato de tabuleiro de xadrez.\n\n")
    }
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancia estat\u00edstica de autocorrela\u00e7\u00e3o espacial (I = %.4f; p = %s).\n",
                moran_i, .fmt_p(p_pseudo)))
    cat("  A distribui\u00e7\u00e3o espacial da vari\u00e1vel \u00e9 compativel com um processo puramente aleat\u00f3rio.\n\n")
  }

  .print_rodape()

  df_lisa <- data.frame(
    Indice_LISA = lisa_i,
    Z_Padronizado = z_pad,
    Spatial_Lag = Wz_pad,
    Quadrante   = quadrante,
    stringsAsFactors = FALSE
  )

  res <- list(
    moran_global  = list(I = moran_i, Esperanca = esperanca_i, EP = ep_i, Z = z_score, p_valor = p_pseudo),
    lisa          = df_lisa,
    quadrantes    = tab_quad,
    matriz_pesos  = W,
    grafico       = p_moran,
    variavel      = var_nome
  )
  class(res) <- c("estatR_autocorrelacao_espacial", "list")
  invisible(res)
}


# ── 2. REGRESSÃO ESPACIAL (SAR E SEM) ────────────────────────────────────────

#' @title Regressão Espacial: Modelos Autorregressivo (SAR) e de Erro (SEM)
#' @description Ajusta modelos lineares com dependência espacial incorporada:
#'   modelo autoregressivo na resposta (SAR / Spatial Lag) ou modelo com correlação
#'   nos erros (SEM / Spatial Error).
#'
#' @param formula Fórmula linear do modelo (ex: \code{y ~ x1 + x2}).
#' @param dados Data frame com os dados da análise.
#' @param coordenadas Matriz ou data frame com coordenadas espaciais.
#' @param matriz_pesos Matriz de vizinhança espacial \eqn{W} opcional.
#' @param modelo Tipo de modelo: \code{"sar"} (Spatial Lag) ou \code{"sem"} (Spatial Error).
#' @param k_vizinhos Número de vizinhos mais próximos caso \code{matriz_pesos} não seja fornecida.
#' @param padronizacao Método de padronização da matriz \eqn{W}: \code{"linha"} (padrão),
#'   \code{"binaria"}, \code{"global"} ou \code{"nenhuma"}.
#' @param ... Argumentos adicionais.
#'
#' @return Objeto invisível da classe \code{estatR_regressao_espacial}.
#' @export
regressao_espacial <- function(formula, dados,
                               coordenadas = NULL,
                               matriz_pesos = NULL,
                               modelo = c("sar", "sem"),
                               k_vizinhos = 5,
                               padronizacao = c("linha", "binaria", "global", "nenhuma"),
                               ...) {

  modelo <- match.arg(modelo)
  padronizacao <- match.arg(padronizacao)

  if (!is.data.frame(dados)) {
    stop("O argumento 'dados' deve ser um data.frame.", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = dados)
  y  <- stats::model.response(mf)
  X  <- stats::model.matrix(formula, data = mf)
  n  <- length(y)

  # Matriz de pesos W
  if (!is.null(matriz_pesos)) {
    W <- as.matrix(matriz_pesos)
  } else {
    coords <- NULL
    if (!is.null(coordenadas)) {
      coords <- as.matrix(coordenadas)
    } else {
      cand_cols <- grep("^(x|y|long|lat|lon|latitude|longitude)$", names(dados), ignore.case = TRUE, value = TRUE)
      if (length(cand_cols) >= 2) {
        coords <- as.matrix(dados[, cand_cols[1:2], drop = FALSE])
      } else {
        coords <- matrix(c(seq_len(n), rep(1, n)), ncol = 2)
      }
    }
    W <- .construir_matriz_pesos(coords, k = k_vizinhos, padronizacao = padronizacao)
  }

  # Defasagem espacial Wy
  Wy <- as.numeric(W %*% y)

  # Estimação pelo método de Máxima Verossimilhança / Mínimos Quadrados em 2 Estágios (2SLS)
  if (modelo == "sar") {
    # Instrumento para Wy: W * X (variáveis exógenas defasadas)
    WX <- W %*% X[, -1, drop = FALSE]
    H <- cbind(X, WX)
    P_H <- H %*% solve(t(H) %*% H) %*% t(H)
    Wy_hat <- P_H %*% Wy

    Z <- cbind(Wy = Wy_hat, X[, -1, drop = FALSE])
    Z_full <- cbind(Intercepto = 1, Z)

    # Coeficientes
    coefs_est <- solve(t(Z_full) %*% Z_full) %*% t(Z_full) %*% y
    nomes_coefs <- c("(Intercepto)", "Rho (Defasagem \u03c1)", colnames(X)[-1])
    rownames(coefs_est) <- nomes_coefs

    residuos <- as.numeric(y - (cbind(1, Wy, X[, -1, drop = FALSE]) %*% coefs_est))
    s2 <- sum(residuos^2) / (n - length(coefs_est))
    cov_mat <- s2 * solve(t(Z_full) %*% Z_full)
    ep <- sqrt(diag(cov_mat))

    z_vals <- coefs_est / ep
    p_vals <- 2 * (1 - stats::pnorm(abs(z_vals)))
    rho_est <- coefs_est["Rho (Defasagem \u03c1)", 1]
    rho_p   <- p_vals["Rho (Defasagem \u03c1)", 1]

    nome_mod_ext <- "Modelo Espacial Autorregressivo - SAR (Spatial Lag)"
  } else {
    # SEM: Modelo de Erro Espacial
    mod_ols <- stats::lm(formula, data = dados)
    e_ols <- stats::residuals(mod_ols)
    We <- as.numeric(W %*% e_ols)

    reg_lambda <- stats::lm(e_ols ~ We - 1)
    lambda_est <- stats::coef(reg_lambda)[1]

    # Transformação GLS: y* = y - lambda*Wy ; X* = X - lambda*WX
    y_star <- y - lambda_est * Wy
    WX <- W %*% X
    X_star <- X - lambda_est * WX

    mod_gls <- stats::lm(y_star ~ X_star - 1)
    coefs_est <- matrix(stats::coef(mod_gls), ncol = 1)
    nomes_coefs <- colnames(X)
    nomes_coefs[nomes_coefs == "(Intercept)"] <- "(Intercepto)"
    rownames(coefs_est) <- nomes_coefs

    residuos <- as.numeric(y - X %*% coefs_est)
    s2 <- sum(residuos^2) / (n - length(coefs_est))
    cov_mat <- s2 * solve(t(X_star) %*% X_star)
    ep <- sqrt(diag(cov_mat))

    z_vals <- coefs_est / ep
    p_vals <- 2 * (1 - stats::pnorm(abs(z_vals)))

    rho_est <- lambda_est
    rho_p   <- summary(reg_lambda)$coefficients[1, 4]
    nome_mod_ext <- "Modelo de Erro Espacial - SEM (Spatial Error)"
  }

  r2 <- 1 - (sum(residuos^2) / sum((y - mean(y))^2))
  aic_val <- n * log(mean(residuos^2)) + 2 * (length(coefs_est) + 1)
  bic_val <- n * log(mean(residuos^2)) + log(n) * (length(coefs_est) + 1)

  # Teste de Moran nos Resíduos do modelo
  z_res <- residuos - mean(residuos)
  Wz_res <- as.numeric(W %*% z_res)
  moran_res <- (n / sum(W)) * (sum(z_res * Wz_res) / sum(z_res^2))

  # ── Saída no Console ───────────────────────────────────────────────────────
  pad <- function(s, w, align = "left") {
    s <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right") return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }
  sep_l <- function(w) paste0("  ", strrep("\u2500", w))

  asterisco <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    ""
  }

  .print_titulo("REGRESS\u00c3O ESPACIAL")

  cat(sprintf("  Modelo Estimado:   %s\n", nome_mod_ext))
  cat(sprintf("  F\u00f3rmula:           %s\n", deparse(formula)))
  cat(sprintf("  N\u00ba Observa\u00e7\u00f5es:     %d\n", n))
  cat(sprintf("  Matriz de Pesos:   Padroniza\u00e7\u00e3o por %s\n\n", padronizacao))

  .print_topico("COEFICIENTES ESTIMADOS")
  w_c <- c(var = 22, est = 12, ep = 12, z = 10, p = 14, sig = 5)
  w_tot <- sum(w_c) + length(w_c) - 1

  hdr_c <- paste0(
    pad("Vari\u00e1vel",        w_c["var"], "left"),
    pad("Estimativa",     w_c["est"], "center"),
    pad("Erro Pad.",      w_c["ep"],  "center"),
    pad("Estat. Z",       w_c["z"],   "center"),
    pad("p-valor",        w_c["p"],   "center"),
    pad("Sig.",           w_c["sig"], "center")
  )
  cat(sep_l(w_tot), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(sep_l(w_tot), "\n")

  for (i in seq_len(nrow(coefs_est))) {
    nm_v <- rownames(coefs_est)[i]
    pv <- p_vals[i, 1]
    cat("  ",
        pad(nm_v,                   w_c["var"], "left"),
        pad(.fmt_num(coefs_est[i, 1], 4), w_c["est"], "center"),
        pad(.fmt_num(ep[i], 4),           w_c["ep"],  "center"),
        pad(.fmt_num(z_vals[i, 1], 2),    w_c["z"],   "center"),
        pad(.fmt_p(pv),                   w_c["p"],   "center"),
        pad(asterisco(pv),                w_c["sig"], "center"),
        "\n", sep = "")
  }
  cat(sep_l(w_tot), "\n\n")

  cat("  Notas:\n")
  cat("  Signific\u00e2ncia: *** p < 0,001   ** p < 0,01   * p < 0,05\n\n")

  .print_topico("MEDIDAS DE AJUSTE E DIAGN\u00d3STICO ESPACIAL")
  cat(sprintf("  R\u00b2 do Modelo         : %s\n", .fmt_pct(r2, 1)))
  cat(sprintf("  AIC                    : %s  |  BIC: %s\n", .fmt_num(aic_val, 1), .fmt_num(bic_val, 1)))
  cat(sprintf("  Par\u00e2metro Espacial    : %s = %.4f (p = %s)\n",
              if (modelo == "sar") "\u03c1 (Rho)" else "\u03bb (Lambda)", rho_est, .fmt_p(rho_p)))
  cat(sprintf("  I de Moran nos Res\u00edduos: %.4f (Pr\u00f3ximo de zero indica aus\u00eancia de depend\u00eancia)\n\n", moran_res))

  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (rho_p < 0.05) {
    cat(sprintf("  O efeito espacial (%s = %.4f) \u00e9 estatisticamente significativo (p < 0,05).\n",
                if (modelo == "sar") "\u03c1" else "\u03bb", rho_est))
    if (modelo == "sar") {
      cat("  Existe um efeito de transbordamento (spillover): o valor da vari\u00e1vel resposta em uma\n  unidade geogr\u00e1fica \u00e9 influenciado diretamente pelos valores de seus vizinhos.\n")
    } else {
      cat("  A autocorrela\u00e7\u00e3o espacial est\u00e1 presente nos termos de erro (vari\u00e1veis n\u00e3o observadas\n  compartilhadas regionalmente).\n")
    }
  } else {
    cat("  O par\u00e2metro espacial n\u00e3o apresentou signific\u00e2ncia estat\u00edstica (p \u2265 0,05).\n  O modelo linear tradicional (OLS) pode ser suficiente para estes dados.\n")
  }
  cat("\n")
  .print_rodape()

  res <- list(
    modelo_tipo   = modelo,
    coeficientes  = coefs_est,
    erros_padrao  = ep,
    p_valores     = p_vals,
    param_espacial = list(parametro = if (modelo == "sar") "rho" else "lambda", valor = rho_est, p_valor = rho_p),
    r2            = r2,
    aic           = aic_val,
    bic           = bic_val,
    moran_residuos = moran_res,
    residuos      = residuos,
    matriz_pesos  = W,
    formula       = formula
  )
  class(res) <- c("estatR_regressao_espacial", "list")
  invisible(res)
}

# ── DIAGNÓSTICO RESIDUAL DE REGRESSÃO ESPACIAL ───────────────────────────────

#' @keywords internal
analise_residual_espacial <- function(modelo, grafico = TRUE) {
  res <- modelo$residuos
  n   <- length(res)
  W   <- modelo$matriz_pesos

  # Testes de Normalidade
  sw <- tryCatch(stats::shapiro.test(res), error = function(e) list(statistic = NA, p.value = NA))

  # Moran nos resíduos
  z <- res - mean(res)
  Wz <- as.numeric(W %*% z)
  moran_res <- (n / sum(W)) * (sum(z * Wz) / sum(z^2))

  .print_titulo("ANÁLISE RESIDUAL: REGRESSÃO ESPACIAL")
  cat(sprintf("  Modelo:                 %s\n", toupper(modelo$modelo_tipo)))
  cat(sprintf("  Total de Resíduos:      %d\n\n", n))

  .print_topico("DIAGNÓSTICO DE PRESSUPOSTOS")
  cat(sprintf("  Normalidade (Shapiro-Wilk): W = %.4f (p = %s)\n", sw$statistic, .fmt_p(sw$p.value)))
  cat(sprintf("  Autocorrelação Residual (Moran): I = %.4f\n\n", moran_res))

  .print_topico("INTERPRETAÇÃO")
  if (!is.na(sw$p.value) && sw$p.value >= 0.05) {
    cat("  ✓ Os resíduos apresentam distribuição compatível com a Normalidade (p ≥ 0,05).\n")
  } else {
    cat("  ! Evidência de desvio da normalidade nos resíduos (p < 0,05).\n")
  }

  if (abs(moran_res) < 0.1) {
    cat("  ✓ A dependência espacial foi eliminada com sucesso pelo modelo ajustado.\n")
  } else {
    cat("  ! Permanece dependência espacial residual relevante (I de Moran distante de zero).\n")
  }
  cat("\n")
  .print_rodape()

  invisible(list(shapiro = sw, moran_residuos = moran_res, residuos = res))
}
