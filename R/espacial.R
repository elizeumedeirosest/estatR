# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: ESTATÍSTICA ESPACIAL
# Autocorrelação Espacial, Matrizes de Vizinhança, Diagnóstico LM e Regressão Espacial
# ─────────────────────────────────────────────────────────────────────────────

# ── HELPERS INTERNOS DE FORMATAÇÃO E OPERAÇÕES ───────────────────────────────

.esp_fmt_num <- function(x, d = 4) {
  if (is.null(x) || length(x) == 0 || is.na(x) || is.nan(x)) return("-")
  formatC(as.numeric(x), format = "f", digits = d, decimal.mark = ",")
}

.esp_fmt_p <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p) || is.nan(p)) return("-")
  if (p < 0.001) return("< 0,001")
  formatC(as.numeric(p), format = "f", digits = 3, decimal.mark = ",")
}

.esp_asterisco <- function(p) {
  if (is.null(p) || is.na(p) || p >= 0.05) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  "*"
}

.esp_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp / 2)), s, strrep(" ", ceiling(sp / 2))))
  paste0(s, strrep(" ", sp))
}

.esp_sep <- function(w) paste0("  ", strrep("\u2500", w))

.esp_nota_sig <- function() {
  cat("  Notas:\n")
  cat("  Signific\u00e2ncia: *** p < 0,001   ** p < 0,01   * p < 0,05\n\n")
}

# ── 1. MATRIZ DE VIZINHANÇA ESPACIAL (W) ─────────────────────────────────────

#' @title Matriz de Vizinhança Espacial (Pesos W)
#' @description Constrói e padroniza matrizes de pesos espaciais por contiguidade
#'   (Rainha e Torre para malhas poligonais) ou por proximidade euclidiana
#'   (K-vizinhos mais próximos e distância inversa para dados tabulares com coordenadas).
#'
#' @param dados Data frame com coordenadas, matriz de coordenadas (X, Y) ou objeto espacial \code{sf}.
#' @param tipo Tipo de vizinhança: \code{"auto"} (padrão), \code{"rainha"}, \code{"torre"},
#'   \code{"knn"} ou \code{"distancia_inversa"}.
#' @param k Número de vizinhos mais próximos para o método \code{"knn"} (padrão: \code{5}) ou \code{"auto"}.
#' @param coordenadas Nomes das colunas de coordenadas ou matriz de coordenadas quando \code{dados} for data.frame tabular.
#' @param padronizacao Método de padronização: \code{"linha"} (padrão, soma das linhas = 1),
#'   \code{"binaria"} (0 e 1), \code{"global"} (soma total = 1) ou \code{"nenhuma"}.
#'
#' @return Matriz quadrada numérica \eqn{n \times n} de pesos espaciais com atributo de classe \code{estatR_matriz_pesos}.
#' @export
matriz_vizinhanca <- function(dados,
                              tipo = c("auto", "rainha", "torre", "knn", "distancia_inversa"),
                              k = 5,
                              coordenadas = NULL,
                              padronizacao = c("linha", "binaria", "global", "nenhuma")) {

  tipo <- match.arg(tipo)
  padronizacao <- match.arg(padronizacao)

  # 1. Objeto sf (Polígonos / Malha vetorial)
  if (inherits(dados, "sf")) {
    n <- nrow(dados)
    if (tipo %in% c("auto", "rainha", "torre")) {
      is_queen <- (tipo != "torre")
      W <- matrix(0, nrow = n, ncol = n)
      if (requireNamespace("sf", quietly = TRUE)) {
        adj <- tryCatch(sf::st_intersects(dados, dados, sparse = TRUE), error = function(e) NULL)
        if (!is.null(adj)) {
          for (i in seq_len(n)) {
            vzs <- adj[[i]]
            vzs <- vzs[vzs != i]
            if (length(vzs) > 0) W[i, vzs] <- 1
          }
        }
      }
      # Se falhar ou não tiver vizinhos, extrai centroides
      if (all(W == 0)) {
        if (requireNamespace("sf", quietly = TRUE)) {
          coords_pt <- sf::st_coordinates(suppressWarnings(sf::st_centroid(dados)))
          return(matriz_vizinhanca(coords_pt, tipo = "knn", k = k, padronizacao = padronizacao))
        }
      }
    } else {
      # Converte para centroides
      if (requireNamespace("sf", quietly = TRUE)) {
        coords_pt <- sf::st_coordinates(suppressWarnings(sf::st_centroid(dados)))
        return(matriz_vizinhanca(coords_pt, tipo = tipo, k = k, padronizacao = padronizacao))
      }
    }
  } else {
    # 2. Dados Tabulares / Coordenadas
    coords <- NULL
    if (!is.null(coordenadas)) {
      if (is.character(coordenadas) && is.data.frame(dados)) {
        coords <- as.matrix(dados[, coordenadas, drop = FALSE])
      } else {
        coords <- as.matrix(coordenadas)
      }
    } else if (is.data.frame(dados)) {
      cand_cols <- grep("^(x|y|long|lat|lon|latitude|longitude|easting|northing)$", names(dados), ignore.case = TRUE, value = TRUE)
      if (length(cand_cols) >= 2) {
        coords <- as.matrix(dados[, cand_cols[1:2], drop = FALSE])
      } else {
        coords <- matrix(c(seq_len(nrow(dados)), rep(1, nrow(dados))), ncol = 2)
      }
    } else {
      coords <- as.matrix(dados)
    }

    n <- nrow(coords)
    dists <- as.matrix(stats::dist(coords))
    diag(dists) <- Inf

    W <- matrix(0, nrow = n, ncol = n)
    k_ef <- if (identical(k, "auto")) min(5, n - 1) else min(as.integer(k), n - 1)

    if (tipo %in% c("auto", "knn", "rainha", "torre")) {
      for (i in seq_len(n)) {
        vzs <- order(dists[i, ])[seq_len(k_ef)]
        W[i, vzs] <- 1
      }
    } else if (tipo == "distancia_inversa") {
      for (i in seq_len(n)) {
        vzs <- order(dists[i, ])[seq_len(k_ef)]
        d_sub <- dists[i, vzs]
        d_sub[d_sub == 0] <- min(d_sub[d_sub > 0], na.rm = TRUE) / 2
        W[i, vzs] <- 1 / d_sub
      }
    }
  }

  # 3. Padronização da Matriz
  if (padronizacao == "linha") {
    rs <- rowSums(W)
    rs[rs == 0] <- 1
    W <- W / rs
  } else if (padronizacao == "global") {
    st <- sum(W)
    if (st > 0) W <- W / st
  } else if (padronizacao == "binaria") {
    W <- (W > 0) * 1
  }

  attr(W, "tipo") <- tipo
  attr(W, "padronizacao") <- padronizacao
  class(W) <- c("estatR_matriz_pesos", "matrix")
  W
}

#' @noRd
.construir_matriz_pesos <- function(coords, k = 5, padronizacao = "linha") {
  matriz_vizinhanca(coords, tipo = "knn", k = k, padronizacao = padronizacao)
}


# ── 2. AUTOCORRELAÇÃO ESPACIAL (MORAN GLOBAL E LISA) ─────────────────────────

#' @title Análise de Autocorrelação Espacial Global e Local (Moran e LISA)
#' @description Calcula o Índice de Moran Global e os Indicadores Locais de Associação
#'   Espacial (LISA) para identificar padrões de dependência e agrupamentos (Alto-Alto,
#'   Baixo-Baixo e outliers espaciais).
#'
#' @param dados Um \code{data.frame} com as variáveis de interesse.
#' @param variavel Nome da variável numérica a analisar (símbolo ou string).
#' @param coordenadas Matriz ou data frame com coordenadas espaciais (X e Y, ou Long e Lat).
#' @param matriz_pesos Matriz de vizinhança espacial \eqn{W} opcional.
#' @param k_vizinhos Número de vizinhos mais próximos ou \code{"auto"} para seleção ótima.
#' @param n_permutacoes Número de permutações Monte Carlo para o pseudo-p-valor (padrão: 999).
#' @param padronizacao Método de padronização da matriz de pesos (padrão: \code{"linha"}).
#' @param grafico Lógico. Se TRUE (padrão), gera o Diagrama de Dispersão de Moran.
#'
#' @return Objeto invisível da classe \code{estatR_autocorrelacao_espacial}.
#' @export
analise_autocorrelacao <- function(dados, variavel,
                                   coordenadas = NULL,
                                   matriz_pesos = NULL,
                                   k_vizinhos = "auto",
                                   n_permutacoes = 999,
                                   padronizacao = c("linha", "binaria", "global", "nenhuma"),
                                   grafico = TRUE) {

  padronizacao <- match.arg(padronizacao)

  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data.frame.", call. = FALSE)

  var_sub  <- substitute(variavel)
  var_nome <- if (is.symbol(var_sub)) as.character(var_sub) else if (is.character(variavel)) variavel else NULL
  if (is.null(var_nome) || !var_nome %in% names(dados)) {
    stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada nos dados.", var_nome), call. = FALSE)
  }

  y <- as.numeric(dados[[var_nome]])
  idx_val <- !is.na(y)
  y <- y[idx_val]
  n <- length(y)
  if (n < 5) stop("S\u00e3o necess\u00e1rias pelo menos 5 observa\u00e7\u00f5es v\u00e1lidas.", call. = FALSE)

  # Coordenadas
  coords <- NULL
  if (!is.null(coordenadas)) {
    coords <- as.matrix(coordenadas)[idx_val, , drop = FALSE]
  } else {
    cand_cols <- grep("^(x|y|long|lat|lon|latitude|longitude)$", names(dados), ignore.case = TRUE, value = TRUE)
    if (length(cand_cols) >= 2) {
      coords <- as.matrix(dados[idx_val, cand_cols[1:2], drop = FALSE])
    } else {
      coords <- matrix(c(seq_len(n), rep(1, n)), ncol = 2)
    }
  }

  # Seleção automática de k se k_vizinhos == "auto"
  k_sel <- 5
  if (identical(k_vizinhos, "auto") || (is.character(k_vizinhos) && k_vizinhos == "auto")) {
    k_cands <- unique(pmin(pmax(3, c(3, 4, 5, 6, 7, 8)), n - 1))
    melhor_moran <- -1
    for (kc in k_cands) {
      W_tmp <- .construir_matriz_pesos(coords, k = kc, padronizacao = padronizacao)
      z_tmp <- y - mean(y)
      Wz_tmp <- as.numeric(W_tmp %*% z_tmp)
      s_tmp <- sum(z_tmp^2)
      if (s_tmp > 0) {
        m_tmp <- abs((n / sum(W_tmp)) * (sum(z_tmp * Wz_tmp) / s_tmp))
        if (m_tmp > melhor_moran) {
          melhor_moran <- m_tmp
          k_sel <- kc
        }
      }
    }
  } else {
    k_sel <- min(as.integer(k_vizinhos), n - 1)
  }

  if (!is.null(matriz_pesos)) {
    W <- as.matrix(matriz_pesos)[idx_val, idx_val, drop = FALSE]
  } else {
    W <- .construir_matriz_pesos(coords, k = k_sel, padronizacao = padronizacao)
  }

  # Cálculo do Moran Global
  z <- y - mean(y)
  s2 <- sum(z^2)
  S0 <- sum(W)
  Wz <- as.numeric(W %*% z)
  I_global <- (n / S0) * (sum(z * Wz) / s2)
  E_I <- -1 / (n - 1)

  # Permutações Monte Carlo
  I_perm <- numeric(n_permutacoes)
  for (p in seq_len(n_permutacoes)) {
    z_p <- sample(z)
    I_perm[p] <- (n / S0) * (sum(z_p * as.numeric(W %*% z_p)) / s2)
  }
  pseudo_p <- (sum(abs(I_perm - E_I) >= abs(I_global - E_I)) + 1) / (n_permutacoes + 1)
  var_I <- stats::var(I_perm)
  z_score <- if (var_I > 0) (I_global - E_I) / sqrt(var_I) else 0

  # LISA (Moran Local)
  m2 <- s2 / n
  I_local <- (z / m2) * Wz
  p_local <- numeric(n)
  for (i in seq_len(n)) {
    obs_ref <- z[i]
    if (m2 > 0) {
      sim_lisa <- numeric(199)
      for (sp in 1:199) {
        w_row <- sample(W[i, ])
        sim_lisa[sp] <- (obs_ref / m2) * sum(w_row * z)
      }
      p_local[i] <- (sum(abs(sim_lisa) >= abs(I_local[i])) + 1) / 200
    } else {
      p_local[i] <- 1
    }
  }

  # Quadrantes
  z_pad <- z / stats::sd(z)
  Wz_pad <- as.numeric(W %*% z_pad)
  quadrantes <- character(n)
  for (i in seq_len(n)) {
    if (p_local[i] >= 0.05) {
      quadrantes[i] <- "N\u00e3o Significativo"
    } else if (z_pad[i] >= 0 && Wz_pad[i] >= 0) {
      quadrantes[i] <- "Alto-Alto (Q1)"
    } else if (z_pad[i] < 0 && Wz_pad[i] < 0) {
      quadrantes[i] <- "Baixo-Baixo (Q2)"
    } else if (z_pad[i] >= 0 && Wz_pad[i] < 0) {
      quadrantes[i] <- "Alto-Baixo (Q3)"
    } else {
      quadrantes[i] <- "Baixo-Alto (Q4)"
    }
  }

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("AUTOCORRELA\u00c7\u00c3O ESPACIAL: MORAN GLOBAL E LOCAL (LISA)")

  cat(sprintf("  Vari\u00e1vel Analisada : %s\n", var_nome))
  cat(sprintf("  N\u00ba Observa\u00e7\u00f5es   : %d\n", n))
  if (identical(k_vizinhos, "auto") || (is.character(k_vizinhos) && k_vizinhos == "auto")) {
    cat(sprintf("  Vizinhan\u00e7a \u00d3tima  : k = %d vizinhos (selecionada automaticamente)\n", k_sel))
  }
  cat(sprintf("  Matriz de Pesos    : Padroniza\u00e7\u00e3o por %s\n\n", padronizacao))

  .print_topico("\u00cdNDICE DE MORAN GLOBAL")
  w_sep <- 65
  cat(.esp_sep(w_sep), "\n")
  cat(sprintf("  \u00cdndice de Moran (I)    : %s\n", .esp_fmt_num(I_global, 4)))
  cat(sprintf("  Valor Esperado E(I)    : %s\n", .esp_fmt_num(E_I, 4)))
  cat(sprintf("  Estat\u00edstica Z (Score) : %s\n", .esp_fmt_num(z_score, 3)))
  cat(sprintf("  p-valor (Monte Carlo)  : %s %s\n", .esp_fmt_p(pseudo_p), .esp_asterisco(pseudo_p)))
  cat(.esp_sep(w_sep), "\n")
  .esp_nota_sig()

  .print_topico("DISTRIBUI\u00c7\u00c3O DOS AGRUPAMENTOS LOCAIS (LISA)")
  tab_q <- as.data.frame(table(quadrantes), stringsAsFactors = FALSE)
  colnames(tab_q) <- c("Padr\u00e3o", "Frequ\u00eancia")
  tab_q[["Percentual (%)"]] <- sapply(tab_q[["Frequ\u00eancia"]] / n * 100, .esp_fmt_num, d = 1)
  print(tab_q, row.names = FALSE)
  cat("\n")

  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (pseudo_p < 0.05) {
    if (I_global > E_I) {
      cat(sprintf("  [\u2713] Autocorrela\u00e7\u00e3o espacial positiva estatisticamente significativa (I = %s; p %s).\n",
                  .esp_fmt_num(I_global, 4), .esp_fmt_p(pseudo_p)))
      cat(sprintf("  Regi\u00f5es com valores similares de '%s' tendem a estar espacialmente agrupadas\n  (presen\u00e7a de clusters espaciais homog\u00eaneos).\n\n", var_nome))
    } else {
      cat(sprintf("  [\u2713] Autocorrela\u00e7\u00e3o espacial negativa estatisticamente significativa (I = %s; p %s).\n",
                  .esp_fmt_num(I_global, 4), .esp_fmt_p(pseudo_p)))
      cat(sprintf("  Padr\u00e3o espacial disperso ou em mosaico para '%s' (regi\u00f5es com valores altos\n  s\u00e3o predominantemente cercadas por valores baixos).\n\n", var_nome))
    }
  } else {
    cat(sprintf("  [!] N\u00e3o h\u00e1 evid\u00eancia estat\u00edstica de autocorrela\u00e7\u00e3o espacial (I = %s; p = %s).\n",
                .esp_fmt_num(I_global, 4), .esp_fmt_p(pseudo_p)))
    cat(sprintf("  A distribui\u00e7\u00e3o geogr\u00e1fica de '%s' comporta-se de forma espacialmente aleat\u00f3ria.\n\n", var_nome))
  }
  .print_rodape()

  # Gráfico
  p_plot <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot <- data.frame(Z = z_pad, WZ = Wz_pad, Quadrante = quadrantes, stringsAsFactors = FALSE)
    cores_q <- c(
      "Alto-Alto (Q1)"        = "#E63946",
      "Baixo-Baixo (Q2)"      = "#1D3557",
      "Alto-Baixo (Q3)"       = "#F4A261",
      "Baixo-Alto (Q4)"       = "#2A9D8F",
      "N\u00e3o Significativo" = "gray75"
    )

    p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Z, y = WZ, color = Quadrante)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      ggplot2::geom_point(size = 2.8, alpha = 0.85) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#0072CE", linewidth = 1) +
      ggplot2::scale_color_manual(values = cores_q) +
      ggplot2::labs(
        title    = sprintf("Diagrama de Dispers\u00e3o de Moran \u2014 %s", var_nome),
        subtitle = sprintf("I de Moran = %s (p = %s)", .esp_fmt_num(I_global, 4), .esp_fmt_p(pseudo_p)),
        x        = sprintf("%s (Padronizado)", var_nome),
        y        = sprintf("Defasagem Espacial W(%s)", var_nome)
      )

    if (exists("meu_tema", mode = "function")) {
      p_plot <- p_plot + meu_tema(grade = "dupla")
    } else {
      p_plot <- p_plot + ggplot2::theme_minimal()
    }
    print(p_plot)
  }

  res <- list(
    moran_global = list(I = I_global, E_I = E_I, z_score = z_score, p_valor = pseudo_p),
    moran_local  = data.frame(
      id         = seq_len(n),
      valor      = y,
      I_local    = I_local,
      p_valor    = p_local,
      quadrante  = quadrantes,
      stringsAsFactors = FALSE
    ),
    matriz_pesos = W,
    k_selecionado = k_sel,
    grafico      = p_plot
  )
  class(res) <- c("estatR_autocorrelacao_espacial", "list")
  invisible(res)
}


# ── 3. TESTES DOS MULTIPLICADORES DE LAGRANGE (ANSELIN) ──────────────────────

#' @title Testes dos Multiplicadores de Lagrange para Dependência Espacial (Anselin)
#' @description Executa os 4 testes canônicos de Multiplicadores de Lagrange (LM de Anselin):
#'   LM-Defasagem (SAR), LM-Erro (SEM), LM Robusto-Defasagem e LM Robusto-Erro a partir
#'   dos resíduos de Mínimos Quadrados Ordinários (MQO). Inclui árvore de decisão automática
#'   com recomendação explícita do modelo adequado.
#'
#' @param formula Fórmula da regressão linear (ex.: \code{y ~ x1 + x2}).
#' @param dados Data frame com os dados da análise.
#' @param coordenadas Coordenadas espaciais opcionais (X, Y ou Long, Lat).
#' @param matriz_pesos Matriz de vizinhança espacial \eqn{W} opcional.
#' @param k_vizinhos Número de vizinhos mais próximos ou \code{"auto"}.
#' @param padronizacao Método de padronização da matriz: \code{"linha"} (padrão).
#'
#' @return Objeto invisível da classe \code{estatR_lm_espacial}.
#' @export
teste_multiplicadores_lagrange_espacial <- function(formula, dados,
                                                   coordenadas = NULL,
                                                   matriz_pesos = NULL,
                                                   k_vizinhos = "auto",
                                                   padronizacao = c("linha", "binaria", "global", "nenhuma")) {

  padronizacao <- match.arg(padronizacao)

  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data.frame.", call. = FALSE)

  # Ajuste OLS base
  mod_ols <- stats::lm(formula, data = dados)
  e <- stats::residuals(mod_ols)
  y <- stats::model.response(stats::model.frame(formula, data = dados))
  X <- stats::model.matrix(formula, data = dados)
  n <- length(y)
  p_vars <- ncol(X)
  s2 <- sum(e^2) / n

  # Coordenadas e W
  if (!is.null(matriz_pesos)) {
    W <- as.matrix(matriz_pesos)
    k_sel <- "informada"
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

    if (identical(k_vizinhos, "auto") || (is.character(k_vizinhos) && k_vizinhos == "auto")) {
      k_cands <- unique(pmin(pmax(3, c(3, 4, 5, 6, 7, 8)), n - 1))
      melhor_k <- k_cands[1]
      melhor_aic <- Inf
      for (kc in k_cands) {
        W_c <- .construir_matriz_pesos(coords, k = kc, padronizacao = padronizacao)
        Wy_c <- as.numeric(W_c %*% y)
        WX_c <- W_c %*% X[, -1, drop = FALSE]
        H_c <- cbind(X, WX_c)
        P_H <- tryCatch(H_c %*% solve(t(H_c) %*% H_c) %*% t(H_c), error = function(e) NULL)
        if (!is.null(P_H)) {
          Wy_h <- P_H %*% Wy_c
          Z_f <- cbind(Intercepto = 1, Wy = Wy_h, X[, -1, drop = FALSE])
          b_est <- tryCatch(solve(t(Z_f) %*% Z_f) %*% t(Z_f) %*% y, error = function(e) NULL)
          if (!is.null(b_est)) {
            res_c <- as.numeric(y - cbind(1, Wy_c, X[, -1, drop = FALSE]) %*% b_est)
            aic_c <- n * log(mean(res_c^2)) + 2 * (ncol(Z_f) + 1)
            if (aic_c < melhor_aic) {
              melhor_aic <- aic_c
              melhor_k <- kc
            }
          }
        }
      }
      k_sel <- melhor_k
    } else {
      k_sel <- min(as.integer(k_vizinhos), n - 1)
    }
    W <- .construir_matriz_pesos(coords, k = k_sel, padronizacao = padronizacao)
  }

  # Defasagens
  We <- as.numeric(W %*% e)
  Wy <- as.numeric(W %*% y)
  y_hat <- stats::fitted(mod_ols)
  Wy_hat <- as.numeric(W %*% y_hat)

  # Traço e Matriz de Projeção
  W2 <- W %*% W
  tr_term <- sum(diag(W2)) + sum(diag(t(W) %*% W))
  T_val <- tr_term

  # Projeção M = I - X(X'X)^(-1)X'
  XtX_inv <- solve(t(X) %*% X)
  M_Wyhat <- Wy_hat - X %*% (XtX_inv %*% (t(X) %*% Wy_hat))
  D_val <- (as.numeric(t(M_Wyhat) %*% M_Wyhat) / s2) + T_val

  # 1. LM-Erro
  LM_erro <- (sum(e * We) / s2)^2 / T_val
  p_erro  <- stats::pchisq(LM_erro, df = 1, lower.tail = FALSE)

  # 2. LM-Defasagem (Lag)
  LM_lag  <- (sum(e * Wy) / s2)^2 / D_val
  p_lag   <- stats::pchisq(LM_lag, df = 1, lower.tail = FALSE)

  # 3. Robust LM-Erro
  num_r_erro <- sum(e * We) / s2 - (T_val / D_val) * (sum(e * Wy) / s2)
  den_r_erro <- T_val * (1 - (T_val / D_val))
  RLM_erro   <- if (den_r_erro > 0) (num_r_erro^2) / den_r_erro else 0
  p_r_erro   <- stats::pchisq(RLM_erro, df = 1, lower.tail = FALSE)

  # 4. Robust LM-Defasagem (Lag)
  den_r_lag <- D_val - T_val
  num_r_lag <- (sum(e * Wy) / s2) - (sum(e * We) / s2)
  RLM_lag   <- if (den_r_lag > 0) (num_r_lag^2) / den_r_lag else 0
  p_r_lag   <- stats::pchisq(RLM_lag, df = 1, lower.tail = FALSE)

  # Decisão de Anselin
  recomendacao <- ""
  detalhe_rec  <- ""

  if (p_erro >= 0.05 && p_lag >= 0.05) {
    recomendacao <- "Modelo Linear Cl\u00e1ssico (MQO / OLS)"
    detalhe_rec  <- "Nenhum teste LM detectou depend\u00eancia espacial estatisticamente significativa (p \u2265 0,05).\n    O modelo linear cl\u00e1ssico sem termos espaciais \u00e9 o mais parcimonioso e adequado."
  } else if (p_lag < 0.05 && p_erro >= 0.05) {
    recomendacao <- "Modelo de Defasagem Espacial (SAR / Spatial Lag)"
    detalhe_rec  <- "Apenas o teste LM-Defasagem foi significativo (p < 0,05).\n    A depend\u00eancia espacial ocorre na vari\u00e1vel resposta y (efeito de transbordamento)."
  } else if (p_erro < 0.05 && p_lag >= 0.05) {
    recomendacao <- "Modelo de Erro Espacial (SEM / Spatial Error)"
    detalhe_rec  <- "Apenas o teste LM-Erro foi significativo (p < 0,05).\n    A depend\u00eancia espacial reside nos res\u00edduos (fatores n\u00e3o observados espacialmente correlacionados)."
  } else {
    # Ambos foram significativos: avalia os robustos
    if (p_r_lag < 0.05 && p_r_erro >= 0.05) {
      recomendacao <- "Modelo de Defasagem Espacial (SAR / Spatial Lag)"
      detalhe_rec  <- "Ambos os testes simples foram significativos, por\u00e9m apenas o LM Robusto-Defasagem\n    manteve signific\u00e2ncia estat\u00edstica (p < 0,05)."
    } else if (p_r_erro < 0.05 && p_r_lag >= 0.05) {
      recomendacao <- "Modelo de Erro Espacial (SEM / Spatial Error)"
      detalhe_rec  <- "Ambos os testes simples foram significativos, por\u00e9m apenas o LM Robusto-Erro\n    manteve signific\u00e2ncia estat\u00edstica (p < 0,05)."
    } else if (p_r_lag < 0.05 && p_r_erro < 0.05) {
      if (RLM_lag > RLM_erro) {
        recomendacao <- "Modelo de Defasagem Espacial (SAR) ou Modelo de Durbin Espacial (MDE)"
        detalhe_rec  <- "Ambos os testes robustos foram significativos, com maior magnitude em Defasagem.\n    Recomenda-se o Modelo SAR ou Modelo de Durbin Espacial (MDE)."
      } else {
        recomendacao <- "Modelo de Erro Espacial (SEM) ou Modelo de Durbin Espacial (MDE)"
        detalhe_rec  <- "Ambos os testes robustos foram significativos, com maior magnitude em Erro.\n    Recomenda-se o Modelo SEM ou Modelo de Durbin Espacial (MDE)."
      }
    } else {
      if (LM_lag > LM_erro) {
        recomendacao <- "Modelo de Defasagem Espacial (SAR / Spatial Lag)"
        detalhe_rec  <- "Testes robustos inconclusivos; pela magnitude de LM-Defasagem recomenda-se SAR."
      } else {
        recomendacao <- "Modelo de Erro Espacial (SEM / Spatial Error)"
        detalhe_rec  <- "Testes robustos inconclusivos; pela magnitude de LM-Erro recomenda-se SEM."
      }
    }
  }

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("DIAGN\u00d3STICO ESPACIAL: TESTES DO MULTIPLICADOR DE LAGRANGE (ANSELIN)")

  cat(sprintf("  F\u00f3rmula Analisada : %s\n", deparse(formula)))
  cat(sprintf("  N\u00ba Observa\u00e7\u00f5es    : %d\n", n))
  if (k_sel != "informada") {
    cat(sprintf("  Vizinhan\u00e7a \u00d3tima   : k = %s vizinhos (selecionada automaticamente)\n", k_sel))
  }
  cat(sprintf("  Matriz de Pesos     : Padroniza\u00e7\u00e3o por %s\n\n", padronizacao))

  .print_topico("RESULTADOS DOS TESTES LM")

  w_tb <- c(teste = 32, chi = 16, gl = 12, p = 16)
  w_tot <- sum(w_tb)

  hdr <- paste0(
    .esp_pad("Teste de Diagn\u00f3stico", w_tb["teste"], "left"),
    .esp_pad("Estat\u00edstica (\u03c7\u00b2)",  w_tb["chi"],   "center"),
    .esp_pad("Graus Lib.",      w_tb["gl"],    "center"),
    .esp_pad("p-valor",         w_tb["p"],     "center")
  )

  cat(.esp_sep(w_tot), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.esp_sep(w_tot), "\n")

  linhas_res <- list(
    list("LM - Defasagem Espacial (Lag)", LM_lag, 1, p_lag),
    list("LM - Erro Espacial (Error)",       LM_erro, 1, p_erro),
    list("LM Robusto - Defasagem Espacial",  RLM_lag, 1, p_r_lag),
    list("LM Robusto - Erro Espacial",      RLM_erro, 1, p_r_erro)
  )

  for (lr in linhas_res) {
    p_str <- paste0(.esp_fmt_p(lr[[4]]), " ", .esp_asterisco(lr[[4]]))
    cat("  ",
        .esp_pad(lr[[1]], w_tb["teste"], "left"),
        .esp_pad(.esp_fmt_num(lr[[2]], 3), w_tb["chi"], "center"),
        .esp_pad(lr[[3]], w_tb["gl"], "center"),
        .esp_pad(p_str, w_tb["p"], "center"),
        "\n", sep = "")
  }
  cat(.esp_sep(w_tot), "\n")
  .esp_nota_sig()

  .print_topico("DIAGN\u00d3STICO E RECOMENDA\u00c7\u00c3O METODOL\u00d3GICA")
  cat(sprintf("  [\u2713] Modelo Recomendado: %s\n", recomendacao))
  cat(sprintf("    %s\n\n", detalhe_rec))

  .print_rodape()

  df_res <- data.frame(
    Teste           = c("LM_Defasagem", "LM_Erro", "Robust_LM_Defasagem", "Robust_LM_Erro"),
    Estatistica_Chi = c(LM_lag, LM_erro, RLM_lag, RLM_erro),
    GL              = c(1, 1, 1, 1),
    p_valor         = c(p_lag, p_erro, p_r_lag, p_r_erro),
    stringsAsFactors = FALSE
  )

  res <- list(
    tabela        = df_res,
    recomendacao  = recomendacao,
    detalhes      = detalhe_rec,
    matriz_pesos  = W,
    k_selecionado = k_sel,
    modelo_ols    = mod_ols
  )
  class(res) <- c("estatR_lm_espacial", "list")
  invisible(res)
}

#' @export
teste_mult_lagrange_espacial <- teste_multiplicadores_lagrange_espacial

#' @export
teste_lm_espacial <- teste_multiplicadores_lagrange_espacial


# ── 4. REGRESSÃO ESPACIAL (DEFASAGEM, ERRO E DURBIN) ─────────────────────────

#' @title Regressão Espacial (Defasagem, Erro e Modelo de Durbin Espacial)
#' @description Ajusta modelos de regressão espacial para dados de área: Modelo de Defasagem Espacial
#'   (\code{"defasagem"} / SAR), Modelo de Erro Espacial (\code{"erro"} / SEM) e Modelo de Durbin Espacial
#'   (\code{"durbin"} / MDE) com seleção de vizinhança automática por AIC e matriz de pesos padronizada por linha.
#'
#' @param formula Fórmula da regressão (ex.: \code{y ~ x1 + x2}).
#' @param dados Data frame contendo as variáveis.
#' @param coordenadas Coordenadas espaciais opcionais (X, Y ou Long, Lat).
#' @param matriz_pesos Matriz de pesos espaciais \eqn{W} opcional.
#' @param modelo Tipo de modelo em português: \code{"defasagem"} (SAR, padrão),
#'   \code{"erro"} (SEM) ou \code{"durbin"} (MDE / SDM). Também aceita \code{"sar"}, \code{"sem"}, \code{"sdm"}, \code{"mde"}.
#' @param k_vizinhos Número de vizinhos mais próximos ou \code{"auto"} para otimização automática por AIC.
#' @param padronizacao Método de padronização da matriz de pesos (padrão: \code{"linha"}).
#' @param ... Argumentos adicionais.
#'
#' @return Objeto invisível da classe \code{estatR_regressao_espacial}.
#' @export
regressao_espacial <- function(formula, dados,
                               coordenadas = NULL,
                               matriz_pesos = NULL,
                               modelo = c("defasagem", "erro", "durbin", "sar", "sem", "sdm", "mde"),
                               k_vizinhos = "auto",
                               padronizacao = c("linha", "binaria", "global", "nenhuma"),
                               ...) {

  modelo_in <- match.arg(modelo)
  padronizacao <- match.arg(padronizacao)

  # Normalização de nomes para português
  modelo_tipo <- switch(modelo_in,
    sar = "defasagem",
    sem = "erro",
    sdm = "durbin",
    mde = "durbin",
    modelo_in
  )

  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data.frame.", call. = FALSE)

  mf <- stats::model.frame(formula, data = dados)
  y  <- stats::model.response(mf)
  X  <- stats::model.matrix(formula, data = mf)
  n  <- length(y)

  # 1. Matriz de pesos W e seleção ótima de k
  k_sel <- 5
  if (!is.null(matriz_pesos)) {
    W <- as.matrix(matriz_pesos)
    k_sel <- "informada"
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

    if (identical(k_vizinhos, "auto") || (is.character(k_vizinhos) && k_vizinhos == "auto")) {
      k_cands <- unique(pmin(pmax(3, c(3, 4, 5, 6, 7, 8)), n - 1))
      melhor_aic <- Inf
      melhor_k <- k_cands[1]
      for (kc in k_cands) {
        W_c <- .construir_matriz_pesos(coords, k = kc, padronizacao = padronizacao)
        Wy_c <- as.numeric(W_c %*% y)
        WX_c <- W_c %*% X[, -1, drop = FALSE]
        H_c  <- cbind(X, WX_c)
        P_H  <- tryCatch(H_c %*% solve(t(H_c) %*% H_c) %*% t(H_c), error = function(e) NULL)
        if (!is.null(P_H)) {
          Wy_h <- P_H %*% Wy_c
          Z_f  <- cbind(Intercepto = 1, Wy = Wy_h, X[, -1, drop = FALSE])
          b_est <- tryCatch(solve(t(Z_f) %*% Z_f) %*% t(Z_f) %*% y, error = function(e) NULL)
          if (!is.null(b_est)) {
            res_c <- as.numeric(y - cbind(1, Wy_c, X[, -1, drop = FALSE]) %*% b_est)
            aic_c <- n * log(mean(res_c^2)) + 2 * (ncol(Z_f) + 1)
            if (aic_c < melhor_aic) {
              melhor_aic <- aic_c
              melhor_k   <- kc
            }
          }
        }
      }
      k_sel <- melhor_k
    } else {
      k_sel <- min(as.integer(k_vizinhos), n - 1)
    }
    W <- .construir_matriz_pesos(coords, k = k_sel, padronizacao = padronizacao)
  }

  Wy <- as.numeric(W %*% y)

  # 2. Estimação dos Modelos
  if (modelo_tipo == "defasagem") {
    # SAR: y = rho*Wy + X*beta + e
    WX  <- W %*% X[, -1, drop = FALSE]
    H   <- cbind(X, WX)
    P_H <- H %*% solve(t(H) %*% H) %*% t(H)
    Wy_hat <- P_H %*% Wy

    Z_full <- cbind(Intercepto = 1, Wy = Wy_hat, X[, -1, drop = FALSE])
    coefs_est <- solve(t(Z_full) %*% Z_full) %*% t(Z_full) %*% y
    nomes_coefs <- c("(Intercepto)", "Defasagem Espacial (\u03c1)", colnames(X)[-1])
    rownames(coefs_est) <- nomes_coefs

    residuos <- as.numeric(y - (cbind(1, Wy, X[, -1, drop = FALSE]) %*% coefs_est))
    s2 <- sum(residuos^2) / (n - length(coefs_est))
    cov_mat <- s2 * solve(t(Z_full) %*% Z_full)
    ep <- sqrt(diag(cov_mat))

    z_vals <- coefs_est / ep
    p_vals <- 2 * (1 - stats::pnorm(abs(z_vals)))
    param_val <- coefs_est["Defasagem Espacial (\u03c1)", 1]
    param_p   <- p_vals["Defasagem Espacial (\u03c1)", 1]
    param_nome <- "\u03c1 (Rho)"
    nome_mod_ext <- "Modelo de Defasagem Espacial (SAR / Spatial Lag)"

  } else if (modelo_tipo == "durbin") {
    # MDE: y = rho*Wy + X*beta + WX*gamma + e
    WX <- W %*% X[, -1, drop = FALSE]
    W2X <- W %*% WX
    H   <- cbind(X, WX, W2X)
    P_H <- H %*% solve(t(H) %*% H) %*% t(H)
    Wy_hat <- P_H %*% Wy

    nomes_wx <- paste0("W \u00d7 ", colnames(X)[-1])
    Z_full <- cbind(Intercepto = 1, Wy = Wy_hat, X[, -1, drop = FALSE], WX)
    coefs_est <- solve(t(Z_full) %*% Z_full) %*% t(Z_full) %*% y
    nomes_coefs <- c("(Intercepto)", "Defasagem Espacial (\u03c1)", colnames(X)[-1], nomes_wx)
    rownames(coefs_est) <- nomes_coefs

    residuos <- as.numeric(y - (cbind(1, Wy, X[, -1, drop = FALSE], WX) %*% coefs_est))
    s2 <- sum(residuos^2) / (n - length(coefs_est))
    cov_mat <- s2 * solve(t(Z_full) %*% Z_full)
    ep <- sqrt(diag(cov_mat))

    z_vals <- coefs_est / ep
    p_vals <- 2 * (1 - stats::pnorm(abs(z_vals)))
    param_val <- coefs_est["Defasagem Espacial (\u03c1)", 1]
    param_p   <- p_vals["Defasagem Espacial (\u03c1)", 1]
    param_nome <- "\u03c1 (Rho)"
    nome_mod_ext <- "Modelo de Durbin Espacial (MDE / Spatial Durbin)"

  } else {
    # SEM: y = X*beta + u ; u = lambda*Wu + e
    mod_ols <- stats::lm(formula, data = dados)
    e_ols <- stats::residuals(mod_ols)
    We <- as.numeric(W %*% e_ols)

    reg_lambda <- stats::lm(e_ols ~ We - 1)
    lambda_est <- stats::coef(reg_lambda)[1]

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
    param_val <- lambda_est
    param_p   <- summary(reg_lambda)$coefficients[1, 4]
    param_nome <- "\u03bb (Lambda)"
    nome_mod_ext <- "Modelo de Erro Espacial (SEM / Spatial Error)"
  }

  r2 <- 1 - (sum(residuos^2) / sum((y - mean(y))^2))
  aic_val <- n * log(mean(residuos^2)) + 2 * (length(coefs_est) + 1)
  bic_val <- n * log(mean(residuos^2)) + log(n) * (length(coefs_est) + 1)
  loglik_val <- -0.5 * n * (log(2 * pi) + log(mean(residuos^2)) + 1)

  # Teste de Moran nos Resíduos
  z_res <- residuos - mean(residuos)
  Wz_res <- as.numeric(W %*% z_res)
  moran_res <- (n / sum(W)) * (sum(z_res * Wz_res) / sum(z_res^2))

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("REGRESS\u00c3O ESPACIAL")

  cat(sprintf("  Modelo Estimado : %s\n", nome_mod_ext))
  cat(sprintf("  F\u00f3rmula         : %s\n", deparse(formula)))
  cat(sprintf("  N\u00ba Observa\u00e7\u00f5es   : %d\n", n))
  if (k_sel != "informada") {
    cat(sprintf("  Vizinhan\u00e7a \u00d3tima : k = %s vizinhos (selecionada automaticamente por AIC)\n", k_sel))
  }
  cat(sprintf("  Matriz de Pesos : Padroniza\u00e7\u00e3o por %s\n\n", padronizacao))

  .print_topico("COEFICIENTES ESTIMADOS")

  # Tabela sem coluna Sig., asteriscos colados ao p-valor
  w_c <- c(var = 26, est = 14, ep = 14, z = 12, p = 16)
  w_tot <- sum(w_c)

  hdr_c <- paste0(
    .esp_pad("Par\u00e2metro",       w_c["var"], "left"),
    .esp_pad("Estimativa",     w_c["est"], "center"),
    .esp_pad("Erro Padr\u00e3o", w_c["ep"],  "center"),
    .esp_pad("Estat\u00edstica Z", w_c["z"], "center"),
    .esp_pad("p-valor",        w_c["p"],   "center")
  )
  cat(.esp_sep(w_tot), "\n")
  cat("  ", hdr_c, "\n", sep = "")
  cat(.esp_sep(w_tot), "\n")

  for (i in seq_len(nrow(coefs_est))) {
    nm_v <- rownames(coefs_est)[i]
    pv   <- p_vals[i, 1]
    p_str <- paste0(.esp_fmt_p(pv), " ", .esp_asterisco(pv))
    cat("  ",
        .esp_pad(nm_v,                   w_c["var"], "left"),
        .esp_pad(.esp_fmt_num(coefs_est[i, 1], 4), w_c["est"], "center"),
        .esp_pad(.esp_fmt_num(ep[i], 4),           w_c["ep"],  "center"),
        .esp_pad(.esp_fmt_num(z_vals[i, 1], 2),    w_c["z"],   "center"),
        .esp_pad(p_str,                            w_c["p"],   "center"),
        "\n", sep = "")
  }
  cat(.esp_sep(w_tot), "\n")
  .esp_nota_sig()

  .print_topico("MEDIDAS DE AJUSTE E DIAGN\u00d3STICO ESPACIAL")
  cat(sprintf("  Pseudo-R\u00b2 do Modelo : %s\n", .esp_fmt_num(r2, 3)))
  cat(sprintf("  AIC                    : %s   |   BIC: %s\n", .esp_fmt_num(aic_val, 1), .esp_fmt_num(bic_val, 1)))
  cat(sprintf("  Par\u00e2metro Espacial    : %s = %s (p = %s)\n",
              param_nome, .esp_fmt_num(param_val, 4), .esp_fmt_p(param_p)))
  cat(sprintf("  I de Moran nos Res\u00edduos: %s (Pr\u00f3ximo de zero indica aus\u00eancia de depend\u00eancia residual)\n\n", .esp_fmt_num(moran_res, 4)))

  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (param_p < 0.05) {
    cat(sprintf("  [\u2713] O par\u00e2metro de depend\u00eancia espacial (%s = %s) \u00e9 estatisticamente significativo (p < 0,05).\n",
                param_nome, .esp_fmt_num(param_val, 4)))
    if (modelo_tipo %in% c("defasagem", "durbin")) {
      cat("  Existe um efeito relevante de transbordamento (spillover): a vari\u00e1vel resposta em uma\n  regi\u00e3o \u00e9 influenciada de forma direta pelo comportamento de seus vizinhos.\n")
      cat("  Utilize impactos_espaciais() para decompor os efeitos Diretos e Indiretos (Spillovers).\n\n")
    } else {
      cat("  A autocorrela\u00e7\u00e3o espacial manifesta-se nos termos de erro (vari\u00e1veis n\u00e3o observadas\n  compartilhadas geograficamente entre regi\u00f5es cont\u00edguas).\n\n")
    }
  } else {
    cat(sprintf("  [!] O par\u00e2metro espacial (%s = %s) n\u00e3o apresentou signific\u00e2ncia estat\u00edstica (p \u2265 0,05).\n",
                param_nome, .esp_fmt_num(param_val, 4)))
    cat("  O modelo linear tradicional (OLS/MQO) pode ser suficiente para estes dados.\n\n")
  }

  .print_rodape()

  res <- list(
    modelo_tipo   = modelo_tipo,
    coeficientes  = coefs_est,
    erros_padrao  = ep,
    p_valores     = p_vals,
    param_espacial = list(parametro = if (modelo_tipo %in% c("defasagem", "durbin")) "rho" else "lambda", valor = param_val, p_valor = param_p),
    r2            = r2,
    r2_pseudo     = r2,
    aic           = aic_val,
    bic           = bic_val,
    loglik        = loglik_val,
    moran_residuos = moran_res,
    residuos      = residuos,
    valores_ajustados = as.numeric(y - residuos),
    matriz_pesos  = W,
    k_selecionado = k_sel,
    formula       = formula,
    dados         = dados
  )
  class(res) <- c("estatR_regressao_espacial", "list")
  invisible(res)
}


# ── 5. DECOMPOSIÇÃO DE IMPACTOS ESPACIAIS (LESAGE & PACE) ────────────────────

#' @title Decomposição de Impactos Espaciais (Efeitos Direto, Indireto e Total)
#' @description Calcula a partição clássica de impactos de LeSage & Pace (2009)
#'   para modelos com defasagem espacial (SAR e MDE): Impacto Direto (efeito local),
#'   Impacto Indireto (efeito de transbordamento / spillover sobre os vizinhos) e Impacto Total.
#'
#' @param modelo Objeto retornado por \code{regressao_espacial()} ajustado sob modelo \code{"defasagem"} ou \code{"durbin"}.
#' @return Objeto invisível da classe \code{estatR_impactos_espaciais}.
#' @export
impactos_espaciais <- function(modelo) {
  if (!inherits(modelo, "estatR_regressao_espacial")) {
    stop("O argumento 'modelo' deve ser um objeto retornado por regressao_espacial().", call. = FALSE)
  }

  if (!modelo$modelo_tipo %in% c("defasagem", "durbin")) {
    stop("A decomposi\u00e7\u00e3o de impactos \u00e9 v\u00e1lida apenas para modelos com defasagem espacial ('defasagem' ou 'durbin').", call. = FALSE)
  }

  W <- modelo$matriz_pesos
  n <- nrow(W)
  rho <- as.numeric(modelo$param_espacial$valor)

  # Inversa de Leontief espacial: S = (I - rho*W)^(-1)
  I_n <- diag(n)
  A_inv <- tryCatch(solve(I_n - rho * W), error = function(e) {
    # Aproximação de série de potências se singular
    I_n + rho * W + (rho^2) * (W %*% W)
  })

  coefs <- modelo$coeficientes
  nms_c <- rownames(coefs)

  # Identifica variáveis explicativas
  idx_x <- which(!nms_c %in% c("(Intercepto)", "Defasagem Espacial (\u03c1)") & !grepl("^W \u00d7 ", nms_c))
  vars_x <- nms_c[idx_x]

  df_impactos <- data.frame(
    Variavel         = vars_x,
    Impacto_Direto   = numeric(length(vars_x)),
    Impacto_Indireto = numeric(length(vars_x)),
    Impacto_Total    = numeric(length(vars_x)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(vars_x)) {
    vx <- vars_x[i]
    beta_k <- as.numeric(coefs[vx, 1])

    # Termo defasado se for Durbin
    nm_wx <- paste0("W \u00d7 ", vx)
    gamma_k <- if (nm_wx %in% nms_c) as.numeric(coefs[nm_wx, 1]) else 0

    # Matriz de impacto do regressor k: S_k = A_inv %*% (beta_k*I + gamma_k*W)
    S_k <- A_inv %*% (beta_k * I_n + gamma_k * W)

    dir_val <- mean(diag(S_k))
    tot_val <- sum(S_k) / n
    ind_val <- tot_val - dir_val

    df_impactos$Impacto_Direto[i]   <- dir_val
    df_impactos$Impacto_Indireto[i] <- ind_val
    df_impactos$Impacto_Total[i]    <- tot_val
  }

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("DECOMPOSI\u00c7\u00c3O DE IMPACTOS ESPACIAIS (LESAGE & PACE)")

  cat(sprintf("  Modelo Base     : %s\n", toupper(modelo$modelo_tipo)))
  cat(sprintf("  Par\u00e2metro \u03c1      : %s\n", .esp_fmt_num(rho, 4)))
  cat("  Metodologia     : Partic\u00e3o de Efeitos em Matriz de Multiplicador Espacial\n\n")

  .print_topico("TABELA DE IMPACTOS (DIRETO, INDIRETO E TOTAL)")

  w_tb <- c(var = 22, dir = 18, ind = 24, tot = 18)
  w_tot <- sum(w_tb)

  hdr <- paste0(
    .esp_pad("Vari\u00e1vel",                 w_tb["var"], "left"),
    .esp_pad("Impacto Direto",           w_tb["dir"], "center"),
    .esp_pad("Indireto (Spillover)",     w_tb["ind"], "center"),
    .esp_pad("Impacto Total",            w_tb["tot"], "center")
  )

  cat(.esp_sep(w_tot), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.esp_sep(w_tot), "\n")

  for (i in seq_len(nrow(df_impactos))) {
    cat("  ",
        .esp_pad(df_impactos$Variavel[i],                    w_tb["var"], "left"),
        .esp_pad(.esp_fmt_num(df_impactos$Impacto_Direto[i], 4),   w_tb["dir"], "center"),
        .esp_pad(.esp_fmt_num(df_impactos$Impacto_Indireto[i], 4), w_tb["ind"], "center"),
        .esp_pad(.esp_fmt_num(df_impactos$Impacto_Total[i], 4),    w_tb["tot"], "center"),
        "\n", sep = "")
  }
  cat(.esp_sep(w_tot), "\n\n")

  .print_topico("INTERPRETA\u00c7\u00c3O DOS EFEITOS DE TRANSBORDAMENTO")
  for (i in seq_len(nrow(df_impactos))) {
    vx <- df_impactos$Variavel[i]
    d_val <- df_impactos$Impacto_Direto[i]
    in_val <- df_impactos$Impacto_Indireto[i]
    pct_spill <- if (abs(df_impactos$Impacto_Total[i]) > 0) (in_val / df_impactos$Impacto_Total[i]) * 100 else 0

    cat(sprintf("  [\u2713] %s: Para cada aumento unit\u00e1rio nesta vari\u00e1vel, o impacto m\u00e9dio local \u00e9 de %s,\n    enquanto o efeito de transbordamento gerado sobre as regi\u00f5es vizinhas \u00e9 de %s\n    (representando %.1f%% do impacto total consolidado).\n\n",
                vx, .esp_fmt_num(d_val, 4), .esp_fmt_num(in_val, 4), pct_spill))
  }

  .print_rodape()

  res <- list(
    impactos = df_impactos,
    modelo   = modelo
  )
  class(res) <- c("estatR_impactos_espaciais", "list")
  invisible(res)
}


# ── 6. DIAGNÓSTICO RESIDUAL INTERNO (POLIMORFISMO PARA analise_residual) ─────

#' @noRd
.analise_residual_espacial <- function(modelo, grafico = TRUE) {
  res <- modelo$residuos
  n   <- length(res)
  W   <- modelo$matriz_pesos

  # Normalidade
  sw <- tryCatch(stats::shapiro.test(res), error = function(e) list(statistic = NA, p.value = NA))

  # Moran nos resíduos
  z <- res - mean(res)
  Wz <- as.numeric(W %*% z)
  moran_res <- if (sum(z^2) > 0) (n / sum(W)) * (sum(z * Wz) / sum(z^2)) else 0

  .print_titulo("DIAGN\u00d3STICO RESIDUAL: REGRESS\u00c3O ESPACIAL")
  cat(sprintf("  Modelo Analisado:       %s\n", toupper(modelo$modelo_tipo)))
  cat(sprintf("  Total de Res\u00edduos:      %d\n\n", n))

  .print_topico("TESTES DE PRESSUPOSTOS DOS RES\u00cdDUOS")
  w_sep <- 65
  cat(.esp_sep(w_sep), "\n")
  cat(sprintf("  Normalidade (Shapiro-Wilk) : W = %s  |  p-valor = %s %s\n",
              .esp_fmt_num(sw$statistic, 3), .esp_fmt_p(sw$p.value), .esp_asterisco(sw$p.value)))
  cat(sprintf("  I de Moran nos Res\u00edduos     : I = %s  (Pr\u00f3ximo de zero = sem depend\u00eancia)\n",
              .esp_fmt_num(moran_res, 4)))
  cat(.esp_sep(w_sep), "\n")
  .esp_nota_sig()

  .print_topico("INTERPRETA\u00c7\u00c3O DO DIAGN\u00d3STICO")
  if (!is.na(sw$p.value) && sw$p.value >= 0.05) {
    cat("  [\u2713] Os res\u00edduos seguem distribui\u00e7\u00e3o normal (p \u2265 0,05).\n")
  } else {
    cat("  [!] Os res\u00edduos apresentaram desvio da normalidade (p < 0,05).\n")
  }

  if (abs(moran_res) < 0.10) {
    cat("  [\u2713] O modelo eliminou a depend\u00eancia espacial dos res\u00edduos (I de Moran residual pr\u00f3ximo de zero).\n\n")
  } else {
    cat(sprintf("  [!] Ainda h\u00e1 autocorrela\u00e7\u00e3o espacial remanescente nos res\u00edduos (I = %s).\n    Considere testar o Modelo de Durbin Espacial (MDE).\n\n", .esp_fmt_num(moran_res, 4)))
  }

  .print_rodape()

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_p <- data.frame(
      Indice = seq_len(n),
      Residuo = res,
      Ajustado = modelo$valores_ajustados,
      stringsAsFactors = FALSE
    )
    p1 <- ggplot2::ggplot(df_p, ggplot2::aes(x = Ajustado, y = Residuo)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_point(color = "#0072CE", size = 2.4, alpha = 0.8) +
      ggplot2::labs(title = "Res\u00edduos vs Valores Ajustados", x = "Valores Ajustados", y = "Res\u00edduos")

    p2 <- ggplot2::ggplot(df_p, ggplot2::aes(sample = Residuo)) +
      ggplot2::stat_qq(color = "#1D3557", size = 2.2) +
      ggplot2::stat_qq_line(color = "#E63946", linewidth = 0.8) +
      ggplot2::labs(title = "Q-Q Plot da Normalidade", x = "Quantis Te\u00f3ricos", y = "Quantis Amostrais")

    if (exists("meu_tema", mode = "function")) {
      p1 <- p1 + meu_tema(grade = "dupla")
      p2 <- p2 + meu_tema(grade = "dupla")
    }

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    } else {
      print(p1)
    }
  }

  invisible(list(shapiro = sw, moran_residuos = moran_res, residuos = res))
}

#' @noRd
.metricas_espacial <- function(modelo) {
  res <- modelo$residuos
  y <- modelo$valores_ajustados + res
  rmse <- sqrt(mean(res^2, na.rm = TRUE))
  mae  <- mean(abs(res), na.rm = TRUE)
  r2   <- if (!is.null(modelo$r2_pseudo)) modelo$r2_pseudo else NA_real_
  aic_val <- if (!is.null(modelo$aic)) modelo$aic else NA_real_
  bic_val <- if (!is.null(modelo$bic)) modelo$bic else NA_real_
  loglik_val <- if (!is.null(modelo$loglik)) modelo$loglik else NA_real_
  moran_res  <- if (!is.null(modelo$moran_residuos)) modelo$moran_residuos else NA_real_
  p_esp      <- if (!is.null(modelo$param_espacial$valor)) modelo$param_espacial$valor else NA_real_

  f2 <- function(x) if (is.na(x)) "-" else formatC(x, format = "f", digits = 2, decimal.mark = ",")
  f4 <- function(x) if (is.na(x)) "-" else formatC(x, format = "f", digits = 4, decimal.mark = ",")

  lw <- 24
  vw <- 14
  sep_m <- paste0("  ", strrep("─", lw + vw + 1), "\n")

  .print_titulo("MÉTRICAS DO MODELO ESPACIAL")
  cat(sprintf("  %-*s %*s\n", lw, "Métrica", vw, "Valor"))
  cat(sep_m)
  cat(sprintf("  %-*s %*s\n", lw, "Pseudo-R² (%)",        vw, f2(r2 * 100)))
  cat(sprintf("  %-*s %*s\n", lw, "RMSE",                 vw, f4(rmse)))
  cat(sprintf("  %-*s %*s\n", lw, "MAE",                  vw, f4(mae)))
  cat(sprintf("  %-*s %*s\n", lw, "AIC",                  vw, f2(aic_val)))
  cat(sprintf("  %-*s %*s\n", lw, "BIC",                  vw, f2(bic_val)))
  cat(sprintf("  %-*s %*s\n", lw, "LogLik",               vw, f2(loglik_val)))
  cat(sprintf("  %-*s %*s\n", lw, "Parâmetro Espacial",   vw, f4(p_esp)))
  cat(sprintf("  %-*s %*s\n", lw, "I de Moran (Resíduos)",vw, f4(moran_res)))
  cat(sep_m, "\n", sep = "")
  .print_rodape()

  invisible(data.frame(
    Metrica = c("Pseudo_R2", "RMSE", "MAE", "AIC", "BIC", "LogLik", "Param_Espacial", "Moran_Residuos"),
    Valor   = c(r2, rmse, mae, aic_val, bic_val, loglik_val, p_esp, moran_res),
    stringsAsFactors = FALSE
  ))
}





# ── 7. SUAVIZADOR BAYESIANO EMPÍRICO ESPACIAL (DADOS DE ÁREA) ──────────────────

#' @title Suavizador Bayesiano Empírico Espacial (Taxas em Pequenas Áreas)
#' @description Aplica os estimadores Bayesianos Empíricos Global (EBG) e Local (EBL)
#'   de Marshall (1991) para suavizar a variabilidade espúria de taxas de eventos raros
#'   ou mortalidade em áreas geográficas com populações pequenas e desiguais (shrinkage).
#'
#' @param dados Data frame ou objeto \code{sf} contendo as colunas de eventos e população.
#' @param eventos Nome da coluna com a contagem de casos/eventos (string ou símbolo).
#' @param populacao Nome da coluna com a população exposta ao risco (string ou símbolo).
#' @param metodo Método de suavização: \code{"auto"} (seleciona local se coordenadas/pesos existirem,
#'   senão global), \code{"local"} (encolhimento em direção à vizinhança geográfica) ou
#'   \code{"global"} (encolhimento em direção à média global).
#' @param matriz_pesos Matriz quadrada de vizinhança espacial \eqn{W} opcional.
#' @param k_vizinhos Número de vizinhos mais próximos se for necessário construir a matriz de pesos (padrão 5).
#' @param coordenadas Vetor com os nomes das colunas de coordenadas para matriz local se \code{matriz_pesos} for \code{NULL}.
#' @param base Fator multiplicativo da taxa (padrão \code{100000} para taxas por 100 mil habitantes, ou \code{1000}, \code{1}).
#' @param grafico Lógico. Se TRUE (padrão), plota gráfico de dispersão com efeito funil (shrinkage) e densidade comparativa.
#'
#' @return Objeto invisível da classe \code{estatR_suavizador_bayesiano}.
#' @export
suavizador_bayesiano_espacial <- function(dados,
                                          eventos = NULL,
                                          populacao = NULL,
                                          metodo = c("auto", "global", "local"),
                                          matriz_pesos = NULL,
                                          k_vizinhos = 5,
                                          coordenadas = NULL,
                                          base = 100000,
                                          grafico = TRUE) {

  metodo <- match.arg(metodo)
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data.frame.", call. = FALSE)
  df <- as.data.frame(dados)

  # Identificação das colunas de eventos e população
  nome_ev <- if (!is.null(eventos)) as.character(substitute(eventos))[1] else NULL
  if (is.null(nome_ev) || !nome_ev %in% names(df)) {
    if (is.character(eventos) && eventos %in% names(df)) {
      nome_ev <- eventos
    } else {
      # Procura colunas como 'casos', 'obitos', 'eventos', 'mortes'
      cands <- grep("^(casos|obitos|eventos|mortes|y|count|n_casos)$", names(df), ignore.case = TRUE, value = TRUE)
      if (length(cands) > 0) nome_ev <- cands[1] else stop("Coluna de 'eventos' não identificada.", call. = FALSE)
    }
  }

  nome_pop <- if (!is.null(populacao)) as.character(substitute(populacao))[1] else NULL
  if (is.null(nome_pop) || !nome_pop %in% names(df)) {
    if (is.character(populacao) && populacao %in% names(df)) {
      nome_pop <- populacao
    } else {
      # Procura colunas como 'pop', 'populacao', 'habitantes', 'n'
      cands <- grep("^(pop|populacao|habitantes|pop_total|n_pop)$", names(df), ignore.case = TRUE, value = TRUE)
      if (length(cands) > 0) nome_pop <- cands[1] else stop("Coluna de 'populacao' não identificada.", call. = FALSE)
    }
  }

  y <- as.numeric(df[[nome_ev]])
  n_pop <- as.numeric(df[[nome_pop]])

  # Validações
  val <- stats::complete.cases(y, n_pop) & n_pop > 0 & y >= 0
  if (sum(val) < 4) stop("São necessárias pelo menos 4 observações válidas com população positiva.", call. = FALSE)
  df_val <- df[val, , drop = FALSE]
  y <- y[val]
  n_pop <- n_pop[val]
  N <- length(y)

  # Taxa bruta observada
  r_bruta <- (y / n_pop)

  # Estimador Global de Marshall (1991)
  # mu = média global ponderada pela população
  mu_global <- sum(y) / sum(n_pop)
  s2_global <- sum(n_pop * (r_bruta - mu_global)^2) / sum(n_pop)
  n_barra   <- sum(n_pop) / N
  sigma2_global <- max(1e-10, s2_global - (mu_global / n_barra))

  # Pesos bayesianos globais
  w_global <- sigma2_global / (sigma2_global + (mu_global / n_pop))
  taxa_ebg <- w_global * r_bruta + (1 - w_global) * mu_global

  # Estimador Local (se solicitado ou auto com vizinhança)
  tem_vizinhança <- !is.null(matriz_pesos) || !is.null(coordenadas) || any(grepl("^(x|y|lon|lat)$", names(df_val), ignore.case = TRUE))
  metodo_efetivo <- metodo
  if (metodo == "auto") {
    metodo_efetivo <- if (tem_vizinhança) "local" else "global"
  }

  taxa_ebl <- NULL
  w_local  <- NULL
  if (metodo_efetivo == "local") {
    W <- NULL
    if (!is.null(matriz_pesos)) {
      W <- as.matrix(matriz_pesos)[val, val, drop = FALSE]
    } else {
      # Cria matriz knn a partir de coordenadas
      coords_df <- NULL
      if (!is.null(coordenadas)) {
        coords_df <- df_val[, coordenadas, drop = FALSE]
      } else {
        cand_c <- grep("^(x|y|lon|lat|longitude|latitude)$", names(df_val), ignore.case = TRUE, value = TRUE)
        if (length(cand_c) >= 2) coords_df <- df_val[, cand_c[1:2], drop = FALSE]
      }
      if (!is.null(coords_df) && ncol(coords_df) == 2) {
        d_mat <- as.matrix(stats::dist(as.matrix(coords_df)))
        diag(d_mat) <- Inf
        k_ef <- min(as.integer(k_vizinhos), N - 1)
        W <- matrix(0, nrow = N, ncol = N)
        for (i in seq_len(N)) {
          vzs <- order(d_mat[i, ])[seq_len(k_ef)]
          W[i, vzs] <- 1
        }
      }
    }

    if (!is.null(W)) {
      taxa_ebl <- numeric(N)
      w_local  <- numeric(N)
      for (i in seq_len(N)) {
        vzs <- which(W[i, ] > 0)
        vzs_com_i <- unique(c(i, vzs))
        y_loc <- y[vzs_com_i]
        n_loc <- n_pop[vzs_com_i]
        mu_loc <- sum(y_loc) / sum(n_loc)
        s2_loc <- sum(n_loc * (r_bruta[vzs_com_i] - mu_loc)^2) / sum(n_loc)
        n_bar_loc <- mean(n_loc)
        sigma2_loc <- max(1e-10, s2_loc - (mu_loc / n_bar_loc))
        w_i <- sigma2_loc / (sigma2_loc + (mu_loc / n_pop[i]))
        w_i <- max(0, min(1, w_i))
        w_local[i]  <- w_i
        taxa_ebl[i] <- w_i * r_bruta[i] + (1 - w_i) * mu_loc
      }
    } else {
      metodo_efetivo <- "global"
    }
  }

  taxa_final <- if (metodo_efetivo == "local") taxa_ebl else taxa_ebg
  w_final    <- if (metodo_efetivo == "local") w_local else w_global

  # Multiplicação pela base
  taxa_bruta_b <- r_bruta * base
  taxa_suav_b  <- taxa_final * base

  # Variâncias para diagnóstico
  var_bruta <- stats::var(taxa_bruta_b)
  var_suav  <- stats::var(taxa_suav_b)
  red_var   <- if (var_bruta > 0) ((var_bruta - var_suav) / var_bruta) * 100 else 0

  # Estatísticas descritivas das taxas
  q_b <- stats::quantile(taxa_bruta_b, probs = c(0, 0.25, 0.5, 0.75, 1))
  q_s <- stats::quantile(taxa_suav_b, probs = c(0, 0.25, 0.5, 0.75, 1))

  # ── Console Output ──────────────────────────────────────────────────────────
  nome_met_pt <- if (metodo_efetivo == "local") "Bayesiano Empírico Local (EBL)" else "Bayesiano Empírico Global (EBG)"
  .print_titulo("SUAVIZADOR BAYESIANO ESPACIAL (TAXAS EM PEQUENAS ÁREAS)")
  cat(sprintf("  Variável de Eventos : %s (Total: %s casos)\n", nome_ev, .esp_fmt_num(sum(y), 0)))
  cat(sprintf("  População sob Risco : %s (Total: %s habitantes)\n", nome_pop, .esp_fmt_num(sum(n_pop), 0)))
  cat(sprintf("  Unidades Espaciais  : %d áreas geográficas\n", N))
  cat(sprintf("  Fator de Base       : por %s habitantes\n", .esp_fmt_num(base, 0)))
  cat(sprintf("  Método de Suavização: %s\n\n", nome_met_pt))

  .print_topico("COMPARAÇÃO DISTRIBUTIVA: TAXA BRUTA VS SUAVIZADA")
  w_col <- c(est = 18, bruta = 18, suav = 20)
  cat(.esp_sep(sum(w_col)), "\n")
  cat("  ",
      .esp_pad("Estatística", w_col["est"]),
      .esp_pad("Taxa Bruta", w_col["bruta"], "right"),
      .esp_pad("Taxa Suavizada", w_col["suav"], "right"),
      "\n", sep = "")
  cat(.esp_sep(sum(w_col)), "\n")

  imprime_linha <- function(rotulo, vb, vs) {
    cat("  ",
        .esp_pad(rotulo, w_col["est"]),
        .esp_pad(.esp_fmt_num(vb, 2), w_col["bruta"], "right"),
        .esp_pad(.esp_fmt_num(vs, 2), w_col["suav"], "right"),
        "\n", sep = "")
  }

  imprime_linha("Mínimo", q_b[1], q_s[1])
  imprime_linha("1º Quartil", q_b[2], q_s[2])
  imprime_linha("Mediana", q_b[3], q_s[3])
  imprime_linha("Média Ponderada", mean(taxa_bruta_b), mean(taxa_suav_b))
  imprime_linha("3º Quartil", q_b[4], q_s[4])
  imprime_linha("Máximo", q_b[5], q_s[5])
  imprime_linha("Desvio Padrão", stats::sd(taxa_bruta_b), stats::sd(taxa_suav_b))
  cat(.esp_sep(sum(w_col)), "\n\n")

  .print_topico("DIAGNÓSTICO DE ENCOLHIMENTO (SHRINKAGE)")
  cat(sprintf("  Variância da Taxa Bruta     : %s\n", .esp_fmt_num(var_bruta, 2)))
  cat(sprintf("  Variância da Taxa Suavizada : %s\n", .esp_fmt_num(var_suav, 2)))
  cat(sprintf("  Redução de Ruído Amostral   : %s%%\n", .esp_fmt_num(red_var, 1)))
  cat(sprintf("  Peso Médio Local (w̄)        : %s  (0 = 100%% encolhimento | 1 = 100%% taxa local)\n\n",
              .esp_fmt_num(mean(w_final), 3)))

  if (red_var >= 10) {
    cat("  [✓] A suavização bayesiana reduziu com sucesso a instabilidade espúria causada\n")
    cat("      por pequenas populações, preservando a variabilidade genuína das áreas maiores.\n\n")
  } else {
    cat("  [✓] As taxas brutas já apresentavam alta estabilidade (populações suficientemente grandes).\n\n")
  }
  .print_rodape()

  # ── Gráficos em ggplot2 com meu_tema(grade = "dupla") ───────────────────────
  p_combined <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot <- data.frame(
      Populacao      = n_pop,
      Taxa_Bruta     = taxa_bruta_b,
      Taxa_Suavizada = taxa_suav_b,
      Peso           = w_final
    )

    # Painel 1: Efeito Funil (População vs Taxa)
    p1 <- ggplot2::ggplot(df_plot) +
      ggplot2::geom_point(ggplot2::aes(x = Populacao, y = Taxa_Bruta),
                          color = "gray60", size = 2.2, alpha = 0.6) +
      ggplot2::geom_point(ggplot2::aes(x = Populacao, y = Taxa_Suavizada),
                          color = "#1D3557", size = 2.4, alpha = 0.85) +
      ggplot2::scale_x_log10(labels = function(x) formatC(x, big.mark = ".", decimal.mark = ",")) +
      ggplot2::labs(
        title    = "Efeito Funil e Encolhimento Bayesiano (Shrinkage)",
        subtitle = "Pontos Cinzas: Taxa Bruta | Pontos Azuis: Taxa Suavizada",
        x        = "População sob Risco (Escala Logarítmica)",
        y        = sprintf("Taxa por %s hab.", .esp_fmt_num(base, 0))
      ) +
      meu_tema(grade = "dupla")

    # Painel 2: Densidades Comparativas
    df_dens <- data.frame(
      Taxa = c(taxa_bruta_b, taxa_suav_b),
      Tipo = factor(rep(c("Taxa Bruta", "Taxa Suavizada"), each = N),
                    levels = c("Taxa Bruta", "Taxa Suavizada"))
    )

    p2 <- ggplot2::ggplot(df_dens, ggplot2::aes(x = Taxa, color = Tipo, fill = Tipo)) +
      ggplot2::geom_density(alpha = 0.25, linewidth = 1.0) +
      ggplot2::scale_color_manual(values = c("Taxa Bruta" = "#E63946", "Taxa Suavizada" = "#1D3557")) +
      ggplot2::scale_fill_manual(values = c("Taxa Bruta" = "#E63946", "Taxa Suavizada" = "#457B9D")) +
      ggplot2::labs(
        title    = "Distribuição de Densidade das Taxas",
        subtitle = sprintf("Redução de variância de %s%% na cauda de outliers", .esp_fmt_num(red_var, 1)),
        x        = sprintf("Taxa por %s hab.", .esp_fmt_num(base, 0)),
        y        = "Densidade",
        color    = "", fill = ""
      ) +
      meu_tema(grade = "dupla") +
      ggplot2::theme(legend.position = "top")

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      p_combined <- list(painel_funil = p1, painel_densidade = p2)
    } else {
      print(p1)
      p_combined <- p1
    }
  }

  res <- list(
    dados_suavizados = data.frame(
      eventos        = y,
      populacao      = n_pop,
      taxa_bruta     = taxa_bruta_b,
      taxa_suavizada = taxa_suav_b,
      peso_bayesiano = w_final
    ),
    metodo           = metodo_efetivo,
    base             = base,
    reducao_var_pct  = red_var,
    var_bruta        = var_bruta,
    var_suavizada    = var_suav,
    n                = N,
    grafico          = p_combined
  )
  class(res) <- c("estatR_suavizador_bayesiano", "list")
  invisible(res)
}

#' @export
taxa_bayesiana_espacial <- suavizador_bayesiano_espacial




# ==============================================================================
# PARTE II: PROCESSOS PONTUAIS ESPACIAIS (POINT PATTERN ANALYSIS)
# ==============================================================================

# ── HELPERS INTERNOS PARA PROCESSOS PONTUAIS ──────────────────────────────────

#' @noRd
.esp_extrair_coords_pontos <- function(dados, coordenadas = NULL) {
  if (inherits(dados, "sf")) {
    if (requireNamespace("sf", quietly = TRUE)) {
      coords <- sf::st_coordinates(dados)[, 1:2, drop = FALSE]
      colnames(coords) <- c("x", "y")
      return(as.data.frame(coords))
    }
  }
  if (!is.null(coordenadas)) {
    if (is.character(coordenadas) && length(coordenadas) == 2) {
      if (!all(coordenadas %in% colnames(dados))) {
        stop(sprintf("Colunas de coordenadas '%s' e '%s' não encontradas nos dados.",
                     coordenadas[1], coordenadas[2]), call. = FALSE)
      }
      res <- as.data.frame(dados[, coordenadas])
      colnames(res) <- c("x", "y")
      return(res)
    } else if (is.numeric(coordenadas) && length(coordenadas) == 2) {
      res <- as.data.frame(dados[, coordenadas])
      colnames(res) <- c("x", "y")
      return(res)
    }
  }
  if ((is.data.frame(dados) || is.matrix(dados)) && ncol(dados) >= 2) {
    cnames <- tolower(colnames(dados))
    idx_x <- which(cnames %in% c("x", "lon", "longitude", "coord_x", "x_coord", "long"))[1]
    idx_y <- which(cnames %in% c("y", "lat", "latitude", "coord_y", "y_coord"))[1]
    if (!is.na(idx_x) && !is.na(idx_y)) {
      res <- as.data.frame(dados[, c(idx_x, idx_y)])
      colnames(res) <- c("x", "y")
      return(res)
    }
    res <- as.data.frame(dados[, 1:2])
    colnames(res) <- c("x", "y")
    return(res)
  }
  stop("Não foi possível identificar as coordenadas dos pontos. Informe 'coordenadas = c(\"x\", \"y\")'.", call. = FALSE)
}

#' @noRd
.esp_definir_janela <- function(coords, janela = "auto") {
  if (is.character(janela) && janela == "auto") {
    rx <- range(coords$x, na.rm = TRUE)
    ry <- range(coords$y, na.rm = TRUE)
    dx <- rx[2] - rx[1]
    dy <- ry[2] - ry[1]
    pad_x <- if (dx > 0) dx * 0.02 else 0.1
    pad_y <- if (dy > 0) dy * 0.02 else 0.1
    return(c(xmin = rx[1] - pad_x, xmax = rx[2] + pad_x,
             ymin = ry[1] - pad_y, ymax = ry[2] + pad_y))
  }
  if (is.numeric(janela) && length(janela) == 4) {
    names(janela) <- c("xmin", "xmax", "ymin", "ymax")
    return(janela)
  }
  stop("O argumento 'janela' deve ser 'auto' ou um vetor numérico c(xmin, xmax, ymin, ymax).", call. = FALSE)
}


# ── 8. TESTE DE QUADRATS (CSR) ────────────────────────────────────────────────

#' @title Teste de Quadrats para Homogeneidade Espacial (CSR)
#' @description Avalia a hipótese de Aleatoriedade Espacial Completa (CSR - Complete Spatial Randomness)
#'   dividindo o domínio geográfico em uma grade de quadrats e comparando as frequências observadas
#'   com uma distribuição teórica homogênea de Poisson via teste Qui-Quadrado (χ²).
#' @param dados Data frame, matriz ou objeto \code{sf} contendo as coordenadas dos pontos.
#' @param coordenadas Vetor de caracteres ou índices com o nome das colunas de coordenadas (ex: \code{c("x", "y")}).
#' @param nx Número de divisões no eixo X ou \code{"auto"} para determinação automática da grade ótima.
#' @param ny Número de divisões no eixo Y ou \code{"auto"} para determinação automática da grade ótima.
#' @param janela Vetor com os limites do retângulo envolvente \code{c(xmin, xmax, ymin, ymax)} ou \code{"auto"}.
#' @param grafico Lógico. Se TRUE, plota o gráfico em ggplot2 com a distribuição de frequências e \code{meu_tema(grade = "dupla")}.
#' @return Objeto invisível da classe \code{estatR_teste_quadrats}.
#' @export
teste_quadrats <- function(dados,
                           coordenadas = NULL,
                           nx = "auto",
                           ny = "auto",
                           janela = "auto",
                           grafico = TRUE) {

  coords <- .esp_extrair_coords_pontos(dados, coordenadas)
  N <- nrow(coords)
  if (N < 4) stop("São necessárias pelo menos 4 observações para o teste de quadrats.", call. = FALSE)

  j <- .esp_definir_janela(coords, janela)
  Lx <- j["xmax"] - j["xmin"]
  Ly <- j["ymax"] - j["ymin"]

  nx_auto <- identical(nx, "auto")
  ny_auto <- identical(ny, "auto")
  if (nx_auto || ny_auto) {
    M_alvo <- max(4, min(64, round(N / 6)))
    razao <- Lx / Ly
    if (razao <= 0) razao <- 1
    if (nx_auto) {
      nx_val <- max(2, min(12, round(sqrt(M_alvo * razao))))
    } else {
      nx_val <- max(2, as.integer(nx))
    }
    if (ny_auto) {
      ny_val <- max(2, min(12, round(M_alvo / nx_val)))
    } else {
      ny_val <- max(2, as.integer(ny))
    }
  } else {
    nx_val <- max(2, as.integer(nx))
    ny_val <- max(2, as.integer(ny))
  }

  m <- nx_val * ny_val
  cortes_x <- seq(j["xmin"], j["xmax"], length.out = nx_val + 1)
  cortes_y <- seq(j["ymin"], j["ymax"], length.out = ny_val + 1)

  fx <- cut(coords$x, breaks = cortes_x, include.lowest = TRUE)
  fy <- cut(coords$y, breaks = cortes_y, include.lowest = TRUE)
  tab_cont <- as.data.frame(table(fx, fy))
  contagens <- tab_cont$Freq

  E_esperado <- N / m
  chi2 <- sum((contagens - E_esperado)^2 / E_esperado)
  gl <- m - 1
  p_valor <- stats::pchisq(chi2, df = gl, lower.tail = FALSE)
  vmr <- stats::var(contagens) / mean(contagens)

  q_vals <- stats::quantile(contagens, probs = c(0, 0.25, 0.5, 0.75, 1))
  dp_vals <- stats::sd(contagens)

  .print_titulo("TESTE DE HOMOGENEIDADE ESPACIAL: QUADRATS (CSR)")
  cat(sprintf("  Nº Observações    : %d pontos\n", N))
  cat(sprintf("  Janela de Estudo  : [%s, %s] × [%s, %s]\n",
              .esp_fmt_num(j["xmin"], 2), .esp_fmt_num(j["xmax"], 2),
              .esp_fmt_num(j["ymin"], 2), .esp_fmt_num(j["ymax"], 2)))
  modo_grade <- if (nx_auto && ny_auto) " (selecionada automaticamente)" else ""
  cat(sprintf("  Grade de Quadrats : %d × %d (%d células)%s\n", nx_val, ny_val, m, modo_grade))
  cat(sprintf("  Média Esperada (E): %s pontos por quadrat\n\n", .esp_fmt_num(E_esperado, 2)))

  .print_topico("DISTRIBUIÇÃO DAS CONTAGENS POR QUADRAT")
  w_q <- c(min = 10, q1 = 12, med = 12, medm = 12, q3 = 12, max = 10, dp = 14)
  w_tot <- sum(w_q)
  cat(.esp_sep(w_tot), "\n")
  cat("  ",
      .esp_pad("Mínimo", w_q["min"], "center"),
      .esp_pad("1º Quartil", w_q["q1"], "center"),
      .esp_pad("Mediana", w_q["med"], "center"),
      .esp_pad("Média", w_q["medm"], "center"),
      .esp_pad("3º Quartil", w_q["q3"], "center"),
      .esp_pad("Máximo", w_q["max"], "center"),
      .esp_pad("Desvio Padrão", w_q["dp"], "center"),
      "\n", sep = "")
  cat(.esp_sep(w_tot), "\n")
  cat("  ",
      .esp_pad(q_vals[1], w_q["min"], "center"),
      .esp_pad(.esp_fmt_num(q_vals[2], 1), w_q["q1"], "center"),
      .esp_pad(.esp_fmt_num(q_vals[3], 1), w_q["med"], "center"),
      .esp_pad(.esp_fmt_num(mean(contagens), 2), w_q["medm"], "center"),
      .esp_pad(.esp_fmt_num(q_vals[4], 1), w_q["q3"], "center"),
      .esp_pad(q_vals[5], w_q["max"], "center"),
      .esp_pad(.esp_fmt_num(dp_vals, 2), w_q["dp"], "center"),
      "\n", sep = "")
  cat(.esp_sep(w_tot), "\n\n")

  .print_topico("RESULTADOS DO TESTE QUI-QUADRADO")
  w_r <- 65
  cat(sprintf("  Estatística Qui-Quadrado (χ²) : %s\n", .esp_fmt_num(chi2, 3)))
  cat(sprintf("  Graus de Liberdade (gl)       : %d\n", gl))
  desc_vmr <- if (vmr > 1.2) " (superdispersão / agrupamento)" else if (vmr < 0.8) " (subdispersão / regularidade)" else " (compatível com Poisson)"
  cat(sprintf("  Razão Variância/Média (VMR)   : %s%s\n", .esp_fmt_num(vmr, 3), desc_vmr))
  p_str <- paste0(.esp_fmt_p(p_valor), " ", .esp_asterisco(p_valor))
  cat(sprintf("  p-valor                       : %s\n", p_str))
  cat(.esp_sep(w_r), "\n")
  .esp_nota_sig()

  .print_topico("INTERPRETAÇÃO")
  if (p_valor < 0.05) {
    if (vmr > 1) {
      cat(sprintf("  [✓] A hipótese de Aleatoriedade Espacial Completa (CSR) foi rejeitada (p = %s).\n", .esp_fmt_p(p_valor)))
      cat(sprintf("  A razão variância/média elevada (VMR = %s) indica um padrão espacial Agrupado (clusters),\n  com concentrações densas de eventos em determinados quadrantes e outros rarefeitos.\n\n",
                  .esp_fmt_num(vmr, 3)))
    } else {
      cat(sprintf("  [✓] A hipótese de Aleatoriedade Espacial Completa (CSR) foi rejeitada (p = %s).\n", .esp_fmt_p(p_valor)))
      cat(sprintf("  A razão variância/média reduzida (VMR = %s) indica um padrão espacial Regular / Uniforme,\n  evidenciando dispersão ou inibição mútua entre os eventos.\n\n",
                  .esp_fmt_num(vmr, 3)))
    }
  } else {
    cat(sprintf("  [!] Não há evidência estatística para rejeitar a Aleatoriedade Espacial Completa (CSR) (p = %s).\n", .esp_fmt_p(p_valor)))
    cat("  O padrão espacial é compatível com um Processo Pontual de Poisson Homogêneo constante.\n\n")
  }
  .print_rodape()

  p_plot <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot <- data.frame(Contagem = contagens)
    k_range <- 0:max(contagens, stats::qpois(0.999, lambda = E_esperado))
    df_pois <- data.frame(
      k = k_range,
      freq_esperada = stats::dpois(k_range, lambda = E_esperado) * m
    )

    p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Contagem)) +
      ggplot2::geom_histogram(binwidth = 1, fill = "#457B9D", color = "white", alpha = 0.85,
                              center = 0, ggplot2::aes(y = ggplot2::after_stat(count))) +
      ggplot2::geom_line(data = df_pois, ggplot2::aes(x = k, y = freq_esperada),
                         color = "#E63946", linewidth = 1.1, linetype = "dashed") +
      ggplot2::geom_point(data = df_pois, ggplot2::aes(x = k, y = freq_esperada),
                          color = "#E63946", size = 2) +
      ggplot2::labs(
        title    = "Teste de Quadrats: Frequência Observada vs Teórica de Poisson",
        subtitle = sprintf("Grade %d×%d | χ² = %s (gl = %d, p = %s) | Linha vermelha: Poisson E(k) = %s",
                           nx_val, ny_val, .esp_fmt_num(chi2, 2), gl, .esp_fmt_p(p_valor), .esp_fmt_num(E_esperado, 1)),
        x        = "Número de Pontos por Quadrat",
        y        = "Frequência de Quadrats"
      ) +
      meu_tema(grade = "dupla")

    suppressMessages(print(p_plot))
  }

  res <- list(
    estatistica   = chi2,
    gl            = gl,
    p_valor       = p_valor,
    vmr           = vmr,
    esperado      = E_esperado,
    contagens     = contagens,
    nx            = nx_val,
    ny            = ny_val,
    janela        = j,
    grafico       = p_plot
  )
  class(res) <- c("estatR_teste_quadrats", "list")
  invisible(res)
}


# ── 9. TESTE DE CLARK-EVANS (CSR) ─────────────────────────────────────────────

#' @title Teste de Agregação de Clark-Evans (CSR)
#' @description Avalia a dependência espacial de processos pontuais comparando a distância média
#'   observada entre cada evento e seu vizinho mais próximo com a distância teórica esperada sob
#'   Aleatoriedade Espacial Completa (CSR - Processo de Poisson Homogêneo).
#' @param dados Data frame, matriz ou objeto \code{sf} com coordenadas dos pontos.
#' @param coordenadas Vetor de caracteres ou índices com as colunas de coordenadas (ex: \code{c("x", "y")}).
#' @param janela Vetor \code{c(xmin, xmax, ymin, ymax)} ou \code{"auto"}.
#' @param correcao Método de correção de borda: \code{"auto"} (Donnelly para janela retangular),
#'   \code{"donnelly"} ou \code{"nenhuma"}.
#' @param grafico Lógico. Se TRUE, gera gráfico em ggplot2 com \code{meu_tema(grade = "dupla")}.
#' @return Objeto invisível da classe \code{estatR_clark_evans}.
#' @export
teste_clark_evans <- function(dados,
                              coordenadas = NULL,
                              janela = "auto",
                              correcao = c("auto", "donnelly", "nenhuma"),
                              grafico = TRUE) {

  correcao <- match.arg(correcao)
  coords <- .esp_extrair_coords_pontos(dados, coordenadas)
  N <- nrow(coords)
  if (N < 4) stop("São necessárias pelo menos 4 observações para o teste de Clark-Evans.", call. = FALSE)

  j <- .esp_definir_janela(coords, janela)
  Lx <- j["xmax"] - j["xmin"]
  Ly <- j["ymax"] - j["ymin"]
  area <- Lx * Ly
  perimetro <- 2 * (Lx + Ly)
  lambda <- N / area

  coords_mat <- as.matrix(coords)
  d_mat <- as.matrix(stats::dist(coords_mat))
  diag(d_mat) <- Inf
  d_min <- apply(d_mat, 1, min)
  r_A <- mean(d_min)

  r_E <- 1 / (2 * sqrt(lambda))

  if (correcao %in% c("auto", "donnelly")) {
    r_teorico <- r_E + (0.0514 + 0.041 / sqrt(N)) * (perimetro / N)
    var_teorica <- 0.0703 * area / (N^2) + 0.037 * perimetro * sqrt(area) / (N^2.5)
    se <- sqrt(var_teorica)
    nome_corr <- "Correção de Borda de Donnelly (1978)"
  } else {
    r_teorico <- r_E
    se <- 0.26136 / sqrt(N * lambda)
    nome_corr <- "Sem correção de borda"
  }

  R_ratio <- r_A / r_teorico
  z_score <- (r_A - r_teorico) / se
  p_valor <- 2 * stats::pnorm(-abs(z_score))

  .print_titulo("TESTE DE AGREGAÇÃO DE CLARK-EVANS (CSR)")
  cat(sprintf("  Nº Observações      : %d pontos\n", N))
  cat(sprintf("  Janela de Estudo    : Área = %s | Perímetro = %s\n",
              .esp_fmt_num(area, 2), .esp_fmt_num(perimetro, 2)))
  cat(sprintf("  Densidade Estimada  : %s pontos/unidade²\n", .esp_fmt_num(lambda, 4)))
  cat(sprintf("  Método de Correção  : %s\n\n", nome_corr))

  .print_topico("RESULTADOS DO TESTE")
  w_r <- 65
  cat(sprintf("  Distância Média Observada (r̄A)   : %s\n", .esp_fmt_num(r_A, 4)))
  cat(sprintf("  Distância Média Esperada (r̄E)    : %s\n", .esp_fmt_num(r_teorico, 4)))
  desc_r <- if (R_ratio < 0.9) " (agregação / cluster)" else if (R_ratio > 1.1) " (regularidade / dispersão)" else " (aleatório)"
  cat(sprintf("  Índice de Clark-Evans (R)        : %s%s\n", .esp_fmt_num(R_ratio, 4), desc_r))
  cat(sprintf("  Estatística Z (Score)            : %s\n", .esp_fmt_num(z_score, 3)))
  p_str <- paste0(.esp_fmt_p(p_valor), " ", .esp_asterisco(p_valor))
  cat(sprintf("  p-valor                          : %s\n", p_str))
  cat(.esp_sep(w_r), "\n")
  .esp_nota_sig()

  .print_topico("INTERPRETAÇÃO")
  if (p_valor < 0.05) {
    if (R_ratio < 1) {
      cat(sprintf("  [✓] O padrão espacial apresenta Agrupamento estatisticamente significativo (R = %s; p = %s).\n",
                  .esp_fmt_num(R_ratio, 4), .esp_fmt_p(p_valor)))
      cat("  A distância média observada entre vizinhos é inferior à esperada sob CSR,\n  indicando atração mútua ou concentração espacial dos eventos.\n\n")
    } else {
      cat(sprintf("  [✓] O padrão espacial apresenta Regularidade / Dispersão significativa (R = %s; p = %s).\n",
                  .esp_fmt_num(R_ratio, 4), .esp_fmt_p(p_valor)))
      cat("  A distância média observada entre vizinhos é superior à esperada sob CSR,\n  evidenciando inibição mútua ou competição espacial.\n\n")
    }
  } else {
    cat(sprintf("  [!] Não há evidência estatística para rejeitar a Aleatoriedade Espacial Completa (R = %s; p = %s).\n",
                .esp_fmt_num(R_ratio, 4), .esp_fmt_p(p_valor)))
    cat("  A distribuição das distâncias entre vizinhos comporta-se de forma espacialmente aleatória.\n\n")
  }
  .print_rodape()

  p_plot <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot <- data.frame(distancia = d_min)
    p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = distancia)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                              bins = 15, fill = "#1D3557", color = "white", alpha = 0.75) +
      ggplot2::geom_density(color = "#1D3557", linewidth = 1) +
      ggplot2::geom_vline(xintercept = r_teorico, color = "#E63946", linetype = "dashed", linewidth = 1) +
      ggplot2::geom_vline(xintercept = r_A, color = "#2A9D8F", linewidth = 1.1) +
      ggplot2::labs(
        title    = "Teste de Clark-Evans: Distâncias ao Vizinho Mais Próximo",
        subtitle = sprintf("Verde: Média Observada (%s) | Vermelho Tracejado: Média Esperada CSR (%s) | R = %s (p = %s)",
                           .esp_fmt_num(r_A, 2), .esp_fmt_num(r_teorico, 2), .esp_fmt_num(R_ratio, 3), .esp_fmt_p(p_valor)),
        x        = "Distância ao Vizinho Mais Próximo",
        y        = "Densidade"
      ) +
      meu_tema(grade = "dupla")

    suppressMessages(print(p_plot))
  }

  res <- list(
    r_observado = r_A,
    r_esperado  = r_teorico,
    R           = R_ratio,
    z_score     = z_score,
    p_valor     = p_valor,
    distancias  = d_min,
    correcao    = nome_corr,
    janela      = j,
    grafico     = p_plot
  )
  class(res) <- c("estatR_clark_evans", "list")
  invisible(res)
}


# ── 10. FUNÇÃO K DE RIPLEY E L DE BESAG ───────────────────────────────────────

#' @title Função K de Ripley e L de Besag com Envelopes de Monte Carlo
#' @description Estima a função espacial de segunda ordem K de Ripley e sua forma linearizada
#'   L de Besag (L(r) = sqrt(K(r)/pi) - r) para avaliar a dependência espacial em múltiplas
#'   escalas de distância (r), com faixas de simulação CSR (Complete Spatial Randomness) via envelopes de Monte Carlo.
#' @param dados Data frame, matriz ou objeto \code{sf} com coordenadas dos pontos.
#' @param coordenadas Vetor de caracteres ou índices com as colunas de coordenadas.
#' @param r_max Distância máxima avaliada ou \code{"auto"} para adotar 1/4 da menor dimensão da janela.
#' @param n_passos Número de passos na malha de distâncias avaliadas ou \code{"auto"} (padrão 50).
#' @param janela Vetor \code{c(xmin, xmax, ymin, ymax)} ou \code{"auto"}.
#' @param envelope Lógico. Se TRUE, constrói envelopes de simulação sob CSR via Monte Carlo.
#' @param n_sim Número de simulações de Monte Carlo (padrão 99).
#' @param nivel Nível de confiança para o envelope de simulação (padrão 0.95).
#' @param grafico Lógico. Se TRUE, plota o gráfico da curva L(r) com \code{meu_tema(grade = "dupla")}.
#' @return Objeto invisível da classe \code{estatR_k_ripley}.
#' @export
funcao_k_ripley <- function(dados,
                            coordenadas = NULL,
                            r_max = "auto",
                            n_passos = "auto",
                            janela = "auto",
                            envelope = TRUE,
                            n_sim = 99,
                            nivel = 0.95,
                            grafico = TRUE) {

  coords <- .esp_extrair_coords_pontos(dados, coordenadas)
  N <- nrow(coords)
  if (N < 5) stop("São necessárias pelo menos 5 observações para a função K de Ripley.", call. = FALSE)

  j <- .esp_definir_janela(coords, janela)
  Lx <- j["xmax"] - j["xmin"]
  Ly <- j["ymax"] - j["ymin"]
  area <- Lx * Ly
  lambda <- N / area

  r_max_val <- if (identical(r_max, "auto")) 0.25 * min(Lx, Ly) else as.numeric(r_max)
  n_passos_val <- if (identical(n_passos, "auto")) 50 else as.integer(n_passos)
  r_seq <- seq(r_max_val / n_passos_val, r_max_val, length.out = n_passos_val)

  coords_mat <- as.matrix(coords)
  d_mat <- as.matrix(stats::dist(coords_mat))
  diag(d_mat) <- Inf

  fator_k <- area / (N^2)
  calc_k_vetor <- function(dmat, rvec) {
    sapply(rvec, function(r) {
      fator_k * sum(dmat <= r)
    })
  }

  K_obs <- calc_k_vetor(d_mat, r_seq)
  L_obs <- sqrt(K_obs / pi) - r_seq

  L_inf <- NULL
  L_sup <- NULL
  if (envelope) {
    alpha <- 1 - nivel
    p_inf <- alpha / 2
    p_sup <- 1 - alpha / 2

    sim_L_mat <- matrix(0, nrow = n_sim, ncol = length(r_seq))
    for (s in seq_len(n_sim)) {
      sim_x <- stats::runif(N, min = j["xmin"], max = j["xmax"])
      sim_y <- stats::runif(N, min = j["ymin"], max = j["ymax"])
      sim_d <- as.matrix(stats::dist(cbind(sim_x, sim_y)))
      diag(sim_d) <- Inf
      sim_K <- calc_k_vetor(sim_d, r_seq)
      sim_L_mat[s, ] <- sqrt(sim_K / pi) - r_seq
    }
    L_inf <- apply(sim_L_mat, 2, stats::quantile, probs = p_inf)
    L_sup <- apply(sim_L_mat, 2, stats::quantile, probs = p_sup)
  }

  faixas_agreg <- c()
  faixas_inib  <- c()
  if (envelope) {
    idx_agreg <- which(L_obs > L_sup)
    idx_inib  <- which(L_obs < L_inf)
    if (length(idx_agreg) >= 2) {
      r_agr_min <- r_seq[min(idx_agreg)]
      r_agr_max <- r_seq[max(idx_agreg)]
      faixas_agreg <- c(r_agr_min, r_agr_max)
    }
    if (length(idx_inib) >= 2) {
      r_in_min <- r_seq[min(idx_inib)]
      r_in_max <- r_seq[max(idx_inib)]
      faixas_inib <- c(r_in_min, r_in_max)
    }
  }

  .print_titulo("FUNÇÃO K DE RIPLEY E L DE BESAG (PROCESSOS PONTUAIS)")
  cat(sprintf("  Nº Observações    : %d pontos\n", N))
  cat(sprintf("  Janela de Estudo  : Área = %s\n", .esp_fmt_num(area, 2)))
  modo_rmax <- if (identical(r_max, "auto")) " (ótimo: 1/4 da menor dimensão)" else ""
  cat(sprintf("  Raio Máximo (rmax): %s%s\n", .esp_fmt_num(r_max_val, 2), modo_rmax))
  cat(sprintf("  Malha de Avaliação: %d passos de distância\n", n_passos_val))
  if (envelope) {
    cat(sprintf("  Envelope CSR      : %d simulações Monte Carlo (IC %.0f%%)\n\n", n_sim, nivel * 100))
  } else {
    cat("\n")
  }

  .print_topico("DIAGNÓSTICO POR ESCALA ESPACIAL")
  if (envelope) {
    if (length(faixas_agreg) > 0) {
      cat(sprintf("  Agrupamento Espacial (Cluster): Detectado nas distâncias r entre %s e %s.\n",
                  .esp_fmt_num(faixas_agreg[1], 2), .esp_fmt_num(faixas_agreg[2], 2)))
      cat("  Nessa escala de raio, a quantidade média de eventos vizinhos supera significativamente\n  o envelope teórico esperado sob aleatoriedade espacial completa.\n\n")
    }
    if (length(faixas_inib) > 0) {
      cat(sprintf("  Dispersão / Inibição Espacial : Detectada nas distâncias r entre %s e %s.\n",
                  .esp_fmt_num(faixas_inib[1], 2), .esp_fmt_num(faixas_inib[2], 2)))
      cat("  Nessa escala, a quantidade média de vizinhos é inferior ao esperado por Poisson,\n  indicando segregação ou competição por espaço.\n\n")
    }
    if (length(faixas_agreg) == 0 && length(faixas_inib) == 0) {
      cat("  Aleatoriedade Espacial Completa (CSR): A curva L(r) permaneceu inteiramente contida\n")
      cat("  dentro dos limites do envelope de simulação em todas as distâncias avaliadas.\n")
      cat("  Não há evidência de desvio em relação ao processo pontual homogêneo de Poisson.\n\n")
    }
  } else {
    cat("  Envelope de simulação não solicitado. Para diagnóstico inferencial de agrupamento,\n  utilize envelope = TRUE.\n\n")
  }
  .print_rodape()

  p_plot <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot <- data.frame(r = r_seq, L_obs = L_obs)
    if (envelope) {
      df_plot$L_inf <- L_inf
      df_plot$L_sup <- L_sup
    }

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = r))
    if (envelope) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = L_inf, ymax = L_sup),
                                   fill = "gray75", alpha = 0.45)
    }
    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#E63946", linewidth = 0.9) +
      ggplot2::geom_line(ggplot2::aes(y = L_obs), color = "#1D3557", linewidth = 1.2) +
      ggplot2::labs(
        title    = "Função L de Besag: L(r) = √(K(r)/π) - r",
        subtitle = "Linha Azul: Observado | Linha Vermelha: CSR Teórico (L = 0) | Faixa Cinza: Envelope Monte Carlo",
        x        = "Distância (r)",
        y        = "L(r)"
      ) +
      meu_tema(grade = "dupla")

    p_plot <- p
    suppressMessages(print(p_plot))
  }

  res <- list(
    r         = r_seq,
    K         = K_obs,
    L         = L_obs,
    L_inf     = L_inf,
    L_sup     = L_sup,
    r_max     = r_max_val,
    janela    = j,
    envelope  = envelope,
    grafico   = p_plot
  )
  class(res) <- c("estatR_k_ripley", "list")
  invisible(res)
}

#' @export
funcao_ripley <- funcao_k_ripley


# ── 11. FUNÇÃO G DO VIZINHO MAIS PRÓXIMO ──────────────────────────────────────

#' @title Função G do Vizinho Mais Próximo com Envelopes de Monte Carlo
#' @description Estima a função de distribuição acumulada G(r) das distâncias de cada evento
#'   ao seu vizinho mais próximo, comparando com a curva teórica de Poisson e envelopes CSR de Monte Carlo.
#' @param dados Data frame, matriz ou objeto \code{sf} com coordenadas dos pontos.
#' @param coordenadas Vetor de caracteres ou índices com as colunas de coordenadas.
#' @param r_max Distância máxima avaliada ou \code{"auto"}.
#' @param n_passos Número de passos na malha de distâncias avaliadas ou \code{"auto"} (padrão 50).
#' @param janela Vetor \code{c(xmin, xmax, ymin, ymax)} ou \code{"auto"}.
#' @param envelope Lógico. Se TRUE, constrói envelopes de simulação sob CSR via Monte Carlo.
#' @param n_sim Número de simulações de Monte Carlo (padrão 99).
#' @param nivel Nível de confiança para o envelope de simulação (padrão 0.95).
#' @param grafico Lógico. Se TRUE, plota o gráfico da curva G(r) com \code{meu_tema(grade = "dupla")}.
#' @return Objeto invisível da classe \code{estatR_g_pontos}.
#' @export
funcao_g_pontos <- function(dados,
                            coordenadas = NULL,
                            r_max = "auto",
                            n_passos = "auto",
                            janela = "auto",
                            envelope = TRUE,
                            n_sim = 99,
                            nivel = 0.95,
                            grafico = TRUE) {

  coords <- .esp_extrair_coords_pontos(dados, coordenadas)
  N <- nrow(coords)
  if (N < 5) stop("São necessárias pelo menos 5 observações para a função G.", call. = FALSE)

  j <- .esp_definir_janela(coords, janela)
  Lx <- j["xmax"] - j["xmin"]
  Ly <- j["ymax"] - j["ymin"]
  area <- Lx * Ly
  lambda <- N / area

  coords_mat <- as.matrix(coords)
  d_mat <- as.matrix(stats::dist(coords_mat))
  diag(d_mat) <- Inf
  d_min <- apply(d_mat, 1, min)

  r_max_val <- if (identical(r_max, "auto")) max(d_min) * 1.05 else as.numeric(r_max)
  n_passos_val <- if (identical(n_passos, "auto")) 50 else as.integer(n_passos)
  r_seq <- seq(0, r_max_val, length.out = n_passos_val)

  G_obs <- sapply(r_seq, function(r) mean(d_min <= r))
  G_csr <- 1 - exp(-lambda * pi * (r_seq^2))

  G_inf <- NULL
  G_sup <- NULL
  if (envelope) {
    alpha <- 1 - nivel
    p_inf <- alpha / 2
    p_sup <- 1 - alpha / 2

    sim_G_mat <- matrix(0, nrow = n_sim, ncol = length(r_seq))
    for (s in seq_len(n_sim)) {
      sim_x <- stats::runif(N, min = j["xmin"], max = j["xmax"])
      sim_y <- stats::runif(N, min = j["ymin"], max = j["ymax"])
      sim_d <- as.matrix(stats::dist(cbind(sim_x, sim_y)))
      diag(sim_d) <- Inf
      sim_d_min <- apply(sim_d, 1, min)
      sim_G_mat[s, ] <- sapply(r_seq, function(r) mean(sim_d_min <= r))
    }
    G_inf <- apply(sim_G_mat, 2, stats::quantile, probs = p_inf)
    G_sup <- apply(sim_G_mat, 2, stats::quantile, probs = p_sup)
  }

  .print_titulo("FUNÇÃO G: DISTRIBUIÇÃO DO VIZINHO MAIS PRÓXIMO (CSR)")
  cat(sprintf("  Nº Observações    : %d pontos\n", N))
  cat(sprintf("  Janela de Estudo  : Área = %s\n", .esp_fmt_num(area, 2)))
  cat(sprintf("  Distância Média   : %s\n", .esp_fmt_num(mean(d_min), 4)))
  cat(sprintf("  Distância Máxima  : %s\n", .esp_fmt_num(max(d_min), 4)))
  if (envelope) {
    cat(sprintf("  Envelope CSR      : %d simulações Monte Carlo (IC %.0f%%)\n\n", n_sim, nivel * 100))
  } else {
    cat("\n")
  }

  .print_topico("DIAGNÓSTICO DA FUNÇÃO G")
  if (envelope) {
    idx_agr <- which(G_obs > G_sup)
    idx_inb <- which(G_obs < G_inf)
    if (length(idx_agr) >= 2) {
      cat("  Curva Observada acima do Envelope: Indica agregação espacial de eventos.\n")
      cat(sprintf("  As distâncias entre vizinhos são significativamente menores que sob Poisson (r entre %s e %s).\n\n",
                  .esp_fmt_num(r_seq[min(idx_agr)], 2), .esp_fmt_num(r_seq[max(idx_agr)], 2)))
    } else if (length(idx_inb) >= 2) {
      cat("  Curva Observada abaixo do Envelope: Indica padrão regular / disperso.\n")
      cat(sprintf("  Os eventos mantêm distâncias mútuas maiores que o esperado sob CSR (r entre %s e %s).\n\n",
                  .esp_fmt_num(r_seq[min(idx_inb)], 2), .esp_fmt_num(r_seq[max(idx_inb)], 2)))
    } else {
      cat("  Curva contida no envelope: O padrão espacial é plenamente consistente com\n  um Processo Pontual de Poisson Homogêneo (CSR).\n\n")
    }
  } else {
    cat("  Envelope de simulação não solicitado.\n\n")
  }
  .print_rodape()

  p_plot <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_plot <- data.frame(r = r_seq, G_obs = G_obs, G_csr = G_csr)
    if (envelope) {
      df_plot$G_inf <- G_inf
      df_plot$G_sup <- G_sup
    }

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = r))
    if (envelope) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = G_inf, ymax = G_sup),
                                   fill = "gray75", alpha = 0.45)
    }
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = G_csr), linetype = "dashed", color = "#E63946", linewidth = 0.9) +
      ggplot2::geom_line(ggplot2::aes(y = G_obs), color = "#1D3557", linewidth = 1.2) +
      ggplot2::labs(
        title    = "Função G: Probabilidade Acumulada do Vizinho Mais Próximo",
        subtitle = "Linha Azul: G(r) Observado | Linha Vermelha: Teórico Poisson | Faixa Cinza: Envelope CSR",
        x        = "Distância (r)",
        y        = "G(r)"
      ) +
      meu_tema(grade = "dupla")

    p_plot <- p
    suppressMessages(print(p_plot))
  }

  res <- list(
    r         = r_seq,
    G_obs     = G_obs,
    G_csr     = G_csr,
    G_inf     = G_inf,
    G_sup     = G_sup,
    distancias = d_min,
    janela    = j,
    envelope  = envelope,
    grafico   = p_plot
  )
  class(res) <- c("estatR_g_pontos", "list")
  invisible(res)
}




# ── 12. ESTIMATIVA DE DENSIDADE DE KERNEL ESPACIAL (KDE 2D) ───────────────────

#' @title Estimativa de Densidade de Kernel Espacial 2D (KDE)
#' @description Gera uma superfície contínua de intensidade espacial de eventos pontuais
#'   (mancha / mapa de calor de primeira ordem) utilizando estimador de Kernel bivariado
#'   com raio ótimo automático (regra de Silverman 2D) ou largura de banda manual.
#'
#' @param dados Data frame, matriz ou objeto \code{sf} contendo coordenadas dos pontos.
#' @param coordenadas Vetor com nomes ou índices das colunas de coordenadas (ex: \code{c("x", "y")}).
#' @param raio Raio de busca / largura de banda (\eqn{\tau}): \code{"auto"} (regra de Silverman 2D)
#'   ou valor numérico positivo em unidades de distância dos dados.
#' @param kernel Função de Kernel utilizada: \code{"quatico"} (biweight esférico, padrão em SIG)
#'   ou \code{"gaussiano"}.
#' @param resolucao Resolução da malha regular: \code{"auto"} (grade 50×50 = 2.500 células) ou inteiro.
#' @param janela Limites espaciais \code{c(xmin, xmax, ymin, ymax)} ou \code{"auto"}.
#' @param grafico Lógico. Se TRUE (padrão), exibe painel duplo com o mapa de calor contínuo
#'   e curvas de nível de iso-densidade via \code{meu_tema(grade = "dupla")}.
#'
#' @return Objeto invisível da classe \code{estatR_kernel_espacial}.
#' @export
densidade_kernel_espacial <- function(dados,
                                      coordenadas = NULL,
                                      raio = "auto",
                                      kernel = c("quatico", "gaussiano"),
                                      resolucao = "auto",
                                      janela = "auto",
                                      grafico = TRUE) {

  kernel <- match.arg(kernel)
  coords <- .esp_extrair_coords_pontos(dados, coordenadas)
  N <- nrow(coords)
  if (N < 4) stop("São necessários pelo menos 4 pontos para a densidade de Kernel.", call. = FALSE)

  j <- .esp_definir_janela(coords, janela)
  Lx <- j["xmax"] - j["xmin"]
  Ly <- j["ymax"] - j["ymin"]

  # 1. Determinação da largura de banda (raio tau)
  # Regra de Silverman bivariada robusta
  raio_auto <- identical(raio, "auto")
  if (raio_auto) {
    sx <- stats::sd(coords$x, na.rm = TRUE)
    sy <- stats::sd(coords$y, na.rm = TRUE)
    iqr_x <- stats::IQR(coords$x, na.rm = TRUE) / 1.34
    iqr_y <- stats::IQR(coords$y, na.rm = TRUE) / 1.34
    
    sigma_x <- if (iqr_x > 0 && iqr_x < sx) iqr_x else sx
    sigma_y <- if (iqr_y > 0 && iqr_y < sy) iqr_y else sy
    
    sigma_ref <- min(sigma_x, sigma_y)
    if (sigma_ref <= 0) sigma_ref <- max(sigma_x, sigma_y, 1)
    
    # Regra de Silverman para 2D: h = 0.9 * A * n^(-1/6)
    tau <- 0.9 * sigma_ref * (N^(-1/6))
    if (tau <= 0) tau <- min(Lx, Ly) * 0.1
  } else {
    tau <- as.numeric(raio)
    if (tau <= 0) stop("O argumento 'raio' deve ser estritamente positivo.", call. = FALSE)
  }

  # 2. Resolução da grade
  res_val <- if (identical(resolucao, "auto")) 50L else as.integer(resolucao)
  if (res_val < 5) res_val <- 5L

  gx <- seq(j["xmin"], j["xmax"], length.out = res_val)
  gy <- seq(j["ymin"], j["ymax"], length.out = res_val)
  grade <- expand.grid(x = gx, y = gy)
  n_grade <- nrow(grade)

  # 3. Estimativa de densidade vetorizada
  coords_mat <- as.matrix(coords[, c("x", "y")])
  grade_mat  <- as.matrix(grade[, c("x", "y")])

  dens_vec <- numeric(n_grade)
  bloco_tam <- 1000L
  n_blocos <- ceiling(n_grade / bloco_tam)

  # Área do pixel para integração
  dx_p <- (j["xmax"] - j["xmin"]) / (res_val - 1)
  dy_p <- (j["ymax"] - j["ymin"]) / (res_val - 1)
  area_pixel <- dx_p * dy_p

  for (b in seq_len(n_blocos)) {
    idx_b <- seq((b - 1) * bloco_tam + 1, min(b * bloco_tam, n_grade))
    sub_g <- grade_mat[idx_b, , drop = FALSE]
    nb <- length(idx_b)

    # dist_ij entre cada ponto de sub_g e coords_mat
    dx_b <- outer(coords_mat[, 1], sub_g[, 1], "-")
    dy_b <- outer(coords_mat[, 2], sub_g[, 2], "-")
    D_b  <- sqrt(dx_b^2 + dy_b^2)
    U_b  <- D_b / tau

    if (kernel == "quatico") {
      # Kernel biweight / quártico: K(u) = (3 / pi) * (1 - u^2)^2 para u <= 1
      mask_u <- (U_b <= 1)
      K_val <- (3 / pi) * ((1 - U_b^2)^2) * mask_u
      dens_vec[idx_b] <- colSums(K_val) / (tau^2)
    } else {
      # Kernel gaussiano 2D: (1 / (2*pi)) * exp(-0.5 * u^2)
      K_val <- (1 / (2 * pi)) * exp(-0.5 * (U_b^2))
      dens_vec[idx_b] <- colSums(K_val) / (tau^2)
    }
  }

  grade$densidade <- dens_vec

  # Estatísticas da superfície
  q_d <- stats::quantile(dens_vec, probs = c(0, 0.25, 0.5, 0.75, 1))
  vol_estimado <- sum(dens_vec) * area_pixel

  # Identificação do Hotspot (ponto de pico máximo)
  idx_max <- which.max(dens_vec)
  hotspot_x <- grade$x[idx_max]
  hotspot_y <- grade$y[idx_max]
  pico_dens <- grade$densidade[idx_max]

  # ── Console Output ──────────────────────────────────────────────────────────
  nome_k_pt <- if (kernel == "quatico") "Quártico / Biweight (Compacto)" else "Gaussiano (Cauda Infinita)"
  .print_titulo("ESTIMATIVA DE DENSIDADE DE KERNEL ESPACIAL (KDE 2D)")
  cat(sprintf("  Nº de Eventos       : %d pontos observados\n", N))
  cat(sprintf("  Janela de Estudo    : [%s, %s] × [%s, %s]\n",
              .esp_fmt_num(j["xmin"], 2), .esp_fmt_num(j["xmax"], 2),
              .esp_fmt_num(j["ymin"], 2), .esp_fmt_num(j["ymax"], 2)))
  modo_tau <- if (raio_auto) " (largura de banda ótima de Silverman 2D)" else ""
  cat(sprintf("  Raio de Busca (\u03c4) : %s unidades%s\n", .esp_fmt_num(tau, 3), modo_tau))
  cat(sprintf("  Kernel Selecionado  : %s\n", nome_k_pt))
  cat(sprintf("  Malha de Resolução  : %d × %d (%d células regulares)\n\n", res_val, res_val, n_grade))

  .print_topico("DISTRIBUIÇÃO DA DENSIDADE ESTIMADA")
  w_kd <- c(est = 18, val = 18)
  cat(.esp_sep(sum(w_kd)), "\n")
  cat("  ", .esp_pad("Estatística", w_kd["est"]), .esp_pad("Densidade", w_kd["val"], "right"), "\n", sep = "")
  cat(.esp_sep(sum(w_kd)), "\n")
  cat("  ", .esp_pad("Mínimo", w_kd["est"]), .esp_pad(.esp_fmt_num(q_d[1], 4), w_kd["val"], "right"), "\n", sep = "")
  cat("  ", .esp_pad("1º Quartil", w_kd["est"]), .esp_pad(.esp_fmt_num(q_d[2], 4), w_kd["val"], "right"), "\n", sep = "")
  cat("  ", .esp_pad("Mediana", w_kd["est"]), .esp_pad(.esp_fmt_num(q_d[3], 4), w_kd["val"], "right"), "\n", sep = "")
  cat("  ", .esp_pad("Média", w_kd["est"]), .esp_pad(.esp_fmt_num(mean(dens_vec), 4), w_kd["val"], "right"), "\n", sep = "")
  cat("  ", .esp_pad("3º Quartil", w_kd["est"]), .esp_pad(.esp_fmt_num(q_d[4], 4), w_kd["val"], "right"), "\n", sep = "")
  cat("  ", .esp_pad("Pico Máximo", w_kd["est"]), .esp_pad(.esp_fmt_num(q_d[5], 4), w_kd["val"], "right"), "\n", sep = "")
  cat(.esp_sep(sum(w_kd)), "\n\n")

  .print_topico("PICO DE INTENSIDADE ESPACIAL (HOTSPOT)")
  cat(sprintf("  Densidade Máxima    : %s eventos/unidade²\n", .esp_fmt_num(pico_dens, 4)))
  cat(sprintf("  Local do Hotspot    : X = %s | Y = %s\n\n",
              .esp_fmt_num(hotspot_x, 2), .esp_fmt_num(hotspot_y, 2)))
  .print_rodape()

  # ── Gráficos em ggplot2 com meu_tema(grade = "dupla") ───────────────────────
  p_combined <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    # Painel 1: Heatmap com pontos discretos sobrepostos
    p1 <- ggplot2::ggplot(grade, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = densidade)) +
      ggplot2::scale_fill_gradientn(
        colours = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594"),
        name = "Densidade"
      ) +
      ggplot2::geom_point(data = coords, ggplot2::aes(x = x, y = y),
                          color = "#D90429", size = 1.6, alpha = 0.65) +
      ggplot2::labs(
        title    = "Mapa de Calor de Intensidade Espacial (KDE 2D)",
        subtitle = sprintf("Kernel %s | \u03c4 = %s | Pontos vermelhos: Eventos", kernel, .esp_fmt_num(tau, 2)),
        x        = "Coordenada X",
        y        = "Coordenada Y"
      ) +
      meu_tema(grade = "dupla")

    # Painel 2: Curvas de nível (Iso-densidade)
    p2 <- ggplot2::ggplot(grade, ggplot2::aes(x = x, y = y, z = densidade)) +
      ggplot2::geom_contour(ggplot2::aes(color = ggplot2::after_stat(level)), linewidth = 0.85) +
      ggplot2::scale_color_gradientn(
        colours = c("#457B9D", "#E63946", "#800026"),
        name = "Iso-Nível"
      ) +
      ggplot2::geom_point(ggplot2::aes(x = hotspot_x, y = hotspot_y),
                          color = "#800026", shape = 4, size = 3.5, stroke = 1.5) +
      ggplot2::annotate("text", x = hotspot_x, y = hotspot_y, label = " Hotspot",
                        color = "#800026", fontface = "bold", hjust = -0.1, vjust = -0.5, size = 3.5) +
      ggplot2::labs(
        title    = "Curvas de Nível de Concentração (Iso-Densidade)",
        subtitle = "Linhas concêntricas delimitam núcleos de alta aglomeração",
        x        = "Coordenada X",
        y        = "Coordenada Y"
      ) +
      meu_tema(grade = "dupla")

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      p_combined <- list(painel_heatmap = p1, painel_contornos = p2)
    } else {
      print(p1)
      p_combined <- p1
    }
  }

  res <- list(
    grade          = grade,
    densidade      = dens_vec,
    raio           = tau,
    kernel         = kernel,
    hotspot        = c(x = hotspot_x, y = hotspot_y, pico = pico_dens),
    resolucao      = c(x = res_val, y = res_val),
    janela         = j,
    pontos         = coords,
    n              = N,
    grafico        = p_combined
  )
  class(res) <- c("estatR_kernel_espacial", "list")
  invisible(res)
}

#' @export
kde_espacial <- densidade_kernel_espacial


# ==============================================================================
# PARTE III: SUPERFÍCIE CONTÍNUA (GEOESTATÍSTICA / KRIGAGEM)
# ==============================================================================

# ── HELPER INTERNO DE COORDENADAS E VALOR ─────────────────────────────────────

#' @noRd
.geo_extrair_coords_valor <- function(dados, variavel = NULL, coordenadas = NULL) {
  # 1. Suporte a objeto sf
  if (inherits(dados, "sf")) {
    if (requireNamespace("sf", quietly = TRUE)) {
      coords_mat <- sf::st_coordinates(dados)[, 1:2, drop = FALSE]
      coords <- as.data.frame(coords_mat)
      colnames(coords) <- c("x", "y")
      dados_df <- as.data.frame(dados)
      # Se variável informada
      var_nome <- if (!is.null(variavel)) as.character(variavel)[1] else NULL
      if (is.null(var_nome)) {
        num_cols <- setdiff(names(dados_df), c("geometry", colnames(coords)))
        num_cols <- num_cols[vapply(dados_df[num_cols], is.numeric, logical(1))]
        if (length(num_cols) == 0) stop("Nenhuma variável numérica encontrada no objeto sf.", call. = FALSE)
        var_nome <- num_cols[1]
      }
      z <- as.numeric(dados_df[[var_nome]])
      return(list(coords = coords, z = z, var_nome = var_nome))
    }
  }

  if (!is.data.frame(dados) && !is.matrix(dados)) {
    stop("O argumento 'dados' deve ser um data frame ou matriz.", call. = FALSE)
  }
  df <- as.data.frame(dados)

  # 2. Identificação das coordenadas
  coords <- NULL
  if (!is.null(coordenadas)) {
    if (is.character(coordenadas) && length(coordenadas) == 2) {
      if (!all(coordenadas %in% names(df))) {
        stop(sprintf("Colunas de coordenadas '%s' e '%s' não encontradas.", coordenadas[1], coordenadas[2]), call. = FALSE)
      }
      coords <- df[, coordenadas, drop = FALSE]
      colnames(coords) <- c("x", "y")
    } else if (is.numeric(coordenadas) && length(coordenadas) == 2) {
      coords <- df[, coordenadas, drop = FALSE]
      colnames(coords) <- c("x", "y")
    }
  }

  if (is.null(coords)) {
    cnames <- tolower(names(df))
    idx_x <- which(cnames %in% c("x", "lon", "longitude", "coord_x", "x_coord", "long", "easting"))[1]
    idx_y <- which(cnames %in% c("y", "lat", "latitude", "coord_y", "y_coord", "northing"))[1]
    if (!is.na(idx_x) && !is.na(idx_y) && idx_x != idx_y) {
      coords <- df[, c(idx_x, idx_y), drop = FALSE]
      colnames(coords) <- c("x", "y")
    } else if (ncol(df) >= 3) {
      # Assume as duas primeiras colunas como coordenadas se numéricas
      if (is.numeric(df[[1]]) && is.numeric(df[[2]])) {
        coords <- df[, 1:2, drop = FALSE]
        colnames(coords) <- c("x", "y")
      }
    }
  }

  if (is.null(coords)) {
    stop("Não foi possível identificar as coordenadas. Informe 'coordenadas = c(\"x\", \"y\")'.", call. = FALSE)
  }

  # 3. Identificação da variável z
  var_nome <- NULL
  if (!is.null(variavel)) {
    if (is.character(variavel)) {
      var_nome <- variavel[1]
    } else {
      var_nome <- deparse(substitute(variavel))
    }
    if (!var_nome %in% names(df)) {
      stop(sprintf("Variável '%s' não encontrada nos dados.", var_nome), call. = FALSE)
    }
  } else {
    # Pega primeira coluna numérica que não seja uma das coordenadas
    nomes_coords <- if (!is.null(coordenadas) && is.character(coordenadas)) coordenadas else colnames(coords)
    cols_restantes <- setdiff(names(df), nomes_coords)
    cols_num <- cols_restantes[vapply(df[cols_restantes], is.numeric, logical(1))]
    if (length(cols_num) == 0) {
      stop("Nenhuma coluna numérica encontrada para a variável de resposta. Especifique 'variavel'.", call. = FALSE)
    }
    var_nome <- cols_num[1]
  }

  z <- as.numeric(df[[var_nome]])
  list(coords = coords, z = z, var_nome = var_nome)
}

# ── 11. VARIOGRAMA EMPÍRICO E AJUSTE TEÓRICO ──────────────────────────────────

#' @title Variograma Empírico e Ajuste de Modelo Teórico
#' @description Calcula a semivariância empírica de uma variável espacial contínua
#'   e ajusta modelos teóricos (Esférico, Exponencial ou Gaussiano) via mínimos quadrados
#'   com seleção automática do modelo ótimo por menor Soma de Quadrados dos Resíduos (SQR).
#'
#' @param dados Data frame, matriz ou objeto \code{sf} contendo coordenadas e valores amostrados.
#' @param variavel Nome da coluna com a variável de resposta (string ou símbolo).
#' @param coordenadas Vetor de 2 caracteres ou índices com os nomes das colunas de coordenadas (ex: \code{c("x", "y")}).
#'   Se \code{NULL}, detecta automaticamente colunas com nomes usuais.
#' @param modelo Modelo teórico a ajustar: \code{"auto"} (padrão, testa esférico, exponencial e gaussiano
#'   e seleciona o de menor SQR), \code{"esferico"}, \code{"exponencial"} ou \code{"gaussiano"}.
#' @param n_classes Número de classes de distância (lags) ou \code{"auto"} para determinação automática.
#' @param dist_max Distância máxima considerada ou \code{"auto"} (metade da extensão máxima amostrada).
#' @param grafico Lógico. Se TRUE (padrão), plota o variograma empírico com a curva ajustada e o histograma de pares.
#'
#' @return Objeto invisível da classe \code{estatR_variograma}.
#' @export
variograma <- function(dados,
                       variavel = NULL,
                       coordenadas = NULL,
                       modelo = c("auto", "esferico", "exponencial", "gaussiano"),
                       n_classes = "auto",
                       dist_max = "auto",
                       grafico = TRUE) {

  modelo <- match.arg(modelo)
  extraido <- .geo_extrair_coords_valor(dados, variavel = variavel, coordenadas = coordenadas)
  coords <- extraido$coords
  z <- extraido$z
  var_nome <- extraido$var_nome

  # Remove NAs
  val <- stats::complete.cases(coords$x, coords$y, z)
  coords <- coords[val, , drop = FALSE]
  z <- z[val]
  n <- length(z)

  if (n < 6) stop("São necessárias pelo menos 6 observações para o cálculo do variograma.", call. = FALSE)

  # 1. Distâncias e semivariâncias empíricas entre todos os pares
  coords_mat <- as.matrix(coords)
  d_dist <- stats::dist(coords_mat)
  d_vec <- as.numeric(d_dist)
  p_pares <- t(utils::combn(n, 2))
  diff_vec <- 0.5 * (z[p_pares[, 1]] - z[p_pares[, 2]])^2

  max_d_total <- max(d_vec, na.rm = TRUE)
  dist_max_val <- if (identical(dist_max, "auto")) max_d_total / 2 else as.numeric(dist_max)
  if (dist_max_val <= 0) dist_max_val <- max_d_total / 2

  n_classes_val <- if (identical(n_classes, "auto")) max(5, min(15, floor(n / 3))) else as.integer(n_classes)
  if (n_classes_val < 3) n_classes_val <- 3

  # Filtrar pares até dist_max
  idx_valid <- which(d_vec > 0 & d_vec <= dist_max_val)
  if (length(idx_valid) < 10) {
    stop("Poucos pares encontrados dentro da distância máxima especificada. Aumente 'dist_max'.", call. = FALSE)
  }
  d_sub <- d_vec[idx_valid]
  diff_sub <- diff_vec[idx_valid]

  # Lags uniformes
  cortes <- seq(0, dist_max_val, length.out = n_classes_val + 1)
  classes <- cut(d_sub, breaks = cortes, include.lowest = TRUE)

  # Agrupamento por classe
  lag_h <- tapply(d_sub, classes, mean, na.rm = TRUE)
  gamma_emp <- tapply(diff_sub, classes, mean, na.rm = TRUE)
  n_pares <- as.numeric(table(classes))

  # Remove classes vazias
  mask_ok <- !is.na(lag_h) & !is.na(gamma_emp) & n_pares > 0
  lag_h <- as.numeric(lag_h[mask_ok])
  gamma_emp <- as.numeric(gamma_emp[mask_ok])
  n_pares <- n_pares[mask_ok]

  # 2. Funções teóricas de semivariância
  # C0: pepita (nugget)
  # C: contribuição (sill = C0 + C)
  # a: alcance (range)
  f_esferico <- function(h, C0, C, a) {
    ifelse(h <= a, C0 + C * (1.5 * (h / a) - 0.5 * (h / a)^3), C0 + C)
  }
  f_exponencial <- function(h, C0, C, a) {
    C0 + C * (1 - exp(-3 * h / a))  # alcance prático a (onde atinge 95% do patamar)
  }
  f_gaussiano <- function(h, C0, C, a) {
    C0 + C * (1 - exp(-3 * (h / a)^2))  # alcance prático a (95% do patamar)
  }

  var_z <- stats::var(z, na.rm = TRUE)
  c0_init <- max(0, min(gamma_emp) * 0.5)
  c_init <- max(1e-4, var_z - c0_init)
  a_init <- dist_max_val * 0.4

  ajustar_mod <- function(tipo) {
    fn_obj <- function(par) {
      c0 <- par[1]
      c  <- par[2]
      a  <- par[3]
      if (tipo == "esferico")   g_pred <- f_esferico(lag_h, c0, c, a)
      if (tipo == "exponencial") g_pred <- f_exponencial(lag_h, c0, c, a)
      if (tipo == "gaussiano")   g_pred <- f_gaussiano(lag_h, c0, c, a)
      sum(n_pares * (gamma_emp - g_pred)^2)
    }

    low_b <- c(0, 1e-6, 1e-6)
    upp_b <- c(var_z * 2, var_z * 4, dist_max_val * 3)
    p0 <- c(c0_init, c_init, a_init)
    p0 <- pmax(low_b + 1e-5, pmin(upp_b - 1e-5, p0))

    opt <- tryCatch({
      stats::optim(p0, fn_obj, method = "L-BFGS-B", lower = low_b, upper = upp_b)
    }, error = function(e) {
      stats::optim(p0, fn_obj, method = "Nelder-Mead")
    })

    c0_est <- max(0, opt$par[1])
    c_est  <- max(1e-6, opt$par[2])
    a_est  <- max(1e-6, opt$par[3])
    sqr    <- opt$value
    list(tipo = tipo, c0 = c0_est, c = c_est, a = a_est, sqr = sqr)
  }

  modelos_testar <- if (modelo == "auto") c("esferico", "exponencial", "gaussiano") else modelo
  ajustes <- lapply(modelos_testar, ajustar_mod)
  sqrs <- vapply(ajustes, function(x) x$sqr, numeric(1))
  melhor_idx <- which.min(sqrs)
  melhor <- ajustes[[melhor_idx]]

  tipo_selecionado <- melhor$tipo
  C0 <- melhor$c0
  C <- melhor$c
  a <- melhor$a
  sill <- C0 + C
  sqr_val <- melhor$sqr

  # Razão de dependência espacial: C0 / (C0 + C)
  razao_pepita <- (C0 / sill) * 100

  # Nomes amigáveis
  nome_mod_pt <- switch(tipo_selecionado,
    "esferico"    = "Esférico",
    "exponencial" = "Exponencial",
    "gaussiano"   = "Gaussiano"
  )

  # ── Saída no console ────────────────────────────────────────────────────────
  .print_titulo("VARIOGRAMA EMPÍRICO E MODELO TEÓRICO")
  cat(sprintf("  Variável            : %s\n", var_nome))
  cat(sprintf("  Nº Observações      : %d pontos\n", n))
  modo_classes <- if (identical(n_classes, "auto")) " (selecionadas automaticamente)" else ""
  cat(sprintf("  Nº Classes (Lag)    : %d%s\n", length(lag_h), modo_classes))
  modo_dist <- if (identical(dist_max, "auto")) " (metade da extensão total)" else ""
  cat(sprintf("  Distância Máxima    : %s%s\n\n", .esp_fmt_num(dist_max_val, 2), modo_dist))

  .print_topico("PARÂMETROS DO MODELO AJUSTADO")
  w_box <- 68
  cat(.esp_sep(w_box), "\n")
  info_auto <- if (modelo == "auto") " (selecionado automaticamente por menor SQR)" else ""
  cat(sprintf("  Modelo Selecionado  : %s%s\n", nome_mod_pt, info_auto))
  cat(sprintf("  Pepita (C₀)         : %s\n", .esp_fmt_num(C0, 4)))
  cat(sprintf("  Contribuição (C)    : %s\n", .esp_fmt_num(C, 4)))
  cat(sprintf("  Patamar (C₀ + C)    : %s\n", .esp_fmt_num(sill, 4)))
  cat(sprintf("  Alcance Prático (a) : %s\n", .esp_fmt_num(a, 4)))
  cat(sprintf("  SQR do Ajuste       : %s\n", .esp_fmt_num(sqr_val, 4)))
  cat(.esp_sep(w_box), "\n\n")

  .print_topico("DIAGNÓSTICO DO VARIOGRAMA")
  cat(sprintf("  Razão Pepita/Patamar: %s%%\n", .esp_fmt_num(razao_pepita, 1)))
  if (razao_pepita <= 25) {
    cat("  [✓] Dependência espacial Forte (razão ≤ 25%): estrutura espacial bem definida.\n")
  } else if (razao_pepita <= 75) {
    cat("  [✓] Dependência espacial Moderada (25% < razão ≤ 75%): estrutura espacial consistente.\n")
  } else {
    cat("  [!] Dependência espacial Fraca ou Efeito Pepita (razão > 75%): alta variabilidade local.\n")
  }
  cat(sprintf("  O modelo %s indica que a correlação espacial dos dados se esgota a ~%s unidades de distância.\n\n",
              nome_mod_pt, .esp_fmt_num(a, 2)))
  .print_rodape()

  # ── Gráficos em ggplot2 com meu_tema(grade = "dupla") ───────────────────────
  p_plot <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_emp <- data.frame(Distancia = lag_h, Semivariancia = gamma_emp, Pares = n_pares)

    h_curva <- seq(0, dist_max_val, length.out = 150)
    g_curva <- if (tipo_selecionado == "esferico") {
      f_esferico(h_curva, C0, C, a)
    } else if (tipo_selecionado == "exponencial") {
      f_exponencial(h_curva, C0, C, a)
    } else {
      f_gaussiano(h_curva, C0, C, a)
    }
    df_teo <- data.frame(Distancia = h_curva, Semivariancia = g_curva)

    p1 <- ggplot2::ggplot(df_emp, ggplot2::aes(x = Distancia, y = Semivariancia)) +
      ggplot2::geom_hline(yintercept = sill, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = C0, linetype = "dotted", color = "gray40", linewidth = 0.7) +
      ggplot2::geom_line(data = df_teo, color = "#1D3557", linewidth = 1.1) +
      ggplot2::geom_point(color = "#457B9D", size = 2.8, alpha = 0.9) +
      ggplot2::labs(
        title    = sprintf("Variograma: %s (%s)", var_nome, nome_mod_pt),
        subtitle = sprintf("Pepita: %s | Patamar: %s | Alcance: %s",
                           .esp_fmt_num(C0, 2), .esp_fmt_num(sill, 2), .esp_fmt_num(a, 2)),
        x        = "Distância (h)",
        y        = "Semivariância γ(h)"
      ) +
      meu_tema(grade = "dupla")

    p2 <- ggplot2::ggplot(df_emp, ggplot2::aes(x = Distancia, y = Pares)) +
      ggplot2::geom_col(fill = "#457B9D", color = "white", alpha = 0.85) +
      ggplot2::labs(
        title    = "Número de Pares de Pontos por Classe",
        subtitle = sprintf("Total: %s pares analisados", .esp_fmt_num(sum(n_pares), 0)),
        x        = "Distância (h)",
        y        = "Contagem de Pares"
      ) +
      meu_tema(grade = "dupla")

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      p_plot <- list(painel_variograma = p1, painel_pares = p2)
    } else {
      print(p1)
      p_plot <- p1
    }
  }

  res <- list(
    gamma_empirico = data.frame(lag = lag_h, gamma = gamma_emp, n_pares = n_pares),
    modelo         = tipo_selecionado,
    modelo_nome    = nome_mod_pt,
    pepita         = C0,
    nugget         = C0,
    contribuicao   = C,
    patamar        = sill,
    sill           = sill,
    alcance        = a,
    range          = a,
    sqr            = sqr_val,
    razao_pepita   = razao_pepita,
    variavel       = var_nome,
    coordenadas    = coords,
    z              = z,
    dist_max       = dist_max_val,
    n              = n,
    fn_teorica     = function(h) {
      if (tipo_selecionado == "esferico")   f_esferico(h, C0, C, a)
      else if (tipo_selecionado == "exponencial") f_exponencial(h, C0, C, a)
      else f_gaussiano(h, C0, C, a)
    },
    grafico        = p_plot
  )
  class(res) <- c("estatR_variograma", "list")
  invisible(res)
}


# ── 12. KRIGAGEM ORDINÁRIA (INTERPOLAÇÃO ESPACIAL) ───────────────────────────

#' @title Interpolação Espacial por Krigagem Ordinária
#' @description Realiza interpolação espacial ótima de uma variável contínua em uma grade
#'   regular através de Krigagem Ordinária (Best Linear Unbiased Estimator - BLUE),
#'   gerando estimativas pontuais e mapas de incerteza (variância de krigagem).
#'
#' @param dados Data frame, matriz ou objeto \code{sf} contendo coordenadas e valores amostrados.
#' @param variavel Nome da coluna da variável de resposta.
#' @param coordenadas Vetor com nomes ou índices das colunas de coordenadas (ex: \code{c("x", "y")}).
#' @param variograma_obj Objeto retornado por \code{\link{variograma}}. Se \code{NULL},
#'   o variograma é ajustado automaticamente.
#' @param resolucao Resolução da malha regular: \code{"auto"} (grade 40×40 = 1.600 pontos)
#'   ou um inteiro definindo o número de divisões em cada eixo (ex: 30, 50).
#' @param janela Limites espaciais \code{c(xmin, xmax, ymin, ymax)} ou \code{"auto"} (envoltória amostrada + 5% padding).
#' @param grafico Lógico. Se TRUE (padrão), exibe mapa da predição e mapa de variância via \code{meu_tema(grade = "dupla")}.
#'
#' @return Objeto invisível da classe \code{estatR_krigagem}.
#' @export
krigagem <- function(dados,
                     variavel = NULL,
                     coordenadas = NULL,
                     variograma_obj = NULL,
                     resolucao = "auto",
                     janela = "auto",
                     grafico = TRUE) {

  # 1. Ajuste do variograma se não fornecido
  if (is.null(variograma_obj)) {
    variograma_obj <- variograma(dados, variavel = variavel, coordenadas = coordenadas, grafico = FALSE)
  }
  if (!inherits(variograma_obj, "estatR_variograma")) {
    stop("O argumento 'variograma_obj' deve ser um objeto gerado por 'variograma()'.", call. = FALSE)
  }

  coords <- variograma_obj$coordenadas
  z <- variograma_obj$z
  var_nome <- variograma_obj$variavel
  n <- length(z)
  fn_semivar <- variograma_obj$fn_teorica

  # 2. Definição da janela de interpolação
  if (is.character(janela) && janela == "auto") {
    rx <- range(coords$x, na.rm = TRUE)
    ry <- range(coords$y, na.rm = TRUE)
    dx <- rx[2] - rx[1]
    dy <- ry[2] - ry[1]
    pad_x <- if (dx > 0) dx * 0.05 else 0.5
    pad_y <- if (dy > 0) dy * 0.05 else 0.5
    janela_vec <- c(xmin = rx[1] - pad_x, xmax = rx[2] + pad_x,
                    ymin = ry[1] - pad_y, ymax = ry[2] + pad_y)
  } else if (is.numeric(janela) && length(janela) == 4) {
    janela_vec <- janela
    names(janela_vec) <- c("xmin", "xmax", "ymin", "ymax")
  } else {
    stop("O argumento 'janela' deve ser 'auto' ou um vetor c(xmin, xmax, ymin, ymax).", call. = FALSE)
  }

  # 3. Resolução da grade
  res_val <- if (identical(resolucao, "auto")) 40L else as.integer(resolucao)
  if (res_val < 5) res_val <- 5L

  gx <- seq(janela_vec["xmin"], janela_vec["xmax"], length.out = res_val)
  gy <- seq(janela_vec["ymin"], janela_vec["ymax"], length.out = res_val)
  grade <- expand.grid(x = gx, y = gy)
  n_grade <- nrow(grade)

  # 4. Construção do sistema de Krigagem Ordinária
  # Matriz Gamma (n+1) x (n+1)
  coords_mat <- as.matrix(coords[, c("x", "y")])
  D_amostral <- as.matrix(stats::dist(coords_mat))
  Gamma_mat <- fn_semivar(D_amostral)
  diag(Gamma_mat) <- 0

  # Adiciona pequena perturbação na diagonal para estabilidade numérica
  diag(Gamma_mat) <- diag(Gamma_mat) + 1e-8

  A <- matrix(1, nrow = n + 1, ncol = n + 1)
  A[1:n, 1:n] <- Gamma_mat
  A[n + 1, n + 1] <- 0

  # Inversão única do sistema A (n+1) x (n+1)
  A_inv <- tryCatch({
    solve(A)
  }, error = function(e) {
    MASS::ginv(A)
  })

  # 5. Predição vetorizada na grade
  grade_mat <- as.matrix(grade[, c("x", "y")])

  # Matriz de distâncias entre pontos amostrais (n) e pontos da grade (n_grade)
  # dx: grade_mat[,1] vs coords_mat[,1], dy: grade_mat[,2] vs coords_mat[,2]
  # Para economizar memória e manter performance alta, calculamos em blocos se necessário
  bloco_tam <- 2000L
  n_blocos <- ceiling(n_grade / bloco_tam)

  z_pred_vec <- numeric(n_grade)
  var_krig_vec <- numeric(n_grade)

  for (b in seq_len(n_blocos)) {
    idx_b <- seq((b - 1) * bloco_tam + 1, min(b * bloco_tam, n_grade))
    sub_grade <- grade_mat[idx_b, , drop = FALSE]
    nb <- length(idx_b)

    # dist_mat_b: dimensão n x nb
    # dist(sub_grade, coords_mat)
    # dist_ij = sqrt((x_g - x_a)^2 + (y_g - y_a)^2)
    dx_b <- outer(coords_mat[, 1], sub_grade[, 1], "-")
    dy_b <- outer(coords_mat[, 2], sub_grade[, 2], "-")
    D_b <- sqrt(dx_b^2 + dy_b^2)

    # Vetor de semivariâncias b0: dimensão (n+1) x nb
    B0 <- matrix(1, nrow = n + 1, ncol = nb)
    B0[1:n, ] <- fn_semivar(D_b)

    # Pesos W: dimensão (n+1) x nb = A_inv %*% B0
    W_b <- A_inv %*% B0

    # Predição: soma dos pesos 1..n multiplicados por z
    z_pred_vec[idx_b] <- as.numeric(crossprod(W_b[1:n, , drop = FALSE], z))

    # Variância de Krigagem: sum(w_i * gamma_i0) + mu = colSums(W_b * B0)
    var_krig_vec[idx_b] <- colSums(W_b * B0)
  }

  var_krig_vec[var_krig_vec < 0] <- 0
  grade$predicao <- z_pred_vec
  grade$var_krigagem <- var_krig_vec
  grade$dp_krigagem <- sqrt(var_krig_vec)

  # ── Saída no console ────────────────────────────────────────────────────────
  .print_titulo("KRIGAGEM ORDINÁRIA (INTERPOLAÇÃO ESPACIAL)")
  cat(sprintf("  Variável            : %s\n", var_nome))
  cat(sprintf("  Nº Observações      : %d pontos amostrados\n", n))
  cat(sprintf("  Modelo de Variograma: %s (pepita = %s | patamar = %s | alcance = %s)\n",
              variograma_obj$modelo_nome,
              .esp_fmt_num(variograma_obj$pepita, 2),
              .esp_fmt_num(variograma_obj$patamar, 2),
              .esp_fmt_num(variograma_obj$alcance, 2)))
  modo_res <- if (identical(resolucao, "auto")) " (selecionada automaticamente)" else ""
  cat(sprintf("  Grade de Predição   : %d × %d (%d locais)%s\n", res_val, res_val, n_grade, modo_res))
  cat(sprintf("  Janela              : X ∈ [%s, %s] | Y ∈ [%s, %s]\n\n",
              .esp_fmt_num(janela_vec["xmin"], 2), .esp_fmt_num(janela_vec["xmax"], 2),
              .esp_fmt_num(janela_vec["ymin"], 2), .esp_fmt_num(janela_vec["ymax"], 2)))

  .print_topico("DISTRIBUIÇÃO DA SUPERFÍCIE PREDITA")
  w_c <- c(est = 18, pred = 18, vkrig = 20)
  cat(.esp_sep(sum(w_c)), "\n")
  cat("  ",
      .esp_pad("Estatística", w_c["est"]),
      .esp_pad("Predição", w_c["pred"], "right"),
      .esp_pad("Variância Krigagem", w_c["vkrig"], "right"),
      "\n", sep = "")
  cat(.esp_sep(sum(w_c)), "\n")

  q_p <- stats::quantile(grade$predicao, probs = c(0, 0.25, 0.5, 0.75, 1))
  q_v <- stats::quantile(grade$var_krigagem, probs = c(0, 0.25, 0.5, 0.75, 1))

  imprime_linha_tab <- function(rotulo, val_p, val_v) {
    cat("  ",
        .esp_pad(rotulo, w_c["est"]),
        .esp_pad(.esp_fmt_num(val_p, 2), w_c["pred"], "right"),
        .esp_pad(.esp_fmt_num(val_v, 2), w_c["vkrig"], "right"),
        "\n", sep = "")
  }

  imprime_linha_tab("Mínimo", q_p[1], q_v[1])
  imprime_linha_tab("1º Quartil", q_p[2], q_v[2])
  imprime_linha_tab("Mediana", q_p[3], q_v[3])
  imprime_linha_tab("Média", mean(grade$predicao), mean(grade$var_krigagem))
  imprime_linha_tab("3º Quartil", q_p[4], q_v[4])
  imprime_linha_tab("Máximo", q_p[5], q_v[5])
  cat(.esp_sep(sum(w_c)), "\n\n")

  cat(sprintf("  Variância média de Krigagem: %s (incerteza média quadrática de predição)\n\n",
              .esp_fmt_num(mean(grade$var_krigagem), 2)))
  .print_rodape()

  # ── Gráficos em ggplot2 com meu_tema(grade = "dupla") ───────────────────────
  p_plot <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_pts <- as.data.frame(coords)
    df_pts$z <- z

    p1 <- ggplot2::ggplot(grade, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = predicao)) +
      ggplot2::scale_fill_gradientn(
        colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"),
        name = var_nome
      ) +
      ggplot2::geom_point(data = df_pts, ggplot2::aes(x = x, y = y),
                          color = "black", fill = "white", shape = 21, size = 1.8, stroke = 0.5) +
      ggplot2::labs(
        title    = sprintf("Superfície Krigada: %s", var_nome),
        subtitle = sprintf("Interpolação BLUE | Grade %d×%d | Pontos: Amostras", res_val, res_val),
        x        = "Coordenada X",
        y        = "Coordenada Y"
      ) +
      meu_tema(grade = "dupla")

    p2 <- ggplot2::ggplot(grade, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = var_krigagem)) +
      ggplot2::scale_fill_gradientn(
        colours = c("#FFFFCC", "#FED976", "#FD8D3C", "#FC4E2A", "#E31A1C", "#800026"),
        name = "Variância"
      ) +
      ggplot2::geom_point(data = df_pts, ggplot2::aes(x = x, y = y),
                          color = "black", shape = 3, size = 1.5, alpha = 0.6) +
      ggplot2::labs(
        title    = "Mapa de Incerteza (Variância de Krigagem)",
        subtitle = "Menores valores (amarelo) indicam maior proximidade das amostras",
        x        = "Coordenada X",
        y        = "Coordenada Y"
      ) +
      meu_tema(grade = "dupla")

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      p_plot <- list(painel_predicao = p1, painel_variancia = p2)
    } else {
      print(p1)
      p_plot <- p1
    }
  }

  res <- list(
    grade          = grade,
    variograma     = variograma_obj,
    resolucao      = c(x = res_val, y = res_val),
    janela         = janela_vec,
    variavel       = var_nome,
    dados_obs      = data.frame(x = coords$x, y = coords$y, z = z),
    n              = n,
    valores_ajustados = z_pred_vec,
    grafico        = p_plot
  )
  class(res) <- c("estatR_krigagem", "list")
  invisible(res)
}


# ── 13. VALIDAÇÃO CRUZADA DA KRIGAGEM (LEAVE-ONE-OUT) ─────────────────────────

#' @title Validação Cruzada Leave-One-Out para Krigagem Ordinária
#' @description Avalia a qualidade preditiva e calibração do modelo de interpolação
#'   espacial via validação cruzada do tipo Leave-One-Out (LOO), calculando resíduos,
#'   erros padronizados, RMSE, MAE, MAPE e R².
#'
#' @param dados Data frame, matriz ou objeto \code{sf} contendo coordenadas e valores amostrados.
#' @param variavel Nome da coluna com a variável de resposta.
#' @param coordenadas Vetor de nomes ou índices das colunas de coordenadas (ex: \code{c("x", "y")}).
#' @param variograma_obj Objeto retornado por \code{\link{variograma}}. Se \code{NULL},
#'   o variograma é ajustado automaticamente.
#' @param grafico Lógico. Se TRUE (padrão), exibe gráfico de dispersão observado vs. predito
#'   e histograma dos erros padronizados com \code{meu_tema(grade = "dupla")}.
#'
#' @return Objeto invisível da classe \code{estatR_validacao_krigagem}.
#' @export
validar_krigagem <- function(dados,
                             variavel = NULL,
                             coordenadas = NULL,
                             variograma_obj = NULL,
                             grafico = TRUE) {

  if (is.null(variograma_obj)) {
    variograma_obj <- variograma(dados, variavel = variavel, coordenadas = coordenadas, grafico = FALSE)
  }
  if (!inherits(variograma_obj, "estatR_variograma")) {
    stop("O argumento 'variograma_obj' deve ser um objeto gerado por 'variograma()'.", call. = FALSE)
  }

  coords <- variograma_obj$coordenadas
  z <- variograma_obj$z
  var_nome <- variograma_obj$variavel
  n <- length(z)
  fn_semivar <- variograma_obj$fn_teorica

  coords_mat <- as.matrix(coords[, c("x", "y")])
  D_total <- as.matrix(stats::dist(coords_mat))

  z_pred_loo <- numeric(n)
  var_krig_loo <- numeric(n)

  # Leave-One-Out Loop
  for (i in seq_len(n)) {
    idx_train <- seq_len(n)[-i]
    m <- length(idx_train)

    coords_train <- coords_mat[idx_train, , drop = FALSE]
    z_train <- z[idx_train]

    # Matriz A_train
    D_train <- D_total[idx_train, idx_train, drop = FALSE]
    G_train <- fn_semivar(D_train)
    diag(G_train) <- 1e-8

    A_tr <- matrix(1, nrow = m + 1, ncol = m + 1)
    A_tr[1:m, 1:m] <- G_train
    A_tr[m + 1, m + 1] <- 0

    A_tr_inv <- tryCatch(solve(A_tr), error = function(e) MASS::ginv(A_tr))

    # Vetor b0 para o ponto i
    d_i <- D_total[idx_train, i]
    b0 <- c(fn_semivar(d_i), 1)

    w_i <- A_tr_inv %*% b0
    z_pred_loo[i] <- sum(w_i[1:m] * z_train)
    v_i <- sum(w_i * b0)
    var_krig_loo[i] <- max(1e-6, v_i)
  }

  erros <- z - z_pred_loo
  dp_krig <- sqrt(var_krig_loo)
  erros_pad <- erros / dp_krig

  # Métricas
  rmse <- sqrt(mean(erros^2, na.rm = TRUE))
  mae  <- mean(abs(erros), na.rm = TRUE)
  mape <- if (all(abs(z) > 1e-5)) mean(abs(erros / z), na.rm = TRUE) * 100 else NA_real_

  ss_res <- sum(erros^2, na.rm = TRUE)
  ss_tot <- sum((z - mean(z, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- if (ss_tot > 0) max(0, 1 - (ss_res / ss_tot)) else 0

  media_ep <- mean(erros_pad, na.rm = TRUE)
  dp_ep <- stats::sd(erros_pad, na.rm = TRUE)

  # ── Saída no console ────────────────────────────────────────────────────────
  .print_titulo("VALIDAÇÃO CRUZADA (LEAVE-ONE-OUT): KRIGAGEM ORDINÁRIA")
  cat(sprintf("  Variável            : %s\n", var_nome))
  cat(sprintf("  Nº Observações      : %d pontos avaliados\n", n))
  cat(sprintf("  Modelo de Variograma: %s (pepita = %s | patamar = %s | alcance = %s)\n\n",
              variograma_obj$modelo_nome,
              .esp_fmt_num(variograma_obj$pepita, 2),
              .esp_fmt_num(variograma_obj$patamar, 2),
              .esp_fmt_num(variograma_obj$alcance, 2)))

  .print_topico("MÉTRICAS DE DESEMPENHO PREDITIVO (LOO)")
  w_m <- c(met = 12, val = 14, desc = 42)
  cat(.esp_sep(sum(w_m)), "\n")
  cat("  ",
      .esp_pad("Métrica", w_m["met"]),
      .esp_pad("Valor", w_m["val"], "right"),
      "  ", .esp_pad("Interpretação", w_m["desc"]),
      "\n", sep = "")
  cat(.esp_sep(sum(w_m)), "\n")

  cat("  ", .esp_pad("RMSE", w_m["met"]), .esp_pad(.esp_fmt_num(rmse, 4), w_m["val"], "right"),
      "    ", .esp_pad("Raiz do Erro Quadrático Médio", w_m["desc"]), "\n", sep = "")
  cat("  ", .esp_pad("MAE", w_m["met"]), .esp_pad(.esp_fmt_num(mae, 4), w_m["val"], "right"),
      "    ", .esp_pad("Erro Absoluto Médio", w_m["desc"]), "\n", sep = "")
  if (!is.na(mape)) {
    cat("  ", .esp_pad("MAPE (%)", w_m["met"]), .esp_pad(paste0(.esp_fmt_num(mape, 2), "%"), w_m["val"], "right"),
        "    ", .esp_pad("Erro Percentual Absoluto Médio", w_m["desc"]), "\n", sep = "")
  }
  cat("  ", .esp_pad("R² (LOO)", w_m["met"]), .esp_pad(.esp_fmt_num(r2, 4), w_m["val"], "right"),
      "    ", .esp_pad("Coeficiente de Determinação (Predito vs Obs)", w_m["desc"]), "\n", sep = "")
  cat(.esp_sep(sum(w_m)), "\n\n")

  .print_topico("DIAGNÓSTICO DOS ERROS PADRONIZADOS")
  cat(sprintf("  Média dos Erros Padronizados : %s  (ideal: próxima de 0)\n", .esp_fmt_num(media_ep, 4)))
  cat(sprintf("  Desvio Padrão dos Erros Pad. : %s  (ideal: próximo de 1)\n", .esp_fmt_num(dp_ep, 4)))

  if (abs(media_ep) <= 0.20 && dp_ep >= 0.80 && dp_ep <= 1.20) {
    cat("  [✓] Variograma bem calibrado: média ≈ 0 e desvio ≈ 1 confirmam acurácia e incerteza realistas.\n\n")
  } else if (dp_ep > 1.20) {
    cat("  [!] Variograma subestimou a incerteza real (desvio padrão > 1,20): o patamar pode estar subestimado.\n\n")
  } else if (dp_ep < 0.80) {
    cat("  [!] Variograma superestimou a incerteza real (desvio padrão < 0,80): o patamar pode estar superestimado.\n\n")
  } else {
    cat("  [!] Erros padronizados apresentam leve viés médio na predição pontual.\n\n")
  }
  .print_rodape()

  # ── Gráficos em ggplot2 com meu_tema(grade = "dupla") ───────────────────────
  p_plot <- NULL
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    df_loo <- data.frame(
      Observado = z,
      Predito = z_pred_loo,
      Erro = erros,
      Erro_Pad = erros_pad
    )

    lims <- range(c(z, z_pred_loo), na.rm = TRUE)

    p1 <- ggplot2::ggplot(df_loo, ggplot2::aes(x = Observado, y = Predito)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#E63946", linewidth = 0.9) +
      ggplot2::geom_point(color = "#1D3557", fill = "#457B9D", shape = 21, size = 2.4, alpha = 0.85) +
      ggplot2::labs(
        title    = "Observado vs. Predito (Validação Cruzada LOO)",
        subtitle = sprintf("R² = %s | RMSE = %s | Linha tracejada: y = x", .esp_fmt_num(r2, 3), .esp_fmt_num(rmse, 2)),
        x        = "Valor Observado",
        y        = "Valor Predito"
      ) +
      ggplot2::coord_fixed(ratio = 1, xlim = lims, ylim = lims) +
      meu_tema(grade = "dupla")

    p2 <- ggplot2::ggplot(df_loo, ggplot2::aes(x = Erro_Pad)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                              bins = max(7, min(20, round(n / 4))),
                              fill = "#457B9D", color = "white", alpha = 0.85) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
      ggplot2::stat_function(fun = stats::dnorm, color = "#1D3557", linewidth = 1.1) +
      ggplot2::labs(
        title    = "Distribuição dos Erros Padronizados",
        subtitle = sprintf("Média = %s | DP = %s | Curva azul: N(0, 1)", .esp_fmt_num(media_ep, 2), .esp_fmt_num(dp_ep, 2)),
        x        = "Erro Padronizado",
        y        = "Densidade"
      ) +
      meu_tema(grade = "dupla")

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p1, p2, ncol = 2)
      p_plot <- list(painel_obs_pred = p1, painel_erros_pad = p2)
    } else {
      print(p1)
      p_plot <- p1
    }
  }

  res <- list(
    observado       = z,
    predito         = z_pred_loo,
    erro            = erros,
    erro_padronizado= erros_pad,
    variancia_krig  = var_krig_loo,
    rmse            = rmse,
    mae             = mae,
    mape            = mape,
    r2              = r2,
    media_erro_pad  = media_ep,
    dp_erro_pad     = dp_ep,
    variograma      = variograma_obj,
    n               = n,
    grafico         = p_plot
  )
  class(res) <- c("estatR_validacao_krigagem", "list")
  invisible(res)
}


# ── 14. POLIMORFISMO PARA SUPERFÍCIE CONTÍNUA ─────────────────────────────────

#' @noRd
.analise_residual_krigagem <- function(modelo, grafico = TRUE) {
  validar_krigagem(
    dados          = modelo$dados_obs,
    variavel       = "z",
    coordenadas    = c("x", "y"),
    variograma_obj = modelo$variograma,
    grafico        = grafico
  )
}

#' @noRd
.metricas_krigagem <- function(modelo) {
  val <- validar_krigagem(
    dados          = modelo$dados_obs,
    variavel       = "z",
    coordenadas    = c("x", "y"),
    variograma_obj = modelo$variograma,
    grafico        = FALSE
  )
  invisible(data.frame(
    Metrica = c("RMSE", "MAE", "MAPE", "R2_LOO", "Media_Erro_Pad", "DP_Erro_Pad"),
    Valor   = c(val$rmse, val$mae, val$mape, val$r2, val$media_erro_pad, val$dp_erro_pad),
    stringsAsFactors = FALSE
  ))
}

#' @noRd
.metricas_variograma <- function(modelo) {
  .print_titulo("MÉTRICAS DO VARIOGRAMA AJUSTADO")
  cat(sprintf("  Variável            : %s\n", modelo$variavel))
  cat(sprintf("  Modelo              : %s\n", modelo$modelo_nome))
  cat(sprintf("  Pepita (C₀)         : %s\n", .esp_fmt_num(modelo$pepita, 4)))
  cat(sprintf("  Patamar (C₀ + C)    : %s\n", .esp_fmt_num(modelo$patamar, 4)))
  cat(sprintf("  Contribuição (C)    : %s\n", .esp_fmt_num(modelo$contribuicao, 4)))
  cat(sprintf("  Alcance Prático (a) : %s\n", .esp_fmt_num(modelo$alcance, 4)))
  cat(sprintf("  Razão Pepita/Patamar: %s%%\n", .esp_fmt_num(modelo$razao_pepita, 2)))
  cat(sprintf("  SQR do Ajuste       : %s\n\n", .esp_fmt_num(modelo$sqr, 6)))
  .print_rodape()

  invisible(data.frame(
    Parametro = c("Pepita", "Patamar", "Contribuicao", "Alcance", "Razao_Pepita_Pct", "SQR"),
    Valor     = c(modelo$pepita, modelo$patamar, modelo$contribuicao, modelo$alcance, modelo$razao_pepita, modelo$sqr),
    stringsAsFactors = FALSE
  ))
}

#' @noRd
.comparar_modelos_variograma <- function(modelos, nomes = NULL) {
  n_mods <- length(modelos)
  if (is.null(nomes)) nomes <- paste0("Modelo_", seq_len(n_mods))

  tipos <- vapply(modelos, function(m) m$modelo_nome, character(1))
  vars  <- vapply(modelos, function(m) m$variavel, character(1))
  peps  <- vapply(modelos, function(m) m$pepita, numeric(1))
  sills <- vapply(modelos, function(m) m$patamar, numeric(1))
  alcs  <- vapply(modelos, function(m) m$alcance, numeric(1))
  sqrs  <- vapply(modelos, function(m) m$sqr, numeric(1))

  df_comp <- data.frame(
    Nome     = nomes,
    Modelo   = tipos,
    Variavel = vars,
    Pepita   = peps,
    Patamar  = sills,
    Alcance  = alcs,
    SQR      = sqrs,
    stringsAsFactors = FALSE
  )
  df_comp <- df_comp[order(df_comp$SQR), ]

  .print_titulo("COMPARAÇÃO DE MODELOS DE VARIOGRAMA")
  w <- c(nome = 14, mod = 14, pep = 10, pat = 10, alc = 10, sqr = 12)
  cat(.esp_sep(sum(w)), "\n")
  cat("  ",
      .esp_pad("Nome", w["nome"]),
      .esp_pad("Modelo", w["mod"]),
      .esp_pad("Pepita", w["pep"], "right"),
      .esp_pad("Patamar", w["pat"], "right"),
      .esp_pad("Alcance", w["alc"], "right"),
      .esp_pad("SQR", w["sqr"], "right"),
      "\n", sep = "")
  cat(.esp_sep(sum(w)), "\n")

  for (i in seq_len(nrow(df_comp))) {
    cat("  ",
        .esp_pad(df_comp$Nome[i], w["nome"]),
        .esp_pad(df_comp$Modelo[i], w["mod"]),
        .esp_pad(.esp_fmt_num(df_comp$Pepita[i], 2), w["pep"], "right"),
        .esp_pad(.esp_fmt_num(df_comp$Patamar[i], 2), w["pat"], "right"),
        .esp_pad(.esp_fmt_num(df_comp$Alcance[i], 2), w["alc"], "right"),
        .esp_pad(.esp_fmt_num(df_comp$SQR[i], 4), w["sqr"], "right"),
        "\n", sep = "")
  }
  cat(.esp_sep(sum(w)), "\n")
  cat("  Notas:\n")
  cat("  Modelos ordenados pelo menor SQR (melhor ajuste ao variograma empírico).\n\n")
  .print_rodape()

  invisible(df_comp)
}
