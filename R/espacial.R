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
