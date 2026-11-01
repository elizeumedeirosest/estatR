# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: PROCESSOS ESTOCÁSTICOS
# Cadeias de Markov, Passeio Aleatório e Processos Pontuais de Poisson
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. CADEIAS DE MARKOV EM TEMPO DISCRETO ───────────────────────────────────

#' @title Análise e Simulação de Cadeias de Markov
#' @description Realiza a análise teórica e a simulação de trajetórias de uma Cadeia
#'   de Markov em Tempo Discreto (DTMC). Calcula a distribuição estacionária (se ergódica),
#'   classifica os estados (absorventes, transientes, recorrentes) e simula o comportamento futuro.
#'
#' @param matriz_transicao Matriz quadrada de probabilidades de transição \eqn{P},
#'   onde cada linha deve somar 1.
#' @param estados Vetor de caracteres com os nomes dos estados. Se \code{NULL}, usa E1, E2...
#' @param estado_inicial Nome ou índice do estado inicial para a simulação.
#' @param passos Número de passos temporais a projetar / simular (padrão: 10).
#' @param n_simulacoes Número de trajetórias a simular (padrão: 1).
#' @param semente Semente para reprodutibilidade da simulação.
#' @param grafico Lógico. Se TRUE (padrão), plota a evolução das probabilidades dos estados.
#'
#' @return Objeto invisível da classe \code{estatR_markov}.
#' @export
cadeia_markov <- function(matriz_transicao,
                          estados = NULL,
                          estado_inicial = 1,
                          passos = 10,
                          n_simulacoes = 1,
                          semente = NULL,
                          grafico = TRUE) {

  if (!is.null(semente)) set.seed(semente)

  P <- as.matrix(matriz_transicao)
  if (nrow(P) != ncol(P)) {
    stop("A matriz de transi\u00e7\u00e3o deve ser estritamente quadrada.", call. = FALSE)
  }

  k <- nrow(P)
  if (k < 2) {
    stop("A cadeia deve conter pelo menos 2 estados.", call. = FALSE)
  }

  if (any(P < 0) || any(P > 1)) {
    stop("Todas as probabilidades de transi\u00e7\u00e3o devem estar no intervalo [0, 1].", call. = FALSE)
  }

  somas_linhas <- rowSums(P)
  if (any(abs(somas_linhas - 1) > 1e-4)) {
    stop("Cada linha da matriz de transi\u00e7\u00e3o deve somar exatamente 1.", call. = FALSE)
  }

  # Nomes dos estados
  if (is.null(estados)) {
    if (!is.null(rownames(P))) {
      estados <- rownames(P)
    } else {
      estados <- paste0("E", seq_len(k))
    }
  } else if (length(estados) != k) {
    stop("O comprimento do vetor de 'estados' deve coincidir com a dimens\u00e3o da matriz.", call. = FALSE)
  }
  rownames(P) <- colnames(P) <- estados

  # Identificação de estados absorventes
  absorventes <- which(diag(P) == 1)
  nomes_abs   <- estados[absorventes]

  # ── Cálculo da Distribuição Estacionária Pi (pi * P = pi) ────────────────────
  pi_estacionario <- rep(NA_real_, k)
  tryCatch({
    A <- t(P) - diag(k)
    A[k, ] <- 1
    b <- c(rep(0, k - 1), 1)
    sol <- solve(A, b)
    if (all(sol >= -1e-5) && abs(sum(sol) - 1) < 1e-3) {
      pi_estacionario <- pmax(0, sol)
      pi_estacionario <- pi_estacionario / sum(pi_estacionario)
    }
  }, error = function(e) NULL)

  # Projeção analítica de P^n
  Pn <- diag(k)
  Pn_lista <- list(Pn)
  for (t in seq_len(passos)) {
    Pn <- Pn %*% P
    Pn_lista[[t + 1]] <- Pn
  }

  # Vetor de estado inicial
  v_inicial <- rep(0, k)
  if (is.character(estado_inicial)) {
    idx_init <- match(estado_inicial, estados)
    if (is.na(idx_init)) stop(sprintf("Estado inicial '%s' n\u00e3o encontrado.", estado_inicial), call. = FALSE)
  } else {
    idx_init <- as.integer(estado_inicial)
    if (idx_init < 1 || idx_init > k) stop("Índice do estado inicial fora dos limites.", call. = FALSE)
  }
  v_inicial[idx_init] <- 1

  # Evolução temporal das probabilidades a partir do estado inicial
  df_evolucao_list <- list()
  v_t <- v_inicial
  for (t in 0:passos) {
    df_evolucao_list[[t + 1]] <- data.frame(
      Tempo       = t,
      Estado      = estados,
      Probabilidade = as.numeric(v_t),
      stringsAsFactors = FALSE
    )
    v_t <- as.numeric(v_t %*% P)
  }
  df_evolucao <- do.call(rbind, df_evolucao_list)

  # ── Simulação de Trajetórias ────────────────────────────────────────────────
  trajetorias <- matrix(0L, nrow = n_simulacoes, ncol = passos + 1)
  for (s in seq_len(n_simulacoes)) {
    traj <- integer(passos + 1)
    traj[1] <- idx_init
    atual <- idx_init
    for (t in seq_len(passos)) {
      atual <- sample(seq_len(k), size = 1, prob = P[atual, ])
      traj[t + 1] <- atual
    }
    trajetorias[s, ] <- traj
  }

  # ── Gráfico de Evolução ─────────────────────────────────────────────────────
  p_markov <- NULL
  if (isTRUE(grafico)) {
    p_markov <- ggplot2::ggplot(df_evolucao, ggplot2::aes(x = Tempo, y = Probabilidade, color = Estado, group = Estado)) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::scale_x_continuous(breaks = 0:passos) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1)) +
      ggplot2::labs(
        title    = "Evolu\u00e7\u00e3o das Probabilidades de Estado: Cadeia de Markov",
        subtitle = sprintf("Estado Inicial: %s  |  Proje\u00e7\u00e3o para %d passos", estados[idx_init], passos),
        x        = "Passo de Tempo (n)",
        y        = "Probabilidade de Ocupa\u00e7\u00e3o P(X_n = j)"
      )

    if (exists("meu_tema", mode = "function")) {
      p_markov <- p_markov + meu_tema()
    } else {
      p_markov <- p_markov + ggplot2::theme_minimal()
    }
    print(p_markov)
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

  .print_titulo("PROCESSOS ESTOC\u00c1STICOS: CADEIA DE MARKOV")

  cat(sprintf("  Espa\u00e7o de Estados    : %d estados (%s)\n", k, paste(estados, collapse = ", ")))
  cat(sprintf("  Estado Inicial       : %s\n", estados[idx_init]))
  cat(sprintf("  Horizonte de Proje\u00e7\u00e3o: %d passos\n", passos))
  if (length(nomes_abs) > 0) {
    cat(sprintf("  Estados Absorventes  : %s\n\n", paste(nomes_abs, collapse = ", ")))
  } else {
    cat("  Estrutura da Cadeia  : Irredut\u00edvel / N\u00e3o absorvente\n\n")
  }

  .print_topico("MATRIZ DE TRANSI\u00c7\u00c3O DE 1 PASSO (P)")
  w_m <- max(10, max(nchar(estados)) + 3)
  hdr_p <- paste0(pad("", w_m, "left"), paste(sapply(estados, function(e) pad(e, w_m, "center")), collapse = ""))
  w_tot_p <- nchar(hdr_p)

  cat(sep_l(w_tot_p), "\n")
  cat("  ", hdr_p, "\n", sep = "")
  cat(sep_l(w_tot_p), "\n")
  for (i in seq_len(k)) {
    linha_txt <- paste(sapply(P[i, ], function(v) pad(.fmt_pct(v, 1), w_m, "center")), collapse = "")
    cat("  ", pad(estados[i], w_m, "left"), linha_txt, "\n", sep = "")
  }
  cat(sep_l(w_tot_p), "\n\n")

  if (all(!is.na(pi_estacionario))) {
    .print_topico("DISTRIBUI\u00c7\u00c3O ESTACION\u00c1RIA DE LONGO PRAZO (\u03c0)")
    w_pi <- c(est = 16, prob = 16)
    hdr_pi <- paste0(pad("Estado", w_pi["est"], "left"), pad("Probabilidade \u03c0", w_pi["prob"], "center"))
    cat(sep_l(sum(w_pi) + 2), "\n")
    cat("  ", hdr_pi, "\n", sep = "")
    cat(sep_l(sum(w_pi) + 2), "\n")
    for (i in seq_len(k)) {
      cat("  ", pad(estados[i], w_pi["est"], "left"),
          pad(sprintf("%.4f (%.1f%%)", pi_estacionario[i], pi_estacionario[i] * 100), w_pi["prob"], "center"), "\n", sep = "")
    }
    cat(sep_l(sum(w_pi) + 2), "\n\n")
  }

  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (all(!is.na(pi_estacionario))) {
    est_mais_provavel <- estados[which.max(pi_estacionario)]
    cat(sprintf("  No equil\u00edbrio assint\u00f3tico de longo prazo, o sistema converge para uma distribui\u00e7\u00e3o est\u00e1vel.\n"))
    cat(sprintf("  O estado com maior probabilidade estacion\u00e1ria \u00e9 '%s' com %.1f%% de perman\u00eancia esperada.\n",
                est_mais_provavel, max(pi_estacionario) * 100))
  }
  if (length(nomes_abs) > 0) {
    cat(sprintf("  A presen\u00e7a de estado(s) absorvente(s) [%s] implica que, com probabilidade 1,\n  o processo acabar\u00e1 sendo retido nesses estados.\n",
                paste(nomes_abs, collapse = ", ")))
  }
  cat("\n")
  .print_rodape()

  res <- list(
    matriz_transicao       = P,
    estados                = estados,
    distribuicao_estacionaria = pi_estacionario,
    evolucao_probabilidades = df_evolucao,
    trajetorias            = trajetorias,
    estados_absorventes    = nomes_abs,
    grafico                = p_markov
  )
  class(res) <- c("estatR_markov", "list")
  invisible(res)
}


# ── 2. PASSEIO ALEATÓRIO (RANDOM WALK) ────────────────────────────────────────

#' @title Simulação de Passeio Aleatório (Random Walk)
#' @description Simula trajetórias de passeio aleatório unidimensional simples,
#'   avaliando limites de absorção, tempos de ruína ou retorno à origem.
#'
#' @param n_passos Número total de passos a simular (padrão: 100).
#' @param p_sucesso Probabilidade de incremento unitário positivo (+1). Padrão: 0.5 (simétrico).
#' @param n_trajetorias Número de trajetórias a simular simultaneamente (padrão: 5).
#' @param barreira_superior Limite superior opcional de parada / absorção.
#' @param barreira_inferior Limite inferior opcional de parada / ruína.
#' @param semente Semente para reprodutibilidade.
#' @param grafico Lógico. Se TRUE (padrão), plota os caminhos percorridos.
#'
#' @return Objeto invisível da classe \code{estatR_passeio_aleatorio}.
#' @export
passeio_aleatorio <- function(n_passos = 100,
                              p_sucesso = 0.5,
                              n_trajetorias = 5,
                              barreira_superior = NULL,
                              barreira_inferior = NULL,
                              semente = NULL,
                              grafico = TRUE) {

  if (!is.null(semente)) set.seed(semente)

  df_list <- list()
  posicoes_finais <- numeric(n_trajetorias)

  for (traj in seq_len(n_trajetorias)) {
    passos_vals <- sample(c(1, -1), size = n_passos, replace = TRUE, prob = c(p_sucesso, 1 - p_sucesso))
    pos <- c(0, cumsum(passos_vals))

    if (!is.null(barreira_superior) || !is.null(barreira_inferior)) {
      for (t in seq_along(pos)) {
        if (!is.null(barreira_superior) && pos[t] >= barreira_superior) {
          pos[t:length(pos)] <- barreira_superior
          break
        }
        if (!is.null(barreira_inferior) && pos[t] <= barreira_inferior) {
          pos[t:length(pos)] <- barreira_inferior
          break
        }
      }
    }

    posicoes_finais[traj] <- pos[length(pos)]
    df_list[[traj]] <- data.frame(
      Passo      = 0:n_passos,
      Posicao    = pos,
      Trajetoria = paste0("Caminho ", traj),
      stringsAsFactors = FALSE
    )
  }

  df_caminhos <- do.call(rbind, df_list)

  p_rw <- NULL
  if (isTRUE(grafico)) {
    p_rw <- ggplot2::ggplot(df_caminhos, ggplot2::aes(x = Passo, y = Posicao, color = Trajetoria)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_line(linewidth = 0.9, alpha = 0.85) +
      ggplot2::labs(
        title    = "Simula\u00e7\u00e3o de Passeio Aleat\u00f3rio (Random Walk 1D)",
        subtitle = sprintf("%d Trajet\u00f3rias | %d Passos | p(+1) = %.2f", n_trajetorias, n_passos, p_sucesso),
        x        = "Passo de Tempo",
        y        = "Posi\u00e7\u00e3o (S_n)"
      )

    if (!is.null(barreira_superior)) {
      p_rw <- p_rw + ggplot2::geom_hline(yintercept = barreira_superior, color = "#D90429", linetype = "dotted", linewidth = 1)
    }
    if (!is.null(barreira_inferior)) {
      p_rw <- p_rw + ggplot2::geom_hline(yintercept = barreira_inferior, color = "#D90429", linetype = "dotted", linewidth = 1)
    }

    if (exists("meu_tema", mode = "function")) {
      p_rw <- p_rw + meu_tema()
    } else {
      p_rw <- p_rw + ggplot2::theme_minimal()
    }
    print(p_rw)
  }

  .print_titulo("SIMULA\u00c7\u00c3O: PASSEIO ALEAT\u00d3RIO")
  cat(sprintf("  Total de Passos      : %d\n", n_passos))
  cat(sprintf("  Probabilidade (+1)   : %.2f (%s)\n", p_sucesso, if (p_sucesso == 0.5) "Sim\u00e9trico" else "Com vi\u00e9s/tend\u00eancia"))
  cat(sprintf("  N\u00ba de Trajet\u00f3rias    : %d\n", n_trajetorias))
  cat(sprintf("  Posi\u00e7\u00e3o Final M\u00e9dia : %.2f  (DP: %.2f)\n\n", mean(posicoes_finais), stats::sd(posicoes_finais)))

  .print_topico("INTERPRETA\u00c7\u00c3O")
  if (p_sucesso == 0.5) {
    cat("  O passeio \u00e9 n\u00e3o-viesado (esperan\u00e7a zero). A dispers\u00e3o t\u00edpica cresce com a raiz quadrada do tempo sqrt(n).\n")
  } else if (p_sucesso > 0.5) {
    cat("  O passeio apresenta deriva positiva (drift para cima com probabilidade tendenciosa).\n")
  } else {
    cat("  O passeio apresenta deriva negativa (drift descendente com risco de ru\u00edna r\u00e1pida).\n")
  }
  cat("\n")
  .print_rodape()

  res <- list(
    caminhos        = df_caminhos,
    posicoes_finais = posicoes_finais,
    grafico         = p_rw
  )
  class(res) <- c("estatR_passeio_aleatorio", "list")
  invisible(res)
}
