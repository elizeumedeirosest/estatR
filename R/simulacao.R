# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: SIMULAÇÃO
# Bootstrap, Teorema Central do Limite e Monte Carlo
# ─────────────────────────────────────────────────────────────────────────────

# ── Helpers internos ──────────────────────────────────────────────────────────
.sim_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
  paste0(s, strrep(" ", sp))
}

.sim_sep <- function(w) paste0("  ", strrep("\u2500", w))
.sim_fmt <- function(x, d = 4) formatC(x, format = "f", digits = d, decimal.mark = ",")


# ─────────────────────────────────────────────────────────────────────────────
#' @title Bootstrap: Estimação por Reamostragem
#' @description Estima a variabilidade de uma estatística usando reamostragem Bootstrap.
#'   Calcula a estimativa original, o viés estimado, o erro padrão Bootstrap e o
#'   Intervalo de Confiança por percentis.
#' @param x Vetor numérico com os dados observados.
#' @param estatistica A estatística de interesse. Pode ser um nome entre
#'   \code{"media"}, \code{"mediana"}, \code{"variancia"}, \code{"dp"} ou
#'   \code{"proporcao"} (proporção de valores > 0), ou uma \strong{função} customizada
#'   (ex: \code{estatistica = function(x) quantile(x, 0.9)}).
#' @param repeticoes Número de amostras Bootstrap a gerar (padrão 1000).
#' @param confianca Nível de confiança para o IC Bootstrap (padrão 0.95).
#' @param semente Semente para reprodutibilidade (padrão \code{NULL}).
#' @param grafico Lógico. Se TRUE, plota o histograma das estimativas Bootstrap.
#' @return Retorna invisivelmente um data frame com as estimativas Bootstrap.
#' @export
bootstrap <- function(x, estatistica = "media", repeticoes = 1000,
                      confianca = 0.95, semente = NULL, grafico = TRUE) {
  if (!is.null(semente)) set.seed(semente)

  x_expr   <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x")) var_nome <- "x"

  x_c <- x[!is.na(x)]
  n   <- length(x_c)
  if (n < 5) stop("S\u00e3o necess\u00e1rias pelo menos 5 observa\u00e7\u00f5es.")

  # Define a função da estatística
  if (is.function(estatistica)) {
    fn_stat   <- estatistica
    stat_nome <- "Customizada"
  } else {
    stat_nome <- switch(estatistica,
      media     = "M\u00e9dia",
      mediana   = "Mediana",
      variancia = "Vari\u00e2ncia",
      dp        = "Desvio Padr\u00e3o",
      proporcao = "Propor\u00e7\u00e3o (x > 0)",
      stop("'estatistica' inv\u00e1lida. Use: media, mediana, variancia, dp, proporcao ou uma fun\u00e7\u00e3o.")
    )
    fn_stat <- switch(estatistica,
      media     = mean,
      mediana   = median,
      variancia = var,
      dp        = sd,
      proporcao = function(v) mean(v > 0)
    )
  }

  # Estimativa original
  est_original <- fn_stat(x_c)

  # Reamostras Bootstrap
  boot_vals <- replicate(repeticoes, fn_stat(sample(x_c, n, replace = TRUE)))

  # Resultados
  est_boot <- mean(boot_vals)
  vies     <- est_boot - est_original
  ep_boot  <- sd(boot_vals)
  alpha    <- 1 - confianca
  linf     <- quantile(boot_vals, alpha / 2)
  lsup     <- quantile(boot_vals, 1 - alpha / 2)

  w_sep <- 67
  cat(sprintf("\n\u2500\u2500 BOOTSTRAP (%d REPETI\u00c7\u00f5ES) %s\n", repeticoes, strrep("\u2500", w_sep - 21 - nchar(as.character(repeticoes)))))
  cat(sprintf("  Vari\u00e1vel: %s   |   N: %d   |   Estat\u00edstica: %s\n\n", var_nome, n, stat_nome))

  cat("  \u25b6 RESULTADOS\n")
  wc <- c(met = 30, val = 20)
  hdr <- paste0(.sim_pad("M\u00e9trica", wc["met"], "left"), .sim_pad("Valor", wc["val"], "center"))
  cat(.sim_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.sim_sep(nchar(hdr)), "\n")
  cat("  ", .sim_pad("Estimativa Original",       wc["met"], "left"), .sim_pad(.sim_fmt(est_original, 4), wc["val"], "center"), "\n", sep = "")
  cat("  ", .sim_pad("M\u00e9dia Bootstrap",       wc["met"], "left"), .sim_pad(.sim_fmt(est_boot,     4), wc["val"], "center"), "\n", sep = "")
  cat("  ", .sim_pad("Vi\u00e9s Estimado",          wc["met"], "left"), .sim_pad(.sim_fmt(vies,         4), wc["val"], "center"), "\n", sep = "")
  cat("  ", .sim_pad("Erro Padr\u00e3o Bootstrap", wc["met"], "left"), .sim_pad(.sim_fmt(ep_boot,      4), wc["val"], "center"), "\n", sep = "")
  cat(.sim_sep(nchar(hdr)), "\n\n")

  cat(sprintf("  \u25b6 INTERVALO DE CONFIAN\u00c7A BOOTSTRAP (%.0f%% — Percentis)\n", confianca * 100))
  cat(.sim_sep(nchar(hdr)), "\n")
  ic_str <- sprintf("[ %s  ;  %s ]", .sim_fmt(linf, 4), .sim_fmt(lsup, 4))
  cat(sprintf("  %s\n", .sim_pad(ic_str, nchar(hdr), "center")))
  cat(.sim_sep(nchar(hdr)), "\n\n")

  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  cat(sprintf("  A partir de %d reamostras Bootstrap, o erro padr\u00e3o\n", repeticoes))
  cat(sprintf("  estimado da %s foi %s.\n\n", tolower(stat_nome), .sim_fmt(ep_boot, 4)))
  vies_txt <- if (abs(vies) < 0.001 * abs(est_original)) "praticamente nulo" else
              if (vies > 0) "positivo (superestima\u00e7\u00e3o)" else "negativo (subestima\u00e7\u00e3o)"
  cat(sprintf("  O vi\u00e9s Bootstrap \u00e9 %s (%s).\n\n", vies_txt, .sim_fmt(vies, 4)))
  cat(sprintf("  Com %.0f%% de confian\u00e7a, a verdadeira %s\n  est\u00e1 entre %s e %s.\n",
              confianca * 100, tolower(stat_nome), .sim_fmt(linf, 4), .sim_fmt(lsup, 4)))
  .print_rodape()

  # Gráfico
  if (grafico && exists("meu_tema")) {
    tryCatch({
      df_boot <- data.frame(vals = boot_vals)
      old_w <- getOption("warn"); options(warn = -1)

      p_plot <- ggplot2::ggplot(df_boot, ggplot2::aes(x = vals)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                bins = min(50, ceiling(sqrt(repeticoes))),
                                fill = "#888888", color = "white", alpha = 0.8) +
        ggplot2::geom_density(color = "#333333", linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = est_original, color = "#D90429",
                            linetype = "dashed", linewidth = 1) +
        ggplot2::geom_vline(xintercept = c(linf, lsup), color = "#1B4F72",
                            linetype = "dotted", linewidth = 0.9) +
        ggplot2::annotate("text", x = est_original, y = Inf,
                          label = sprintf(" Original\n %s", .sim_fmt(est_original, 3)),
                          hjust = -0.05, vjust = 1.5, size = 4,
                          fontface = "bold", color = "#D90429") +
        ggplot2::labs(
          title    = sprintf("Distribui\u00e7\u00e3o Bootstrap — %s de %s", stat_nome, var_nome),
          subtitle = sprintf("%d repeti\u00e7\u00f5es   |   IC %.0f%%: [%s ; %s]",
                             repeticoes, confianca * 100, .sim_fmt(linf, 3), .sim_fmt(lsup, 3)),
          x = stat_nome, y = "Densidade", caption = "estatR"
        ) +
        meu_tema()

      suppressMessages(print(p_plot))
      options(warn = old_w)
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }

  invisible(data.frame(
    estimativa_original = est_original, media_bootstrap = est_boot,
    vies = vies, ep_bootstrap = ep_boot, li = linf, ls = lsup
  ))
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Teorema Central do Limite (TCL) — Simulação Visual
#' @description Demonstra visualmente o Teorema Central do Limite: extrai amostras
#'   da população fornecida e mostra que a distribuição das médias tende à Normal.
#' @param x Vetor numérico representando a população (ou uma grande amostra dela).
#' @param n_amostra Tamanho de cada amostra a extrair (padrão 30). Teste com 5, 10, 50 para ver o efeito.
#' @param repeticoes Número de amostras a extrair (padrão 1000).
#' @param semente Semente para reprodutibilidade (padrão \code{NULL}).
#' @return Retorna invisivelmente o vetor das médias simuladas.
#' @export
simular_tcl <- function(x, n_amostra = 30, repeticoes = 1000, semente = NULL) {
  if (!is.null(semente)) set.seed(semente)

  x_expr   <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x")) var_nome <- "x"

  x_c <- x[!is.na(x)]
  n   <- length(x_c)
  if (n < n_amostra) stop("A popula\u00e7\u00e3o deve ter pelo menos 'n_amostra' observa\u00e7\u00f5es.")

  # Simular
  medias <- replicate(repeticoes, mean(sample(x_c, n_amostra, replace = TRUE)))

  mu_pop  <- mean(x_c)
  dp_pop  <- sd(x_c)
  ep_teo  <- dp_pop / sqrt(n_amostra)
  mu_boot <- mean(medias)
  dp_boot <- sd(medias)

  # Normalidade das médias simuladas
  st <- shapiro.test(sample(medias, min(length(medias), 5000)))

  w_sep <- 67
  cat(sprintf("\n\u2500\u2500 TEOREMA CENTRAL DO LIMITE — SIMULA\u00c7\u00c3O %s\n", strrep("\u2500", w_sep - 42)))
  cat(sprintf("  Popula\u00e7\u00e3o: %s   |   N pop: %d   |   Repeti\u00e7\u00f5es: %d\n", var_nome, n, repeticoes))
  cat(sprintf("  Tamanho de cada amostra (n): %d\n\n", n_amostra))

  cat("  \u25b6 COMPARA\u00c7\u00c3O: POPULA\u00c7\u00c3O vs DISTRIBUI\u00c7\u00c3O DAS M\u00c9DIAS\n")
  wc <- c(met = 30, pop = 16, med = 16)
  hdr <- paste0(.sim_pad("M\u00e9trica",              wc["met"], "left"),
                .sim_pad("Popula\u00e7\u00e3o",        wc["pop"], "center"),
                .sim_pad("M\u00e9dias Simuladas", wc["med"], "center"))
  cat(.sim_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.sim_sep(nchar(hdr)), "\n")
  cat("  ", .sim_pad("M\u00e9dia",           wc["met"], "left"), .sim_pad(.sim_fmt(mu_pop, 3),  wc["pop"], "center"), .sim_pad(.sim_fmt(mu_boot, 3), wc["med"], "center"), "\n", sep = "")
  cat("  ", .sim_pad("Desvio Padr\u00e3o",   wc["met"], "left"), .sim_pad(.sim_fmt(dp_pop, 3),  wc["pop"], "center"), .sim_pad(.sim_fmt(dp_boot, 3), wc["med"], "center"), "\n", sep = "")
  cat("  ", .sim_pad("EP Te\u00f3rico (\u03c3/\u221an)", wc["met"], "left"), .sim_pad("—",  wc["pop"], "center"), .sim_pad(.sim_fmt(ep_teo, 3), wc["med"], "center"), "\n", sep = "")
  cat(.sim_sep(nchar(hdr)), "\n\n")

  cat("  \u25b6 NORMALIDADE DAS M\u00c9DIAS SIMULADAS (Shapiro-Wilk)\n")
  formata_p <- function(p) if (p < 0.001) "< 0.001" else formatC(p, format="f", digits=3, decimal.mark=",")
  cat(sprintf("  W = %s   |   p = %s\n", .sim_fmt(st$statistic, 4), formata_p(st$p.value)))
  if (st$p.value >= 0.05) {
    cat(sprintf("  [\u2713] As m\u00e9dias simuladas seguem distribui\u00e7\u00e3o Normal (n = %d confirma o TCL).\n\n", n_amostra))
  } else {
    cat(sprintf("  [!] As m\u00e9dias ainda n\u00e3o s\u00e3o perfeitamente normais. Tente aumentar 'n_amostra'.\n\n"))
  }
  .print_rodape()

  # Gráfico: painel 2 lados
  if (exists("meu_tema")) {
    tryCatch({
      old_w <- getOption("warn"); options(warn = -1)

      df_pop  <- data.frame(x = x_c)
      df_med  <- data.frame(x = medias)

      p1 <- ggplot2::ggplot(df_pop, ggplot2::aes(x = x)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                bins = min(30, ceiling(sqrt(n))),
                                fill = "#888888", color = "white", alpha = 0.8) +
        ggplot2::geom_density(color = "#D90429", linewidth = 0.9) +
        ggplot2::labs(title = sprintf("Popula\u00e7\u00e3o Original\n(%s)", var_nome),
                      x = var_nome, y = "Densidade") +
        meu_tema()

      x_med_seq <- seq(min(medias), max(medias), length.out = 300)
      df_norm   <- data.frame(x = x_med_seq, y = dnorm(x_med_seq, mean = mu_boot, sd = ep_teo))

      p2 <- ggplot2::ggplot(df_med, ggplot2::aes(x = x)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                bins = min(50, ceiling(sqrt(repeticoes))),
                                fill = "#555555", color = "white", alpha = 0.8) +
        ggplot2::geom_line(data = df_norm, ggplot2::aes(x = x, y = y),
                           color = "#D90429", linewidth = 1.1) +
        ggplot2::labs(title = sprintf("Distribui\u00e7\u00e3o das M\u00e9dias\n(%d amostras, n = %d)", repeticoes, n_amostra),
                      x = sprintf("M\u00e9dia de %s", var_nome), y = "Densidade",
                      caption = "estatR") +
        meu_tema()

      if (requireNamespace("patchwork", quietly = TRUE)) {
        print(p1 | p2)
      } else {
        suppressMessages(print(p1))
        suppressMessages(print(p2))
      }

      options(warn = old_w)
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }

  invisible(medias)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Simulação de Monte Carlo
#' @description Estima a probabilidade de um evento por simulação repetida.
#'   O usuário define um experimento aleatório como uma função que retorna TRUE/FALSE.
#' @param experimento Função sem argumentos que representa um experimento e retorna um
#'   valor lógico (\code{TRUE}/\code{FALSE}) ou numérico (1/0). Exemplo:
#'   \code{function() mean(rnorm(30)) > 1}.
#' @param repeticoes Número de repetições da simulação (padrão 10000).
#' @param confianca Nível de confiança para o IC da proporção (padrão 0.95).
#' @param semente Semente para reprodutibilidade (padrão \code{NULL}).
#' @param grafico Lógico. Se TRUE, plota a convergência da probabilidade estimada.
#' @return Retorna invisivelmente a probabilidade estimada.
#' @export
monte_carlo <- function(experimento, repeticoes = 10000, confianca = 0.95,
                        semente = NULL, grafico = TRUE) {
  if (!is.function(experimento)) stop("'experimento' deve ser uma fun\u00e7\u00e3o.")
  if (!is.null(semente)) set.seed(semente)

  resultados <- as.numeric(replicate(repeticoes, experimento()))

  p_hat  <- mean(resultados)
  ep     <- sqrt(p_hat * (1 - p_hat) / repeticoes)
  zc     <- qnorm(1 - (1 - confianca) / 2)
  linf   <- max(0, p_hat - zc * ep)
  lsup   <- min(1, p_hat + zc * ep)
  margem <- zc * ep

  # Convergência (médias acumuladas)
  p_acum <- cumsum(resultados) / seq_along(resultados)

  w_sep <- 67
  cat(sprintf("\n\u2500\u2500 SIMULA\u00c7\u00c3O DE MONTE CARLO %s\n", strrep("\u2500", w_sep - 26)))
  cat(sprintf("  Repeti\u00e7\u00f5es: %d   |   Confian\u00e7a: %.0f%%\n\n", repeticoes, confianca * 100))

  cat("  \u25b6 RESULTADOS\n")
  wc <- c(met = 30, val = 20)
  hdr <- paste0(.sim_pad("M\u00e9trica", wc["met"], "left"), .sim_pad("Valor", wc["val"], "center"))
  cat(.sim_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.sim_sep(nchar(hdr)), "\n")
  cat("  ", .sim_pad("Probabilidade Estimada (p\u0302)", wc["met"], "left"),
      .sim_pad(sprintf("%s%%", .sim_fmt(p_hat * 100, 2)), wc["val"], "center"), "\n", sep = "")
  cat("  ", .sim_pad("Erro Padr\u00e3o",        wc["met"], "left"), .sim_pad(.sim_fmt(ep, 5),     wc["val"], "center"), "\n", sep = "")
  cat("  ", .sim_pad("Margem de Erro",      wc["met"], "left"), .sim_pad(.sim_fmt(margem, 5), wc["val"], "center"), "\n", sep = "")
  cat("  ", .sim_pad("Sucessos Observados", wc["met"], "left"), .sim_pad(sum(resultados),   wc["val"], "center"), "\n", sep = "")
  cat(.sim_sep(nchar(hdr)), "\n\n")

  cat(sprintf("  \u25b6 INTERVALO DE CONFIAN\u00c7A (%.0f%% — Wald)\n", confianca * 100))
  cat(.sim_sep(nchar(hdr)), "\n")
  ic_str <- sprintf("[ %s%%  ;  %s%% ]", .sim_fmt(linf * 100, 2), .sim_fmt(lsup * 100, 2))
  cat(sprintf("  %s\n", .sim_pad(ic_str, nchar(hdr), "center")))
  cat(.sim_sep(nchar(hdr)), "\n\n")

  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  cat(sprintf("  Com base em %d simula\u00e7\u00f5es, a probabilidade estimada\n", repeticoes))
  cat(sprintf("  do evento \u00e9 de %s%% (IC %.0f%%: %s%% a %s%%).\n",
              .sim_fmt(p_hat * 100, 2), confianca * 100,
              .sim_fmt(linf * 100, 2), .sim_fmt(lsup * 100, 2)))
  cat(sprintf("  A margem de erro da estimativa \u00e9 de %s%%.\n", .sim_fmt(margem * 100, 3)))
  .print_rodape()

  # Gráfico de convergência
  if (grafico && exists("meu_tema")) {
    tryCatch({
      old_w <- getOption("warn"); options(warn = -1)

      df_conv <- data.frame(n = seq_along(p_acum), p = p_acum)

      p_plot <- ggplot2::ggplot(df_conv, ggplot2::aes(x = n, y = p * 100)) +
        ggplot2::geom_line(color = "#555555", linewidth = 0.7, alpha = 0.8) +
        ggplot2::geom_hline(yintercept = p_hat * 100, color = "#D90429",
                            linetype = "dashed", linewidth = 1) +
        ggplot2::geom_hline(yintercept = c(linf * 100, lsup * 100), color = "#1B4F72",
                            linetype = "dotted", linewidth = 0.8) +
        ggplot2::annotate("text", x = repeticoes, y = p_hat * 100,
                          label = sprintf("  p\u0302 = %s%%", .sim_fmt(p_hat * 100, 2)),
                          hjust = 1, vjust = -0.5, size = 4.5,
                          fontface = "bold", color = "#D90429") +
        ggplot2::scale_x_continuous(labels = scales::comma_format(big.mark = ".", decimal.mark = ",")) +
        ggplot2::labs(
          title    = "Converg\u00eancia da Simula\u00e7\u00e3o de Monte Carlo",
          subtitle = sprintf("Probabilidade estabilizando em %s%%  |  IC %.0f%%: [%s%% ; %s%%]",
                             .sim_fmt(p_hat * 100, 2), confianca * 100,
                             .sim_fmt(linf * 100, 2), .sim_fmt(lsup * 100, 2)),
          x = "N\u00famero de Simula\u00e7\u00f5es",
          y = "Probabilidade Acumulada (%)",
          caption = "estatR"
        ) +
        meu_tema()

      suppressMessages(print(p_plot))
      options(warn = old_w)
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }

  invisible(p_hat)
}
