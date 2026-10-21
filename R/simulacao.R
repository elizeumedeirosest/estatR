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
#' @title Bootstrap: Estimação por Reamostragem (Simples, Estratificado, Residual e por Bloco)
#' @description Estima a variabilidade, viés, erro padrão e intervalos de confiança não-paramétricos
#'   usando reamostragem Bootstrap. Suporta quatro modalidades principais:
#'   \itemize{
#'     \item \code{"simples"}: Reamostragem padrão independente com reposição sobre vetor.
#'     \item \code{"estratificado"}: Reamostragem balanceada dentro de subgrupos/estratos.
#'     \item \code{"residual"}: Reamostragem de resíduos para modelos de regressão (\code{lm}),
#'       gerando distribuições e ICs robustos para os coeficientes \eqn{\beta} e \eqn{R^2}.
#'     \item \code{"bloco"}: Moving Block Bootstrap para dados com dependência temporal.
#'   }
#'
#' @param x Vetor numérico ou série temporal (para \code{"simples"}, \code{"estratificado"} ou \code{"bloco"}).
#' @param estatistica A estatística de interesse (padrão \code{"media"}). Pode ser
#'   \code{"media"}, \code{"mediana"}, \code{"variancia"}, \code{"dp"}, \code{"proporcao"}
#'   ou uma função customizada.
#' @param tipo Modalidade de Bootstrap: \code{"simples"} (padrão), \code{"estratificado"},
#'   \code{"residual"} (para regressão) ou \code{"bloco"} (para séries temporais).
#' @param estrato Vetor de fatores ou caracteres definindo os estratos (obrigatório se \code{tipo = "estratificado"}).
#' @param modelo Objeto da classe \code{lm} (obrigatório se \code{tipo = "residual"}).
#' @param bloco Tamanho do bloco para reamostragem temporal (se \code{NULL}, calculado automaticamente).
#' @param repeticoes Número de amostras Bootstrap a gerar (padrão 1000).
#' @param confianca Nível de confiança para o IC Bootstrap (padrão 0.95).
#' @param semente Semente para reprodutibilidade (padrão \code{NULL}).
#' @param grafico Lógico. Se TRUE (padrão), plota gráfico com \code{grade = "dupla"}.
#'
#' @return Retorna invisivelmente uma lista com as estimativas Bootstrap e métricas.
#' @export
# Complete bootstrap function definition for test
bootstrap <- function(x = NULL,
                          estatistica = "media",
                          tipo = c("simples", "estratificado", "residual", "bloco"),
                          estrato = NULL,
                          modelo = NULL,
                          bloco = NULL,
                          repeticoes = 1000,
                          confianca = 0.95,
                          semente = NULL,
                          grafico = TRUE) {
  if (!is.null(semente)) set.seed(semente)
  tipo <- match.arg(tipo)

  # Suporte polimórfico: se x for lm e tipo não foi forçado, vira residual
  if (inherits(x, "lm") && is.null(modelo)) {
    modelo <- x
    tipo <- "residual"
    x <- NULL
  }

  alpha <- 1 - confianca

  # ── 1. BOOTSTRAP RESIDUAL PARA REGRESSÃO ──────────────────────────────────
  if (tipo == "residual") {
    if (is.null(modelo) || !inherits(modelo, "lm")) {
      stop("Para 'tipo = residual', deve ser fornecido um modelo ajustado da classe 'lm'.")
    }

    mf <- stats::model.frame(modelo)
    X  <- stats::model.matrix(modelo)
    y  <- stats::model.response(mf)
    y_hat <- stats::fitted(modelo)
    res_orig <- stats::residuals(modelo)
    n <- length(res_orig)
    p <- ncol(X)
    termos_nomes <- colnames(X)

    # Resíduos centrados (premissa de média zero)
    res_centrados <- res_orig - mean(res_orig)
    coef_orig <- stats::coef(modelo)
    ss_tot <- sum((y - mean(y))^2)
    r2_orig <- if (ss_tot > 0) 1 - sum(res_orig^2) / ss_tot else 0

    # Matriz para armazenar coeficientes bootstrap
    mat_boot_coefs <- matrix(NA, nrow = repeticoes, ncol = p)
    vec_boot_r2    <- numeric(repeticoes)

    XtX_inv_Xt <- solve(t(X) %*% X) %*% t(X)

    for (b in seq_len(repeticoes)) {
      e_star <- sample(res_centrados, size = n, replace = TRUE)
      y_star <- y_hat + e_star
      beta_star <- as.vector(XtX_inv_Xt %*% y_star)
      mat_boot_coefs[b, ] <- beta_star

      res_star <- y_star - (X %*% beta_star)
      ss_tot_star <- sum((y_star - mean(y_star))^2)
      vec_boot_r2[b] <- if (ss_tot_star > 0) max(0, 1 - sum(res_star^2) / ss_tot_star) else 0
    }

    # Métricas para cada coeficiente + R²
    df_res_list <- list()
    for (j in seq_len(p)) {
      b_vals <- mat_boot_coefs[, j]
      b_orig <- coef_orig[j]
      b_media <- mean(b_vals)
      b_vies <- b_media - b_orig
      b_ep <- stats::sd(b_vals)
      b_ic <- stats::quantile(b_vals, probs = c(alpha / 2, 1 - alpha / 2))

      df_res_list[[j]] <- data.frame(
        Parâmetro = termos_nomes[j],
        Original  = b_orig,
        Bootstrap = b_media,
        Viés      = b_vies,
        Erro_Padrao = b_ep,
        IC_Inf    = b_ic[1],
        IC_Sup    = b_ic[2],
        stringsAsFactors = FALSE
      )
    }

    # Adiciona R²
    r2_media <- mean(vec_boot_r2)
    r2_ic <- stats::quantile(vec_boot_r2, probs = c(alpha / 2, 1 - alpha / 2))
    df_res_list[[p + 1]] <- data.frame(
      Parâmetro = "R²",
      Original  = r2_orig,
      Bootstrap = r2_media,
      Viés      = r2_media - r2_orig,
      Erro_Padrao = stats::sd(vec_boot_r2),
      IC_Inf    = r2_ic[1],
      IC_Sup    = r2_ic[2],
      stringsAsFactors = FALSE
    )

    tabela_final <- do.call(rbind, df_res_list)

    .print_titulo(sprintf("BOOTSTRAP RESIDUAL (%d REPETIÇÕES)", repeticoes))
    cat(sprintf("  Modelo Linear:   %s (%d observações)\n", deparse(stats::formula(modelo)), n))
    cat(sprintf("  Nível Confiança: %.0f%% (Intervalos por Percentis)\n\n", confianca * 100))

    .print_topico("ESTIMATIVAS E INTERVALOS BOOTSTRAP")
    df_print <- data.frame(
      "Parâmetro"    = tabela_final$Parâmetro,
      "Original"     = sapply(tabela_final$Original, .fmt_num, decimais = 3),
      "Média Boot"   = sapply(tabela_final$Bootstrap, .fmt_num, decimais = 3),
      "Viés"         = sapply(tabela_final$Viés, .fmt_num, decimais = 4),
      "Erro Padrão"  = sapply(tabela_final$Erro_Padrao, .fmt_num, decimais = 4),
      "IC Inf"       = sapply(tabela_final$IC_Inf, .fmt_num, decimais = 3),
      "IC Sup"       = sapply(tabela_final$IC_Sup, .fmt_num, decimais = 3),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    .print_tabela_estatR(df_print, align = c("left", "center", "center", "center", "center", "center", "center"))
    .print_rodape()

    resultado <- list(
      tipo            = "residual",
      repeticoes      = repeticoes,
      confianca       = confianca,
      tabela          = tabela_final,
      coefs_bootstrap = mat_boot_coefs,
      r2_bootstrap    = vec_boot_r2
    )
    return(invisible(resultado))
  }

  # ── 2. PREPARAÇÃO PARA DEMAIS MODALIDADES (Simples, Estratificado, Bloco) ──
  if (is.null(x)) stop("O argumento 'x' deve ser fornecido para esta modalidade de bootstrap.")
  x_c <- na.omit(x)
  n <- length(x_c)
  if (n < 5) stop("São necessárias pelo menos 5 observações válidas.")

  x_expr   <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x") || length(var_nome) == 0) var_nome <- "Variável"

  # Define função estatística
  if (is.function(estatistica)) {
    fn_stat   <- estatistica
    stat_nome <- "Customizada"
  } else {
    stat_nome <- switch(estatistica,
      media     = "Média",
      mediana   = "Mediana",
      variancia = "Variância",
      dp        = "Desvio Padrão",
      proporcao = "Proporção (x > 0)",
      stop("'estatistica' inválida. Use: media, mediana, variancia, dp, proporcao ou uma função.")
    )
    fn_stat <- switch(estatistica,
      media     = mean,
      mediana   = stats::median,
      variancia = stats::var,
      dp        = stats::sd,
      proporcao = function(v) mean(v > 0)
    )
  }

  est_original <- fn_stat(x_c)

  # ── 3. EXECUÇÃO DA REAMOSTRAGEM ───────────────────────────────────────────
  boot_vals <- numeric(repeticoes)

  if (tipo == "simples") {
    boot_vals <- replicate(repeticoes, fn_stat(sample(x_c, n, replace = TRUE)))

  } else if (tipo == "estratificado") {
    if (is.null(estrato)) {
      stop("Para 'tipo = estratificado', é necessário informar o argumento 'estrato'.")
    }
    if (length(estrato) != length(x)) {
      stop("O argumento 'estrato' deve ter o mesmo comprimento que 'x'.")
    }
    # Filtra NAs conjuntos
    idx_val <- which(!is.na(x) & !is.na(estrato))
    x_val   <- x[idx_val]
    est_val <- factor(estrato[idx_val])
    n_val   <- length(x_val)

    split_idx <- split(seq_len(n_val), est_val)
    n_estratos <- length(split_idx)

    for (b in seq_len(repeticoes)) {
      idx_boot <- unlist(lapply(split_idx, function(ids) {
        sample(ids, size = length(ids), replace = TRUE)
      }), use.names = FALSE)
      boot_vals[b] <- fn_stat(x_val[idx_boot])
    }

  } else if (tipo == "bloco") {
    # Moving Block Bootstrap
    b_len <- if (!is.null(bloco)) as.integer(bloco) else {
      if (stats::is.ts(x)) max(2, stats::frequency(x)) else max(2, floor(n^(1/3)))
    }
    n_blocos_totais <- n - b_len + 1
    if (n_blocos_totais < 2) stop("Tamanho do bloco muito grande para o comprimento da série.")

    num_blocos_necessarios <- ceiling(n / b_len)

    for (b in seq_len(repeticoes)) {
      start_pts <- sample(seq_len(n_blocos_totais), size = num_blocos_necessarios, replace = TRUE)
      series_b  <- unlist(lapply(start_pts, function(st) x_c[st:(st + b_len - 1)]), use.names = FALSE)
      boot_vals[b] <- fn_stat(series_b[1:n])
    }
  }

  est_boot <- mean(boot_vals)
  vies     <- est_boot - est_original
  ep_boot  <- stats::sd(boot_vals)
  linf     <- stats::quantile(boot_vals, alpha / 2)
  lsup     <- stats::quantile(boot_vals, 1 - alpha / 2)

  # ── SAÍDA NO CONSOLE ──────────────────────────────────────────────────────
  rotulo_tipo <- switch(
    tipo,
    simples       = "BOOTSTRAP SIMPLES",
    estratificado = sprintf("BOOTSTRAP ESTRATIFICADO (%d ESTRATOS)", length(unique(estrato))),
    bloco         = sprintf("BOOTSTRAP DE BLOCOS MÓVEIS (BLOCO = %d)", b_len)
  )

  .print_titulo(sprintf("%s — %d REPETIÇÕES", rotulo_tipo, repeticoes))
  cat(sprintf("  Variável:     %s (N = %d)\n", var_nome, n))
  cat(sprintf("  Estatística:  %s\n", stat_nome))
  cat(sprintf("  Confiança IC: %.0f%% (Percentis)\n\n", confianca * 100))

  .print_topico("RESULTADOS DA REAMOSTRAGEM")
  df_saida <- data.frame(
    "Métrica" = c("Estimativa Original", "Média Bootstrap", "Viés Estimado", "Erro Padrão Bootstrap", "IC Limite Inferior", "IC Limite Superior"),
    "Valor"   = c(.fmt_num(est_original, 4), .fmt_num(est_boot, 4), .fmt_num(vies, 4), .fmt_num(ep_boot, 4), .fmt_num(linf, 4), .fmt_num(lsup, 4)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  .print_tabela_estatR(df_saida, align = c("left", "center"))

  vies_txt <- if (abs(vies) < 0.001 * pmax(1e-10, abs(est_original))) "praticamente nulo" else
              if (vies > 0) "positivo (tendência à superestimação)" else "negativo (tendência à subestimação)"

  cat(sprintf("  * O viés da reamostragem é %s (%s).\n", vies_txt, .fmt_num(vies, 4)))
  cat(sprintf("  * Com %.0f%% de confiança, a estatística verdadeira estimada situa-se em [%s ; %s].\n\n",
              confianca * 100, .fmt_num(linf, 4), .fmt_num(lsup, 4)))

  .print_rodape()

  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch({
      tema_padrao <- if (exists("meu_tema", mode = "function")) {
        meu_tema(grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }

      df_boot_g <- data.frame(vals = boot_vals)
      p_plot <- ggplot2::ggplot(df_boot_g, ggplot2::aes(x = vals)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                bins = min(50, ceiling(sqrt(repeticoes))),
                                fill = "#888888", color = "white", alpha = 0.8) +
        ggplot2::geom_density(color = "#2B5C8F", linewidth = 1.0) +
        ggplot2::geom_vline(xintercept = est_original, color = "#D90429",
                            linetype = "dashed", linewidth = 1.0) +
        ggplot2::geom_vline(xintercept = c(linf, lsup), color = "#2B5C8F",
                            linetype = "dotted", linewidth = 0.9) +
        ggplot2::labs(
          title    = sprintf("Distribuição Bootstrap — %s (%s)", stat_nome, var_nome),
          subtitle = sprintf("%d repetições | IC %.0f%%: [%s ; %s]",
                             repeticoes, confianca * 100, .fmt_num(linf, 3), .fmt_num(lsup, 3)),
          x = stat_nome, y = "Densidade"
        ) +
        tema_padrao

      print(p_plot)
    }, error = function(e) {
      message("[Aviso] Não foi possível plotar o gráfico bootstrap: ", e$message)
    })
  }

  resultado <- list(
    tipo                = tipo,
    estatistica         = stat_nome,
    estimativa_original = est_original,
    media_bootstrap     = est_boot,
    vies                = vies,
    ep_bootstrap        = ep_boot,
    ic_inf              = linf,
    ic_sup              = lsup,
    valores_bootstrap   = boot_vals
  )
  invisible(resultado)
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
  .print_titulo("TEOREMA CENTRAL DO LIMITE — SIMULA\u00c7\u00c3O %s")
  cat(sprintf("  Popula\u00e7\u00e3o: %s   |   N pop: %d   |   Repeti\u00e7\u00f5es: %d\n", var_nome, n, repeticoes))
  cat(sprintf("  Tamanho de cada amostra (n): %d\n\n", n_amostra))

  .print_topico("COMPARAÇÃO: POPULAÇÃO vs DISTRIBUIÇÃO DAS MÉDIAS")
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

  .print_topico("NORMALIDADE DAS MÉDIAS SIMULADAS (Shapiro-Wilk)")
  formata_p <- function(p) if (p < 0.001) "< 0.001" else formatC(p, format="f", digits=3, decimal.mark=",")
  cor_ok   <- "\033[1;32m"
  cor_warn <- "\033[1;33m"
  reset    <- "\033[0m"

  cat(sprintf("  W = %s   |   p = %s\n", .sim_fmt(st$statistic, 4), formata_p(st$p.value)))
  if (st$p.value >= 0.05) {
    cat(sprintf("  %s[\u2713]%s As m\u00e9dias simuladas seguem distribui\u00e7\u00e3o Normal (n = %d confirma o TCL).\n\n", cor_ok, reset, n_amostra))
  } else {
    cat(sprintf("  %s[!]%s As m\u00e9dias ainda n\u00e3o s\u00e3o perfeitamente normais. Tente aumentar 'n_amostra'.\n\n", cor_warn, reset))
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
  .print_titulo("SIMULA\u00c7\u00c3O DE MONTE CARLO %s")
  cat(sprintf("  Repeti\u00e7\u00f5es: %d   |   Confian\u00e7a: %.0f%%\n\n", repeticoes, confianca * 100))

  .print_topico("RESULTADOS")
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

  .print_topico(sprintf("INTERVALO DE CONFIANÇA (%.0f%% — Wald)", confianca * 100))
  cat(.sim_sep(nchar(hdr)), "\n")
  ic_str <- sprintf("[ %s%%  ;  %s%% ]", .sim_fmt(linf * 100, 2), .sim_fmt(lsup * 100, 2))
  cat(sprintf("  %s\n", .sim_pad(ic_str, nchar(hdr), "center")))
  cat(.sim_sep(nchar(hdr)), "\n\n")

  .print_topico("INTERPRETAÇÃO")
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
