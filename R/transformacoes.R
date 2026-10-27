# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: TRANSFORMAÇÕES DE VARIÁVEIS (BOX-COX, YEO-JOHNSON E FAMÍLIA DE POTÊNCIAS)
# ─────────────────────────────────────────────────────────────────────────────

#' @title Transformação de Variáveis e Seleção Ótima de Potência
#' @description Aplica transformações estatísticas clássicas e modernas sobre vetores
#'   ou conjuntos de dados numéricos para aproximação da normalidade, estabilização
#'   de variância ou reescalonamento. Suporta estimação ótima de Box-Cox e Yeo-Johnson,
#'   além de seleção automatizada do melhor método por teste de normalidade.
#'
#' @param dados Vetor numérico ou \code{data.frame} com variáveis a transformar.
#' @param metodo Método de transformação: \code{"auto"} (seleção automática pelo
#'   maior W de Shapiro-Wilk), \code{"box_cox"}, \code{"yeo_johnson"}, \code{"log"},
#'   \code{"log10"}, \code{"raiz_quadrada"}, \code{"raiz_cubica"}, \code{"arco_seno"},
#'   \code{"logit"}, \code{"inversa"}, \code{"padronizacao"} (escore Z) ou
#'   \code{"minmax"} (escala 0-1).
#' @param lambda Valor de potência fixado manualmente (opcional, para Box-Cox ou Yeo-Johnson).
#' @param colunas Vetor de nomes de colunas a transformar quando \code{dados} for um \code{data.frame}.
#'   Se \code{NULL} (padrão), transforma todas as colunas numéricas.
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o gráfico apropriado: curva de
#'   perfil de verossimilhança para Box-Cox ou painel comparativo Antes/Depois para os demais.
#'
#' @return Retorna um objeto da classe \code{estatR_transformacao} contendo:
#' \itemize{
#'   \item \code{dados_transformados}: Vetor ou \code{data.frame} com os dados transformados.
#'   \item \code{valores}: Alias direto para \code{dados_transformados}.
#'   \item \code{tabela}: Tabela resumo contendo método, parâmetros e ganho de normalidade.
#'   \item \code{lambda}: Valor de \eqn{\lambda} estimado ou utilizado (quando aplicável).
#'   \item \code{metodo}: Nome do método aplicado.
#'   \item \code{grafico}: Objeto gráfico do \code{ggplot2}.
#' }
#' @export
transformar <- function(dados,
                        metodo  = c("auto", "box_cox", "yeo_johnson", "log", "log10",
                                    "raiz_quadrada", "raiz_cubica", "arco_seno",
                                    "logit", "inversa", "padronizacao", "minmax"),
                        lambda  = NULL,
                        colunas = NULL,
                        grafico = TRUE) {

  metodo <- match.arg(metodo)

  # ── Caso 1: Entrada como Data Frame ──────────────────────────────────────────
  if (is.data.frame(dados)) {
    return(.transformar_df(dados, metodo = metodo, lambda = lambda,
                           colunas = colunas, grafico = grafico))
  }

  # ── Caso 2: Entrada como Vetor ───────────────────────────────────────────────
  if (!is.numeric(dados)) {
    stop("O argumento 'dados' deve ser um vetor num\u00e9rico ou data.frame.", call. = FALSE)
  }

  nome_var <- deparse(substitute(dados))
  nome_var <- sub(".*\\$", "", nome_var)
  if (nchar(nome_var) > 25 || grepl("c\\(", nome_var)) nome_var <- "Variavel"

  res_vetor <- .transformar_vetor(dados, metodo = metodo, lambda = lambda,
                                  nome_var = nome_var, grafico = grafico)
  return(invisible(res_vetor))
}


# ── Função interna para transformar vetor ──────────────────────────────────────
.transformar_vetor <- function(x, metodo, lambda = NULL, nome_var = "Y", grafico = TRUE) {
  x_orig <- x
  idx_val <- !is.na(x)
  y <- x[idx_val]
  n_total <- length(x)
  n_validos <- length(y)

  if (n_validos < 4) {
    stop("S\u00e3o necess\u00e1rias pelo menos 4 observa\u00e7\u00f5es v\u00e1lidas.", call. = FALSE)
  }

  # Teste Shapiro antes
  sw_antes <- tryCatch(stats::shapiro.test(y), error = function(e) list(statistic = NA, p.value = NA))

  lambda_usado <- NULL
  ic_lambda <- NULL
  perfil_bc <- NULL
  formula_txt <- ""
  metodo_aplicado <- metodo

  # ── Seleção Automática ───────────────────────────────────────────────────────
  if (metodo == "auto") {
    candidatos <- c("box_cox", "yeo_johnson", "log", "raiz_quadrada", "inversa")
    melhor_w <- if (!is.na(sw_antes$statistic)) sw_antes$statistic else -Inf
    melhor_m <- "original"
    melhor_y <- y
    melhor_lam <- NULL

    for (cand in candidatos) {
      tentativa <- tryCatch(.executar_transf(y, cand, lambda = NULL), error = function(e) NULL)
      if (!is.null(tentativa)) {
        sw_t <- tryCatch(stats::shapiro.test(tentativa$y_transf), error = function(e) NULL)
        if (!is.null(sw_t) && !is.na(sw_t$statistic) && sw_t$statistic > melhor_w) {
          melhor_w <- sw_t$statistic
          melhor_m <- cand
          melhor_y <- tentativa$y_transf
          melhor_lam <- tentativa$lambda
        }
      }
    }

    if (melhor_m == "original") {
      metodo_aplicado <- "original"
      y_transf <- y
      formula_txt <- "y (sem altera\u00e7\u00e3o)"
    } else {
      metodo_aplicado <- melhor_m
      lambda_usado <- melhor_lam
      exec <- .executar_transf(y, metodo_aplicado, lambda = lambda_usado)
      y_transf <- exec$y_transf
      formula_txt <- exec$formula_txt
      perfil_bc <- exec$perfil_bc
      ic_lambda <- exec$ic_lambda
    }
  } else {
    exec <- .executar_transf(y, metodo, lambda = lambda)
    y_transf <- exec$y_transf
    lambda_usado <- exec$lambda
    formula_txt <- exec$formula_txt
    perfil_bc <- exec$perfil_bc
    ic_lambda <- exec$ic_lambda
  }

  # Teste Shapiro depois
  sw_depois <- tryCatch(stats::shapiro.test(y_transf), error = function(e) list(statistic = NA, p.value = NA))

  # Montar vetor de saída preservando NAs
  res_vec <- rep(NA_real_, n_total)
  res_vec[idx_val] <- y_transf

  # ── Tabela de Diagnóstico ────────────────────────────────────────────────────
  tb_resumo <- data.frame(
    Variavel   = nome_var,
    Metodo     = metodo_aplicado,
    Expressao  = formula_txt,
    W_Antes    = as.numeric(sw_antes$statistic),
    p_Antes    = as.numeric(sw_antes$p.value),
    W_Depois   = as.numeric(sw_depois$statistic),
    p_Depois   = as.numeric(sw_depois$p.value),
    stringsAsFactors = FALSE
  )

  # ── Saída no Console ─────────────────────────────────────────────────────────
  .print_titulo("TRANSFORMA\u00c7\u00c3O DE VARI\u00c1VEIS")

  cat(sprintf("  Vari\u00e1vel: %s  |  Observa\u00e7\u00f5es: %d (V\u00e1lidas: %d)\n", nome_var, n_total, n_validos))
  cat(sprintf("  M\u00e9todo Selecionado: %s\n", toupper(metodo_aplicado)))
  cat(sprintf("  Express\u00e3o Aplicada:  y* = %s\n\n", formula_txt))

  if (!is.null(lambda_usado) && !is.null(ic_lambda)) {
    cat(sprintf("  Estimativa de \u03bb (\u00d3timo):  %.3f\n", lambda_usado))
    cat(sprintf("  Intervalo de Confian\u00e7a 95%%:   [%.3f; %.3f]\n\n", ic_lambda[1], ic_lambda[2]))
  }

  .print_topico("AVALIA\u00c7\u00c3O DA ADER\u00caNCIA \u00c0 NORMALIDADE (SHAPIRO-WILK)")

  w_cols <- c(18, 14, 14, 14, 14)
  hdr_txt <- paste0(
    "  ",
    .th_pad("Momento",       w_cols[1], "left"),
    .th_pad("Estat\u00edstica W", w_cols[2], "center"),
    .th_pad("p-valor",       w_cols[3], "center"),
    .th_pad("M\u00e9dia",        w_cols[4], "center"),
    .th_pad("Desvio-Padr\u00e3o", w_cols[5], "center")
  )
  sep_line <- paste0("  ", strrep("\u2500", sum(w_cols)))

  cat(sep_line, "\n")
  cat(hdr_txt, "\n")
  cat(sep_line, "\n")

  cat(
    "  ",
    .th_pad("Original (Antes)", w_cols[1], "left"),
    .th_pad(formatC(sw_antes$statistic, format = "f", digits = 4, decimal.mark = ","), w_cols[2], "center"),
    .th_pad(formatC(sw_antes$p.value, format = "f", digits = 4, decimal.mark = ","), w_cols[3], "center"),
    .th_pad(formatC(mean(y), format = "f", digits = 2, decimal.mark = ","), w_cols[4], "center"),
    .th_pad(formatC(stats::sd(y), format = "f", digits = 2, decimal.mark = ","), w_cols[5], "center"),
    "\n", sep = ""
  )

  cat(
    "  ",
    .th_pad("Transformada",     w_cols[1], "left"),
    .th_pad(formatC(sw_depois$statistic, format = "f", digits = 4, decimal.mark = ","), w_cols[2], "center"),
    .th_pad(formatC(sw_depois$p.value, format = "f", digits = 4, decimal.mark = ","), w_cols[3], "center"),
    .th_pad(formatC(mean(y_transf), format = "f", digits = 2, decimal.mark = ","), w_cols[4], "center"),
    .th_pad(formatC(stats::sd(y_transf), format = "f", digits = 2, decimal.mark = ","), w_cols[5], "center"),
    "\n", sep = ""
  )
  cat(sep_line, "\n\n")

  # Conclusão pedagógica
  ganho <- as.numeric(sw_depois$statistic) - as.numeric(sw_antes$statistic)
  if (!is.na(ganho) && ganho > 0) {
    cat(sprintf("  \u2713 A transforma\u00e7\u00e3o aumentou o coeficiente W em +%.4f, aproximando a distribui\u00e7\u00e3o da Normal.\n", ganho))
  } else if (!is.na(ganho) && abs(ganho) <= 0.001) {
    cat("  * A distribui\u00e7\u00e3o permaneceu estatisticamente est\u00e1vel ap\u00f3s o procedimento.\n")
  }
  .print_rodape()

  # ── Geração do Gráfico ───────────────────────────────────────────────────────
  p_out <- NULL
  if (grafico) {
    if (metodo_aplicado == "box_cox" && !is.null(perfil_bc)) {
      # Gráfico de perfil de log-verossimilhança com bandas de confiança
      df_bc <- perfil_bc$grid
      l_opt <- lambda_usado
      ic_inf <- ic_lambda[1]
      ic_sup <- ic_lambda[2]
      lim_y <- perfil_bc$lim_loglik

      p_out <- ggplot2::ggplot(df_bc, ggplot2::aes(x = lambda, y = loglik)) +
        # Faixa de confiança horizontal e vertical
        ggplot2::annotate("rect", xmin = ic_inf, xmax = ic_sup,
                          ymin = -Inf, ymax = lim_y,
                          fill = "#85b4d4", alpha = 0.25) +
        ggplot2::geom_line(color = "#2B5C8F", linewidth = 1.1) +
        ggplot2::geom_hline(yintercept = lim_y, linetype = "dashed", color = "#555555", linewidth = 0.7) +
        ggplot2::geom_vline(xintercept = l_opt, color = "#D90429", linewidth = 1.1) +
        ggplot2::geom_vline(xintercept = c(ic_inf, ic_sup), color = "#85b4d4", linetype = "dotted", linewidth = 0.8) +
        ggplot2::annotate("point", x = l_opt, y = max(df_bc$loglik), color = "#D90429", size = 3) +
        ggplot2::annotate("text", x = l_opt, y = max(df_bc$loglik),
                          label = sprintf("  \u03bb = %.2f", l_opt),
                          color = "#D90429", fontface = "bold", hjust = -0.1, vjust = 0) +
        ggplot2::labs(
          title    = sprintf("Perfil de Log-Verossimilhan\u00e7a de Box-Cox (%s)", nome_var),
          subtitle = sprintf("\u03bb \u00d3timo: %.2f  |  IC 95%%: [%.2f; %.2f]", l_opt, ic_inf, ic_sup),
          x = "Par\u00e2metro de Pot\u00eancia (\u03bb)",
          y = "Log-Verossimilhan\u00e7a"
        )
    } else {
      # Painel comparativo Antes vs Depois (Densidades)
      df_comp <- data.frame(
        Valor = c(y, y_transf),
        Etapa = factor(rep(c("Antes (Original)", "Depois (Transformada)"), each = length(y)),
                       levels = c("Antes (Original)", "Depois (Transformada)"))
      )

      p_out <- ggplot2::ggplot(df_comp, ggplot2::aes(x = Valor, fill = Etapa)) +
        ggplot2::geom_density(alpha = 0.45, color = "#2B5C8F", linewidth = 0.8) +
        ggplot2::scale_fill_manual(values = c("Antes (Original)" = "#D90429", "Depois (Transformada)" = "#2B5C8F")) +
        ggplot2::facet_wrap(~ Etapa, scales = "free") +
        ggplot2::labs(
          title    = sprintf("Compara\u00e7\u00e3o da Distribui\u00e7\u00e3o: %s", nome_var),
          subtitle = sprintf("M\u00e9todo: %s  |  y* = %s", toupper(metodo_aplicado), formula_txt),
          x = "Valores", y = "Densidade Te\u00f3rica"
        ) +
        ggplot2::theme(legend.position = "none")
    }

    if (exists("meu_tema", mode = "function")) {
      p_out <- p_out + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p_out <- p_out + tema_estatR(estilo = 1, grade = "dupla")
    } else {
      p_out <- p_out + ggplot2::theme_minimal()
    }

    suppressMessages(print(p_out))
  }

  res <- list(
    dados_transformados = res_vec,
    valores             = res_vec,
    tabela              = tb_resumo,
    lambda              = lambda_usado,
    metodo              = metodo_aplicado,
    grafico             = p_out
  )
  class(res) <- c("estatR_transformacao", "list")
  return(invisible(res))
}


# ── Execução matemática das transformações ─────────────────────────────────────
.executar_transf <- function(y, metodo, lambda = NULL) {
  perfil_bc <- NULL
  ic_lambda <- NULL
  formula_txt <- ""
  y_transf <- NULL

  if (metodo == "box_cox") {
    # Requer positivos
    shift <- 0
    if (min(y) <= 0) {
      shift <- abs(min(y)) + 1
      y_est <- y + shift
    } else {
      y_est <- y
    }

    if (is.null(lambda)) {
      grid_lam <- seq(-2.5, 2.5, by = 0.02)
      n <- length(y_est)
      geom_mean <- exp(mean(log(y_est)))

      llik <- sapply(grid_lam, function(l) {
        yt <- if (abs(l) < 1e-4) log(y_est) else (y_est^l - 1) / l
        s2 <- mean((yt - mean(yt))^2)
        if (s2 <= 0) return(-Inf)
        -0.5 * n * log(s2) + (l - 1) * sum(log(y_est))
      })

      idx_max <- which.max(llik)
      lambda_opt <- grid_lam[idx_max]
      max_llik <- llik[idx_max]
      lim_llik <- max_llik - 0.5 * stats::qchisq(0.95, df = 1)

      ic_vals <- grid_lam[llik >= lim_llik]
      ic_lambda <- c(min(ic_vals), max(ic_vals))

      perfil_bc <- list(
        grid = data.frame(lambda = grid_lam, loglik = llik),
        lim_loglik = lim_llik
      )
      lambda <- round(lambda_opt, 2)
    }

    if (abs(lambda) < 1e-3) {
      y_transf <- log(y_est)
      formula_txt <- if (shift == 0) "ln(y)" else sprintf("ln(y + %.1f)", shift)
    } else {
      y_transf <- (y_est^lambda - 1) / lambda
      formula_txt <- if (shift == 0) sprintf("(y^%.2f - 1) / %.2f", lambda, lambda) else sprintf("((y + %.1f)^%.2f - 1) / %.2f", shift, lambda, lambda)
    }

  } else if (metodo == "yeo_johnson") {
    # Aceita qualquer valor real
    if (is.null(lambda)) {
      grid_lam <- seq(-2.5, 2.5, by = 0.02)
      n <- length(y)
      llik <- sapply(grid_lam, function(l) {
        yt <- .calc_yj(y, l)
        s2 <- mean((yt - mean(yt))^2)
        if (s2 <= 0) return(-Inf)
        -0.5 * n * log(s2) + (l - 1) * sum(sign(y) * log(abs(y) + 1))
      })
      lambda <- round(grid_lam[which.max(llik)], 2)
    }
    y_transf <- .calc_yj(y, lambda)
    formula_txt <- sprintf("Yeo-Johnson(\u03bb = %.2f)", lambda)

  } else if (metodo == "log") {
    if (min(y) <= 0) {
      c_add <- abs(min(y)) + 1
      y_transf <- log(y + c_add)
      formula_txt <- sprintf("ln(y + %.2f)", c_add)
    } else {
      y_transf <- log(y)
      formula_txt <- "ln(y)"
    }

  } else if (metodo == "log10") {
    if (min(y) <= 0) {
      c_add <- abs(min(y)) + 1
      y_transf <- log10(y + c_add)
      formula_txt <- sprintf("log10(y + %.2f)", c_add)
    } else {
      y_transf <- log10(y)
      formula_txt <- "log10(y)"
    }

  } else if (metodo == "raiz_quadrada") {
    if (min(y) < 0) {
      c_add <- abs(min(y))
      y_transf <- sqrt(y + c_add)
      formula_txt <- sprintf("\u221a(y + %.2f)", c_add)
    } else {
      y_transf <- sqrt(y)
      formula_txt <- "\u221a(y)"
    }

  } else if (metodo == "raiz_cubica") {
    y_transf <- sign(y) * (abs(y)^(1/3))
    formula_txt <- "\u221b(y)"

  } else if (metodo == "arco_seno") {
    if (min(y) < 0 || max(y) > 1) {
      stop("A transforma\u00e7\u00e3o arco-seno exige dados estritamente em [0, 1].", call. = FALSE)
    }
    y_transf <- asin(sqrt(y))
    formula_txt <- "arcsen(\u221ay)"

  } else if (metodo == "logit") {
    if (min(y) <= 0 || max(y) >= 1) {
      stop("A transforma\u00e7\u00e3o logit exige dados estritamente em (0, 1).", call. = FALSE)
    }
    y_transf <- log(y / (1 - y))
    formula_txt <- "ln(y / (1 - y))"

  } else if (metodo == "inversa") {
    if (any(y == 0)) {
      c_add <- 1
      y_transf <- 1 / (y + c_add)
      formula_txt <- "1 / (y + 1)"
    } else {
      y_transf <- 1 / y
      formula_txt <- "1 / y"
    }

  } else if (metodo == "padronizacao") {
    s_val <- stats::sd(y)
    if (s_val == 0) stop("Vari\u00e2ncia nula para padroniza\u00e7\u00e3o.", call. = FALSE)
    y_transf <- (y - mean(y)) / s_val
    formula_txt <- "(y - \u03bc) / \u03c3"

  } else if (metodo == "minmax") {
    rng <- max(y) - min(y)
    if (rng == 0) stop("Amplitude nula para min-max.", call. = FALSE)
    y_transf <- (y - min(y)) / rng
    formula_txt <- "(y - min) / (max - min)"
  }

  return(list(
    y_transf    = y_transf,
    lambda      = lambda,
    formula_txt = formula_txt,
    perfil_bc   = perfil_bc,
    ic_lambda   = ic_lambda
  ))
}

.calc_yj <- function(y, l) {
  res <- numeric(length(y))
  pos <- y >= 0
  neg <- !pos

  if (any(pos)) {
    yp <- y[pos]
    if (abs(l) < 1e-4) {
      res[pos] <- log(yp + 1)
    } else {
      res[pos] <- ((yp + 1)^l - 1) / l
    }
  }

  if (any(neg)) {
    yn <- y[neg]
    l2 <- 2 - l
    if (abs(l2) < 1e-4) {
      res[neg] <- -log(-yn + 1)
    } else {
      res[neg] <- -((-yn + 1)^l2 - 1) / l2
    }
  }
  return(res)
}


# ── Caso Data Frame ────────────────────────────────────────────────────────────
.transformar_df <- function(dados, metodo, lambda, colunas, grafico) {
  df_transf <- dados
  num_cols <- names(dados)[sapply(dados, is.numeric)]

  if (!is.null(colunas)) {
    num_cols <- intersect(colunas, num_cols)
    if (length(num_cols) == 0) stop("Nenhuma coluna num\u00e9rica v\u00e1lida selecionada.", call. = FALSE)
  }

  tabelas_list <- list()
  lista_plots <- list()

  for (col in num_cols) {
    vec <- dados[[col]]
    res_v <- .transformar_vetor(vec, metodo = metodo, lambda = lambda, nome_var = col, grafico = FALSE)
    df_transf[[col]] <- res_v$dados_transformados
    tabelas_list[[col]] <- res_v$tabela
    if (!is.null(res_v$grafico)) lista_plots[[col]] <- res_v$grafico
  }

  tb_geral <- do.call(rbind, tabelas_list)
  rownames(tb_geral) <- NULL

  p_panel <- NULL
  if (grafico && length(num_cols) > 0) {
    # Gráfico facetado de densidade para todas as variáveis transformadas
    long_list <- list()
    for (col in num_cols) {
      yo <- stats::na.omit(dados[[col]])
      yt <- stats::na.omit(df_transf[[col]])
      long_list[[col]] <- data.frame(
        Variavel = col,
        Valor    = c(yo, yt),
        Etapa    = factor(rep(c("Antes", "Depois"), times = c(length(yo), length(yt))),
                          levels = c("Antes", "Depois"))
      )
    }
    df_long <- do.call(rbind, long_list)

    p_panel <- ggplot2::ggplot(df_long, ggplot2::aes(x = Valor, fill = Etapa)) +
      ggplot2::geom_density(alpha = 0.45, color = "#2B5C8F", linewidth = 0.7) +
      ggplot2::scale_fill_manual(values = c("Antes" = "#D90429", "Depois" = "#2B5C8F")) +
      ggplot2::facet_wrap(~ Variavel, scales = "free") +
      ggplot2::labs(
        title    = "Compara\u00e7\u00e3o das Distribui\u00e7\u00f5es: Banco de Dados",
        subtitle = sprintf("M\u00e9todo Aplicado: %s", toupper(metodo)),
        x = "Valores das Vari\u00e1veis", y = "Densidade Te\u00f3rica"
      )

    if (exists("meu_tema", mode = "function")) {
      p_panel <- p_panel + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p_panel <- p_panel + tema_estatR(estilo = 1, grade = "dupla")
    } else {
      p_panel <- p_panel + ggplot2::theme_minimal()
    }
    suppressMessages(print(p_panel))
  }

  res <- list(
    dados_transformados = df_transf,
    valores             = df_transf,
    tabela              = tb_geral,
    lambda              = if (length(num_cols) == 1) tb_geral$Expressao else NULL,
    metodo              = metodo,
    grafico             = p_panel
  )
  class(res) <- c("estatR_transformacao", "list")
  return(invisible(res))
}

#' @export
print.estatR_transformacao <- function(x, ...) {
  if (!is.null(x$tabela)) {
    # Se chamado diretamente no print
    invisible(x)
  }
}
