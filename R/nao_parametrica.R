# ==============================================================================
# MÓDULO: ESTATÍSTICA NÃO-PARAMÉTRICA & CONCORDÂNCIA
# Kruskal-Wallis (Dunn), Friedman, Concordância Kappa (Cohen/Fleiss) e Bland-Altman
# ==============================================================================

# ── 1. TESTE DE KRUSKAL-WALLIS COM POST-HOC DE DUNN ─────────────────────────

#' @title Teste de Kruskal-Wallis com Comparações Múltiplas (Dunn)
#' @description Alternativa não-paramétrica à ANOVA a um fator para comparar a
#'   distribuição de \eqn{k \ge 2}{k >= 2} grupos independentes. Caso haja significância estatística,
#'   executa automaticamente o teste post-hoc de Dunn com ajuste de p-valor
#'   (Holm, Bonferroni ou FDR) e gera gráfico de boxplot padronizado.
#'
#' @param dados Data frame contendo os dados.
#' @param formula Fórmula no formato \code{y ~ grupo}.
#' @param y Variável numérica resposta (símbolo ou string), se fórmula não fornecida.
#' @param grupo Variável categórica / fator de agrupamento, se fórmula não fornecida.
#' @param post_hoc Método de comparações múltiplas caso \eqn{p < 0{,}05}{p < 0.05}: \code{"dunn"} (padrão),
#'   \code{"bonferroni"} ou \code{"nenhum"}.
#' @param ajuste_p Método de correção para múltiplos testes: \code{"holm"} (padrão),
#'   \code{"bonferroni"} ou \code{"fdr"}.
#' @param confianca Nível de confiança (padrão: 0.95).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o boxplot com medianas.
#'
#' @return Objeto invisível com estatísticas globais e comparações múltiplas.
#' @export
teste_kruskal <- function(dados, formula = NULL, y = NULL, grupo = NULL,
                          post_hoc = c("dunn", "bonferroni", "nenhum"),
                          ajuste_p = c("holm", "bonferroni", "fdr"),
                          confianca = 0.95, grafico = TRUE) {

  post_hoc <- match.arg(post_hoc)
  ajuste_p <- match.arg(ajuste_p)

  # Resolução de argumentos: Fórmula vs colunas
  if (!is.null(formula)) {
    if (!inherits(formula, "formula")) stop("O argumento 'formula' deve ser uma f\u00f3rmula v\u00e1lida (ex: y ~ grupo).")
    vars <- all.vars(formula)
    if (length(vars) < 2) stop("A f\u00f3rmula deve conter vari\u00e1vel resposta e fator (ex: y ~ grupo).")
    y_name <- vars[1]
    g_name <- vars[2]
  } else {
    y_sub <- substitute(y)
    g_sub <- substitute(grupo)
    y_name <- if (is.symbol(y_sub)) as.character(y_sub) else if (is.character(y)) y else NULL
    g_name <- if (is.symbol(g_sub)) as.character(g_sub) else if (is.character(grupo)) grupo else NULL

    if (is.null(y_name) || is.null(g_name)) {
      stop("Forne\u00e7a uma 'formula' (y ~ grupo) ou especifique 'y' e 'grupo'.")
    }
  }

  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  if (!y_name %in% names(dados)) stop(sprintf("Vari\u00e1vel resposta '%s' n\u00e3o encontrada nos dados.", y_name))
  if (!g_name %in% names(dados)) stop(sprintf("Vari\u00e1vel de grupo '%s' n\u00e3o encontrada nos dados.", g_name))

  df_clean <- dados[!is.na(dados[[y_name]]) & !is.na(dados[[g_name]]), c(y_name, g_name), drop = FALSE]
  y_vec <- as.numeric(df_clean[[y_name]])
  g_vec <- as.factor(df_clean[[g_name]])

  n_total <- length(y_vec)
  niveis <- levels(g_vec)
  k <- length(niveis)

  if (k < 2) stop("O fator de agrupamento deve conter pelo menos 2 n\u00edveis.")

  # 1. Teste de Kruskal-Wallis oficial via stats::kruskal.test
  kw_res <- stats::kruskal.test(y_vec ~ g_vec)
  stat_h <- as.numeric(kw_res$statistic)
  df_h   <- as.numeric(kw_res$parameter)
  p_val  <- as.numeric(kw_res$p.value)

  # Formatação p-valor
  fmt_p <- function(p, com_sinal = FALSE) {
    if (is.na(p)) return("-")
    if (p < 0.001) return(if (com_sinal) "< 0,001" else "< 0,001")
    val <- formatC(p, format = "f", digits = 3, decimal.mark = ",")
    if (com_sinal) paste0("= ", val) else val
  }

  # Tabela global (SEM coluna de decisão, conforme solicitado)
  df_global <- data.frame(
    "Estat\u00edstica H" = .fmt_num(stat_h, 2),
    "Graus de Liberdade" = as.character(df_h),
    "p-valor" = fmt_p(p_val),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # 2. Comparações Múltiplas Post-Hoc (Dunn ou Bonferroni)
  df_posthoc <- NULL
  if (p_val < (1 - confianca) && post_hoc != "nenhum" && k > 2) {
    # Cálculo manual do teste de Dunn para manter independência e robustez sem pacotes externos
    ranks <- rank(y_vec)
    n_i <- table(g_vec)
    mean_ranks <- tapply(ranks, g_vec, mean)

    # Fator de correção para empates
    ties <- table(ranks)
    C_ties <- 1 - sum(ties^3 - ties) / (n_total^3 - n_total)
    var_base <- (n_total * (n_total + 1) / 12) * C_ties

    pares <- utils::combn(niveis, 2, simplify = FALSE)
    n_comp <- length(pares)

    diff_vet <- numeric(n_comp)
    z_vet    <- numeric(n_comp)
    p_raw    <- numeric(n_comp)
    rotulos_pares <- character(n_comp)

    for (i in seq_along(pares)) {
      g1 <- pares[[i]][1]
      g2 <- pares[[i]][2]
      rotulos_pares[i] <- paste(g1, "vs", g2)

      diff_val <- mean_ranks[g1] - mean_ranks[g2]
      diff_vet[i] <- abs(diff_val)

      se_diff <- sqrt(var_base * (1 / n_i[g1] + 1 / n_i[g2]))
      z_val <- diff_val / se_diff
      z_vet[i] <- abs(z_val)
      p_raw[i] <- 2 * stats::pnorm(-abs(z_val))
    }

    p_adj <- stats::p.adjust(p_raw, method = ajuste_p)

    sig_ast <- sapply(p_adj, function(p) {
      if (is.na(p) || p >= (1 - confianca)) "" else "*"
    })

    df_posthoc <- data.frame(
      "Compara\u00e7\u00e3o" = rotulos_pares,
      "Diferen\u00e7a Postos" = sapply(diff_vet, .fmt_num, decimais = 2),
      "Estat\u00edstica Z" = sapply(z_vet, .fmt_num, decimais = 2),
      "p-ajustado" = sapply(p_adj, fmt_p),
      "Sig." = sig_ast,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("TESTE DE KRUSKAL-WALLIS")

  cat(sprintf("  Vari\u00e1vel Resposta: %s\n", y_name))
  cat(sprintf("  Fator / Grupos:     %s (%d grupos: %s)\n", g_name, k, paste(niveis, collapse = ", ")))
  cat(sprintf("  Observa\u00e7\u00f5es:        %d v\u00e1lidas\n\n", n_total))

  .print_topico("RESULTADO DO TESTE GLOBAL")

  if (exists(".print_tabela_estatR", mode = "function")) {
    # Todas as colunas centralizadas rigorosamente
    .print_tabela_estatR(df_global, align = c("center", "center", "center"))
  } else {
    print(df_global, row.names = FALSE)
  }
  cat("\n")

  if (!is.null(df_posthoc)) {
    .print_topico(sprintf("COMPARA\u00c7\u00d5ES M\u00daLTIPLAS (TESTE DE DUNN - AJUSTE %s)", toupper(ajuste_p)))
    if (exists(".print_tabela_estatR", mode = "function")) {
      .print_tabela_estatR(df_posthoc, align = c("center", "center", "center", "center", "center"))
    } else {
      print(df_posthoc, row.names = FALSE)
    }
    cat("  * Diferen\u00e7a estatisticamente significativa ao n\u00edvel de signific\u00e2ncia de 5%.\n\n")
  }

  .print_topico("INTERPRETA\u00c7\u00c3O")

  cor_ok   <- "\033[1;32m"
  cor_erro <- "\033[1;31m"
  reset    <- "\033[0m"

  if (p_val < (1 - confianca)) {
    cat(sprintf("  %s[!]%s Diferen\u00e7a estatisticamente significativa detectada entre os grupos (H = %s; p %s).\n",
                cor_erro, reset, .fmt_num(stat_h, 2), fmt_p(p_val, com_sinal = TRUE)))
    cat("      Rejeita-se a hip\u00f3tese nula de igualdade: as distribui\u00e7\u00f5es / medianas dos grupos diferem.\n")
    if (!is.null(df_posthoc)) {
      cat("      Consulte a tabela de compara\u00e7\u00f5es m\u00faltiplas acima para identificar os pares com diferen\u00e7a.\n")
    }
  } else {
    cat(sprintf("  %s[\u2713]%s N\u00e3o h\u00e1 evid\u00eancia de diferen\u00e7a significativa entre os grupos (H = %s; p %s).\n",
                cor_ok, reset, .fmt_num(stat_h, 2), fmt_p(p_val, com_sinal = TRUE)))
    cat("      N\u00e3o se rejeita a hip\u00f3tese nula: as medianas populacionais s\u00e3o consideradas equivalentes.\n")
  }

  .print_rodape()

  # ── Gráfico Padronizado ────────────────────────────────────────────────────
  p_plot <- NULL
  if (isTRUE(grafico)) {
    df_plot <- data.frame(
      .y = y_vec,
      .g = g_vec
    )

    p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .g, y = .y, fill = .g)) +
      ggplot2::geom_boxplot(alpha = 0.75, width = 0.5, outlier.shape = 16, outlier.size = 1.8, outlier.alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 1.5, color = "gray25") +
      ggplot2::stat_summary(fun = stats::median, geom = "point", shape = 18, size = 3, color = "#E31A1C") +
      ggplot2::labs(
        title = paste("Distribui\u00e7\u00e3o de", y_name, "por", g_name),
        subtitle = sprintf("Kruskal-Wallis: H = %s | p = %s", .fmt_num(stat_h, 2), fmt_p(p_val)),
        x = g_name,
        y = y_name
      ) +
      ggplot2::guides(fill = "none") +
      tema_estatR(estilo = 2, grade = "dupla")

    print(p_plot)
  }

  invisible(list(
    teste_global = kw_res,
    post_hoc = df_posthoc,
    grafico = p_plot
  ))
}


# ── 2. TESTE DE FRIEDMAN COM COMPARAÇÕES PAREADAS ────────────────────────────

#' @title Teste de Friedman para Medidas Repetidas / Blocos
#' @description Alternativa não-paramétrica à ANOVA de medidas repetidas ou blocos casualizados (DBC).
#'   Avalia se as distribuições diferem entre \eqn{k \ge 3}{k >= 3} condições ou períodos pareados.
#'   Caso significativo, executa comparações pareadas com ajuste de p-valor.
#'
#' @param dados Data frame contendo os dados.
#' @param formula Fórmula no formato \code{y ~ grupo | bloco}.
#' @param y Variável numérica resposta, se fórmula não fornecida.
#' @param grupo Fator de grupos / tratamentos / tempos, se fórmula não fornecida.
#' @param bloco Fator de identificação do sujeito / bloco, se fórmula não fornecida.
#' @param post_hoc Lógico. Se \code{TRUE} (padrão), executa comparações pareadas caso \eqn{p < 0{,}05}{p < 0.05}.
#' @param ajuste_p Método de ajuste de p-valor: \code{"holm"} (padrão), \code{"bonferroni"} ou \code{"fdr"}.
#' @param confianca Nível de confiança (padrão: 0.95).
#' @param grafico Lógico. Se \code{TRUE} (padrão), gera gráfico de perfil temporal / tratamentos.
#'
#' @return Objeto invisível com o teste global e comparações pareadas.
#' @export
teste_friedman <- function(dados, formula = NULL, y = NULL, grupo = NULL, bloco = NULL,
                           post_hoc = TRUE, ajuste_p = c("holm", "bonferroni", "fdr"),
                           confianca = 0.95, grafico = TRUE) {

  ajuste_p <- match.arg(ajuste_p)

  # Resolução de argumentos
  if (!is.null(formula)) {
    if (!inherits(formula, "formula")) stop("O argumento 'formula' deve ser no formato y ~ grupo | bloco.")
    # Extração de y ~ grupo | bloco
    f_str <- deparse(formula)
    if (!grepl("\\|", f_str)) {
      stop("A f\u00f3rmula para o teste de Friedman deve conter o bloco no formato y ~ grupo | bloco.")
    }
    partes <- strsplit(f_str, "\\|")[[1]]
    esquerda_dir <- strsplit(partes[1], "~")[[1]]
    y_name <- trimws(esquerda_dir[1])
    g_name <- trimws(esquerda_dir[2])
    b_name <- trimws(partes[2])
  } else {
    y_sub <- substitute(y)
    g_sub <- substitute(grupo)
    b_sub <- substitute(bloco)
    y_name <- if (is.symbol(y_sub)) as.character(y_sub) else if (is.character(y)) y else NULL
    g_name <- if (is.symbol(g_sub)) as.character(g_sub) else if (is.character(grupo)) grupo else NULL
    b_name <- if (is.symbol(b_sub)) as.character(b_sub) else if (is.character(bloco)) bloco else NULL

    if (is.null(y_name) || is.null(g_name) || is.null(b_name)) {
      stop("Forne\u00e7a 'formula' (y ~ grupo | bloco) ou informe 'y', 'grupo' e 'bloco'.")
    }
  }

  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  if (!y_name %in% names(dados)) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", y_name))
  if (!g_name %in% names(dados)) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", g_name))
  if (!b_name %in% names(dados)) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", b_name))

  df_clean <- dados[!is.na(dados[[y_name]]) & !is.na(dados[[g_name]]) & !is.na(dados[[b_name]]), , drop = FALSE]

  # Teste oficial via stats::friedman.test
  f_formula <- stats::as.formula(sprintf("%s ~ %s | %s", y_name, g_name, b_name))
  fr_res <- stats::friedman.test(f_formula, data = df_clean)

  stat_chi <- as.numeric(fr_res$statistic)
  df_chi   <- as.numeric(fr_res$parameter)
  p_val    <- as.numeric(fr_res$p.value)

  fmt_p <- function(p, com_sinal = FALSE) {
    if (is.na(p)) return("-")
    if (p < 0.001) return(if (com_sinal) "< 0,001" else "< 0,001")
    val <- formatC(p, format = "f", digits = 3, decimal.mark = ",")
    if (com_sinal) paste0("= ", val) else val
  }

  # Tabela global (SEM coluna de decisão)
  df_global <- data.frame(
    "Qui-Quadrado (\u03c7\u00b2)" = .fmt_num(stat_chi, 2),
    "Graus de Liberdade" = as.character(df_chi),
    "p-valor" = fmt_p(p_val),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  niveis_g <- levels(as.factor(df_clean[[g_name]]))
  n_blocos <- length(unique(df_clean[[b_name]]))

  # Pós-teste de comparações pareadas via Wilcoxon pareado
  df_posthoc <- NULL
  if (p_val < (1 - confianca) && isTRUE(post_hoc) && length(niveis_g) > 2) {
    pares <- utils::combn(niveis_g, 2, simplify = FALSE)
    n_comp <- length(pares)

    rotulos_pares <- character(n_comp)
    stat_v        <- numeric(n_comp)
    p_raw         <- numeric(n_comp)

    # Reestrutura em matriz bloco x grupo
    mat_wide <- stats::reshape(df_clean[, c(y_name, g_name, b_name)],
                               idvar = b_name, timevar = g_name, direction = "wide")

    for (i in seq_along(pares)) {
      g1 <- pares[[i]][1]
      g2 <- pares[[i]][2]
      rotulos_pares[i] <- paste(g1, "vs", g2)

      col1 <- mat_wide[[paste0(y_name, ".", g1)]]
      col2 <- mat_wide[[paste0(y_name, ".", g2)]]

      wt <- tryCatch(
        suppressWarnings(stats::wilcox.test(col1, col2, paired = TRUE, exact = FALSE)),
        error = function(e) list(statistic = NA_real_, p.value = NA_real_)
      )

      stat_v[i] <- as.numeric(wt$statistic)
      p_raw[i]  <- as.numeric(wt$p.value)
    }

    p_adj <- stats::p.adjust(p_raw, method = ajuste_p)

    sig_ast <- sapply(p_adj, function(p) {
      if (is.na(p) || p >= (1 - confianca)) "" else "*"
    })

    df_posthoc <- data.frame(
      "Pares Comparados" = rotulos_pares,
      "Estat\u00edstica V" = sapply(stat_v, .fmt_num, decimais = 1),
      "p-ajustado" = sapply(p_adj, fmt_p),
      "Sig." = sig_ast,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("TESTE DE FRIEDMAN (MEDIDAS REPETIDAS / BLOCOS)")

  cat(sprintf("  Vari\u00e1vel Resposta: %s\n", y_name))
  cat(sprintf("  Fator / Condi\u00e7\u00f5es: %s (%d n\u00edveis: %s)\n", g_name, length(niveis_g), paste(niveis_g, collapse = ", ")))
  cat(sprintf("  Sujeitos / Blocos:  %s (%d blocos avaliados)\n\n", b_name, n_blocos))

  .print_topico("RESULTADO DO TESTE GLOBAL")

  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_global, align = c("center", "center", "center"))
  } else {
    print(df_global, row.names = FALSE)
  }
  cat("\n")

  if (!is.null(df_posthoc)) {
    .print_topico(sprintf("COMPARA\u00c7\u00d5ES M\u00daLTIPLAS PAREADAS (WILCOXON - AJUSTE %s)", toupper(ajuste_p)))
    if (exists(".print_tabela_estatR", mode = "function")) {
      .print_tabela_estatR(df_posthoc, align = c("center", "center", "center", "center"))
    } else {
      print(df_posthoc, row.names = FALSE)
    }
    cat("  * Diferen\u00e7a pareada estatisticamente significativa ao n\u00edvel de 5%.\n\n")
  }

  .print_topico("INTERPRETA\u00c7\u00c3O")

  cor_ok   <- "\033[1;32m"
  cor_erro <- "\033[1;31m"
  reset    <- "\033[0m"

  if (p_val < (1 - confianca)) {
    cat(sprintf("  %s[!]%s Diferen\u00e7a significativa detectada entre as condi\u00e7\u00f5es/tempos (\u03c7\u00b2 = %s; p %s).\n",
                cor_erro, reset, .fmt_num(stat_chi, 2), fmt_p(p_val, com_sinal = TRUE)))
    cat("      Rejeita-se a hip\u00f3tese de igualdade: o efeito de tratamento/tempo \u00e9 estatisticamente relevante.\n")
  } else {
    cat(sprintf("  %s[\u2713]%s N\u00e3o h\u00e1 evid\u00eancia de efeito significativo entre as condi\u00e7\u00f5es (\u03c7\u00b2 = %s; p %s).\n",
                cor_ok, reset, .fmt_num(stat_chi, 2), fmt_p(p_val, com_sinal = TRUE)))
    cat("      As distribui\u00e7\u00f5es das medidas repetidas s\u00e3o consideradas equivalentes.\n")
  }

  .print_rodape()

  # ── Gráfico Padronizado ────────────────────────────────────────────────────
  p_plot <- NULL
  if (isTRUE(grafico)) {
    df_plot <- data.frame(
      .y = as.numeric(df_clean[[y_name]]),
      .g = as.factor(df_clean[[g_name]]),
      .b = as.factor(df_clean[[b_name]])
    )

    p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .g, y = .y, group = .b)) +
      ggplot2::geom_line(color = "gray70", alpha = 0.5, linewidth = 0.5) +
      ggplot2::geom_point(color = "#377EB8", size = 1.6, alpha = 0.7) +
      ggplot2::stat_summary(ggplot2::aes(group = 1), fun = stats::median,
                            geom = "line", color = "#E31A1C", linewidth = 1.05) +
      ggplot2::stat_summary(ggplot2::aes(group = 1), fun = stats::median,
                            geom = "point", shape = 18, size = 3.5, color = "#E31A1C") +
      ggplot2::labs(
        title = paste("Evolu\u00e7\u00e3o de", y_name, "por", g_name),
        subtitle = sprintf("Friedman: \u03c7\u00b2 = %s | p = %s | Linha vermelha: Mediana", .fmt_num(stat_chi, 2), fmt_p(p_val)),
        x = g_name,
        y = y_name
      ) +
      tema_estatR(estilo = 2, grade = "dupla")

    print(p_plot)
  }

  invisible(list(
    teste_global = fr_res,
    post_hoc = df_posthoc,
    grafico = p_plot
  ))
}


# ── 3. ANÁLISE DE CONCORDÂNCIA KAPPA (COHEN & FLEISS) ────────────────────────

#' @title Análise de Concordância (Kappa de Cohen & Fleiss)
#' @description Avalia a concordância e confiabilidade entre avaliadores.
#'   Calcula automaticamente o Kappa de Cohen (simples ou ponderado linear/quadrático)
#'   para 2 avaliadores, ou o Kappa de Fleiss caso haja 3 ou mais avaliadores.
#'   Classifica a qualidade da concordância segundo os critérios de Landis & Koch.
#'
#' @param dados Data frame ou matriz de classificações.
#' @param avaliador1 Nome ou índice da coluna do primeiro avaliador (se 2 avaliadores).
#' @param avaliador2 Nome ou índice da coluna do segundo avaliador (se 2 avaliadores).
#' @param ponderado Tipo de ponderação para o Kappa de Cohen: \code{"nenhum"} (padrão),
#'   \code{"linear"} ou \code{"quadratico"} (Fleiss-Cohen, ideal para categorias ordinais).
#' @param confianca Nível de confiança para o intervalo (padrão: 0.95).
#'
#' @return Retorna invisivelmente uma lista com a matriz de concordância, coeficientes e testes.
#' @export
concordancia_kappa <- function(dados, avaliador1 = NULL, avaliador2 = NULL,
                               ponderado = c("nenhum", "linear", "quadratico"),
                               confianca = 0.95) {

  ponderado <- match.arg(ponderado)

  # Caso 1: 2 avaliadores especificados
  sub_av1 <- substitute(avaliador1)
  sub_av2 <- substitute(avaliador2)
  av1_name <- if (is.symbol(sub_av1)) as.character(sub_av1) else if (is.character(avaliador1)) avaliador1 else NULL
  av2_name <- if (is.symbol(sub_av2)) as.character(sub_av2) else if (is.character(avaliador2)) avaliador2 else NULL

  # Se fornecido dataframe com 2 colunas e sem nomes especificados
  if (is.null(av1_name) && is.null(av2_name) && is.data.frame(dados) && ncol(dados) == 2) {
    av1_name <- names(dados)[1]
    av2_name <- names(dados)[2]
  }

  # Classificação de Landis & Koch
  classificar_kappa <- function(k) {
    if (is.na(k)) return("Indeterminado")
    if (k < 0.00) return("Sem Concord\u00e2ncia (< 0,00)")
    if (k <= 0.20) return("Pobre / Leve (0,00 - 0,20)")
    if (k <= 0.40) return("Razo\u00e1vel (0,21 - 0,40)")
    if (k <= 0.60) return("Moderada (0,41 - 0,60)")
    if (k <= 0.80) return("Substancial (0,61 - 0,80)")
    return("Quase Perfeita (0,81 - 1,00)")
  }

  fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }

  # ── Fluxo A: 2 Avaliadores -> Kappa de Cohen ──────────────────────────────
  if (!is.null(av1_name) && !is.null(av2_name)) {
    v1 <- as.factor(dados[[av1_name]])
    v2 <- as.factor(dados[[av2_name]])

    # Alinha categorias em comum
    niveis_comuns <- sort(unique(c(levels(v1), levels(v2))))
    v1 <- factor(v1, levels = niveis_comuns)
    v2 <- factor(v2, levels = niveis_comuns)

    tab_mat <- table(v1, v2)
    n_tot <- sum(tab_mat)
    k_cat <- length(niveis_comuns)

    # Matriz de pesos
    W <- matrix(1, nrow = k_cat, ncol = k_cat)
    if (ponderado == "nenhum") {
      diag(W) <- 0
    } else if (ponderado == "linear") {
      for (i in 1:k_cat) {
        for (j in 1:k_cat) {
          W[i, j] <- abs(i - j) / (k_cat - 1)
        }
      }
    } else if (ponderado == "quadratico") {
      for (i in 1:k_cat) {
        for (j in 1:k_cat) {
          W[i, j] <- ((i - j)^2) / ((k_cat - 1)^2)
        }
      }
    }

    P_mat <- tab_mat / n_tot
    r_marg <- rowSums(P_mat)
    c_marg <- colSums(P_mat)
    E_mat  <- outer(r_marg, c_marg)

    po <- 1 - sum(W * P_mat)
    pe <- 1 - sum(W * E_mat)

    kappa_val <- if (pe == 1) 1 else (po - pe) / (1 - pe)

    # Variância assintótica e Erro Padrão
    se_k <- sqrt((po * (1 - po)) / (n_tot * (1 - pe)^2))
    z_score <- kappa_val / se_k
    p_val   <- 2 * stats::pnorm(-abs(z_score))

    z_crit <- stats::qnorm((1 + confianca) / 2)
    ic_inf <- max(-1, kappa_val - z_crit * se_k)
    ic_sup <- min(1, kappa_val + z_crit * se_k)

    # Formatação das tabelas
    df_coefs <- data.frame(
      "M\u00e9trica" = c("Concord\u00e2ncia Observada (Po)", "Concord\u00e2ncia Esperada (Pe)", "Kappa de Cohen"),
      "Estimativa" = c(paste0(.fmt_num(po * 100, 1), "%"), paste0(.fmt_num(pe * 100, 1), "%"), .fmt_num(kappa_val, 3)),
      "IC Inferior" = c("-", "-", .fmt_num(ic_inf, 3)),
      "IC Superior" = c("-", "-", .fmt_num(ic_sup, 3)),
      "p-valor"     = c("-", "-", fmt_p(p_val)),
      stringsAsFactors = FALSE, check.names = FALSE
    )

    # Matriz para exibição no console
    df_mat_print <- as.data.frame.matrix(tab_mat)
    df_mat_print <- cbind("Avaliador" = rownames(df_mat_print), df_mat_print)

    # Saída no Console
    .print_titulo("AN\u00c1LISE DE CONCORD\u00c2NCIA (KAPPA DE COHEN)")
    cat(sprintf("  Avaliador 1:   %s\n", av1_name))
    cat(sprintf("  Avaliador 2:   %s\n", av2_name))
    cat(sprintf("  Pondera\u00e7\u00e3o:     %s\n", switch(ponderado, nenhum = "N\u00e3o ponderado (Simples)", linear = "Linear (Ordinais)", quadratico = "Quadr\u00e1tica (Fleiss-Cohen)")))
    cat(sprintf("  Observa\u00e7\u00f5es:   %d unidades classificadas\n\n", n_tot))

    .print_topico("MATRIZ DE CONCORD\u00c2NCIA")
    if (exists(".print_tabela_estatR", mode = "function")) {
      .print_tabela_estatR(df_mat_print, align = rep("center", ncol(df_mat_print)))
    } else {
      print(tab_mat)
    }
    cat("\n")

    .print_topico("COEFICIENTES E TESTE DE SIGNIFIC\u00c2NCIA")
    if (exists(".print_tabela_estatR", mode = "function")) {
      .print_tabela_estatR(df_coefs, align = c("center", "center", "center", "center", "center"))
    } else {
      print(df_coefs, row.names = FALSE)
    }
    cat("\n")

    .print_topico("CLASSIFICA\u00c7\u00c3O DA CONCORD\u00c2NCIA (LANDIS & KOCH)")
    cor_ok   <- "\033[1;32m"
    cor_warn <- "\033[1;33m"
    reset    <- "\033[0m"

    if (kappa_val >= 0.60) {
      cat(sprintf("  %s[\u2713]%s Concord\u00e2ncia %s (Kappa = %s; p = %s).\n",
                  cor_ok, reset, toupper(classificar_kappa(kappa_val)), .fmt_num(kappa_val, 3), fmt_p(p_val)))
      cat("      A concord\u00e2ncia entre os avaliadores \u00e9 substancial/consistente al\u00e9m do acaso.\n")
    } else {
      cat(sprintf("  %s[!]%s Concord\u00e2ncia %s (Kappa = %s; p = %s).\n",
                  cor_warn, reset, toupper(classificar_kappa(kappa_val)), .fmt_num(kappa_val, 3), fmt_p(p_val)))
      cat("      Concord\u00e2ncia moderada ou fraca. H\u00e1 diverg\u00eancia consider\u00e1vel entre os crit\u00e9rios.\n")
    }

    .print_rodape()

    return(invisible(list(
      tipo = "cohen",
      matriz = tab_mat,
      kappa = kappa_val,
      ic = c(ic_inf, ic_sup),
      p_valor = p_val
    )))
  }

  # ── Fluxo B: >= 3 Avaliadores -> Kappa de Fleiss ──────────────────────────
  # Matriz de sujeitos x categorias onde cada célula é o número de avaliadores que escolheram a categoria
  mat_fleiss <- as.matrix(dados)
  N_suj <- nrow(mat_fleiss)
  k_cats <- ncol(mat_fleiss)
  n_raters <- sum(mat_fleiss[1, ])

  p_j <- colSums(mat_fleiss) / (N_suj * n_raters)
  pe_fleiss <- sum(p_j^2)

  P_i <- (rowSums(mat_fleiss^2) - n_raters) / (n_raters * (n_raters - 1))
  po_fleiss <- mean(P_i)

  kappa_fleiss <- (po_fleiss - pe_fleiss) / (1 - pe_fleiss)

  var_fleiss <- (2 / (N_suj * n_raters * (n_raters - 1) * (1 - pe_fleiss)^2)) *
    (pe_fleiss^2 - sum(p_j^2 * (1 - p_j + sum(p_j^2))))
  se_fleiss <- sqrt(max(0, var_fleiss))
  z_score <- if (se_fleiss > 0) kappa_fleiss / se_fleiss else 0
  p_val <- 2 * stats::pnorm(-abs(z_score))

  z_crit <- stats::qnorm((1 + confianca) / 2)
  ic_inf <- max(-1, kappa_fleiss - z_crit * se_fleiss)
  ic_sup <- min(1, kappa_fleiss + z_crit * se_fleiss)

  df_coefs_fl <- data.frame(
    "M\u00e9trica" = c("Concord\u00e2ncia Observada (Po)", "Concord\u00e2ncia Esperada (Pe)", "Kappa de Fleiss"),
    "Estimativa" = c(paste0(.fmt_num(po_fleiss * 100, 1), "%"), paste0(.fmt_num(pe_fleiss * 100, 1), "%"), .fmt_num(kappa_fleiss, 3)),
    "IC Inferior" = c("-", "-", .fmt_num(ic_inf, 3)),
    "IC Superior" = c("-", "-", .fmt_num(ic_sup, 3)),
    "p-valor"     = c("-", "-", fmt_p(p_val)),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  .print_titulo("AN\u00c1LISE DE CONCORD\u00c2NCIA (KAPPA DE FLEISS)")
  cat(sprintf("  Avaliadores:   %d avaliadores por unidade\n", n_raters))
  cat(sprintf("  Unidades:      %d sujeitos/itens avaliados\n", N_suj))
  cat(sprintf("  Categorias:    %d op\u00e7\u00f5es dispon\u00edveis\n\n", k_cats))

  .print_topico("COEFICIENTES E TESTE DE SIGNIFIC\u00c2NCIA")
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_coefs_fl, align = c("center", "center", "center", "center", "center"))
  } else {
    print(df_coefs_fl, row.names = FALSE)
  }
  cat("\n")

  .print_topico("CLASSIFICA\u00c7\u00c3O DA CONCORD\u00c2NCIA (LANDIS & KOCH)")
  cor_ok   <- "\033[1;32m"
  cor_warn <- "\033[1;33m"
  reset    <- "\033[0m"

  if (kappa_fleiss >= 0.60) {
    cat(sprintf("  %s[\u2713]%s Concord\u00e2ncia %s (Kappa = %s; p = %s).\n",
                cor_ok, reset, toupper(classificar_kappa(kappa_fleiss)), .fmt_num(kappa_fleiss, 3), fmt_p(p_val)))
  } else {
    cat(sprintf("  %s[!]%s Concord\u00e2ncia %s (Kappa = %s; p = %s).\n",
                cor_warn, reset, toupper(classificar_kappa(kappa_fleiss)), .fmt_num(kappa_fleiss, 3), fmt_p(p_val)))
  }

  .print_rodape()

  invisible(list(
    tipo = "fleiss",
    kappa = kappa_fleiss,
    ic = c(ic_inf, ic_sup),
    p_valor = p_val
  ))
}


# ── 4. CONCORDÂNCIA CONTÍNUA DE BLAND-ALTMAN ─────────────────────────────────

#' @title Gráfico e Diagnóstico de Bland-Altman
#' @description Avalia a concordância entre dois métodos contínuos de medição
#'   (ex: dois instrumentos, ensaio laboratorial vs teste rápido).
#'   Estima o viés médio (\eqn{\bar{d}}{\code{d_barra}}), limites de concordância de 95\% (\eqn{\bar{d} \pm 1{,}96 \cdot DP}{\code{d_barra +- 1.96*DP}}),
#'   intervalos de confiança dos parâmetros e gera gráfico no padrão visual do estatR.
#'
#' @param dados Data frame contendo as medições.
#' @param metodo1 Coluna numérica com as medidas do primeiro método.
#' @param metodo2 Coluna numérica com as medidas do segundo método.
#' @param confianca Nível de confiança para os intervalos (padrão: 0.95).
#' @param grafico Lógico. Se \code{TRUE} (padrão), plota o gráfico de Bland-Altman.
#'
#' @return Retorna invisivelmente uma lista com parâmetros calculados e o gráfico ggplot.
#' @export
grafico_bland_altman <- function(dados, metodo1, metodo2, confianca = 0.95, grafico = TRUE) {

  m1_sub <- substitute(metodo1)
  m2_sub <- substitute(metodo2)
  m1_name <- if (is.symbol(m1_sub)) as.character(m1_sub) else if (is.character(metodo1)) metodo1 else NULL
  m2_name <- if (is.symbol(m2_sub)) as.character(m2_sub) else if (is.character(metodo2)) metodo2 else NULL

  if (is.null(m1_name) || is.null(m2_name)) {
    stop("Informe as duas vari\u00e1veis num\u00e9ricas: 'metodo1' e 'metodo2'.")
  }

  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  if (!m1_name %in% names(dados)) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", m1_name))
  if (!m2_name %in% names(dados)) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", m2_name))

  v1 <- as.numeric(dados[[m1_name]])
  v2 <- as.numeric(dados[[m2_name]])

  valid <- !is.na(v1) & !is.na(v2)
  v1 <- v1[valid]
  v2 <- v2[valid]
  n_obs <- length(v1)

  if (n_obs < 4) stop("Quantidade insuficiente de pares v\u00e1lidos (m\u00ednimo 4).")

  diferenca <- v1 - v2
  media_par <- (v1 + v2) / 2

  vies_medio <- mean(diferenca)
  dp_dif     <- stats::sd(diferenca)

  # Limites de concordância de 95%
  k_sd <- 1.96
  lim_inf <- vies_medio - k_sd * dp_dif
  lim_sup <- vies_medio + k_sd * dp_dif

  # Erros padrão assintóticos segundo Bland & Altman (1986)
  se_vies <- sqrt(dp_dif^2 / n_obs)
  se_lim  <- sqrt((3 * dp_dif^2) / n_obs)

  t_crit <- stats::qt((1 + confianca) / 2, df = n_obs - 1)

  # Intervalos de confiança de 95%
  ic_vies_inf <- vies_medio - t_crit * se_vies
  ic_vies_sup <- vies_medio + t_crit * se_vies

  ic_lim_inf_baixo <- lim_inf - t_crit * se_lim
  ic_lim_inf_alto  <- lim_inf + t_crit * se_lim

  ic_lim_sup_baixo <- lim_sup - t_crit * se_lim
  ic_lim_sup_alto  <- lim_sup + t_crit * se_lim

  # Teste t para significância do viés (H0: vies = 0)
  t_vies <- stats::t.test(diferenca, mu = 0)
  p_vies <- t_vies$p.value

  # Pontos fora dos limites
  fora_limites <- sum(diferenca < lim_inf | diferenca > lim_sup)
  pct_fora <- (fora_limites / n_obs) * 100

  fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }

  df_params <- data.frame(
    "Par\u00e2metro" = c("Vi\u00e9s M\u00e9dio (d)", "Desvio Padr\u00e3o das Dif. (s)",
                       "Limite Inferior (d - 1,96s)", "Limite Superior (d + 1,96s)"),
    "Valor" = c(.fmt_num(vies_medio, 2), .fmt_num(dp_dif, 2), .fmt_num(lim_inf, 2), .fmt_num(lim_sup, 2)),
    "IC Inferior" = c(.fmt_num(ic_vies_inf, 2), "-", .fmt_num(ic_lim_inf_baixo, 2), .fmt_num(ic_lim_sup_baixo, 2)),
    "IC Superior" = c(.fmt_num(ic_vies_sup, 2), "-", .fmt_num(ic_lim_inf_alto, 2), .fmt_num(ic_lim_sup_alto, 2)),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("AN\u00c1LISE DE CONCORD\u00c2NCIA DE BLAND-ALTMAN")

  cat(sprintf("  M\u00e9todo 1: %s\n", m1_name))
  cat(sprintf("  M\u00e9todo 2: %s\n", m2_name))
  cat(sprintf("  Pares avaliados: %d observa\u00e7\u00f5es\n\n", n_obs))

  .print_topico("PAR\u00c2METROS DE CONCORD\u00c2NCIA")

  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_params, align = c("center", "center", "center", "center"))
  } else {
    print(df_params, row.names = FALSE)
  }
  cat("\n")

  .print_topico("DIAGN\u00d3STICO E CONCORD\u00c2NCIA")
  cat(sprintf("  * Pontos fora dos limites de 95%%: %d de %d (%s - esperado \u2264 5,0%%).\n\n",
              fora_limites, n_obs, paste0(.fmt_num(pct_fora, 1), "%")))

  .print_topico("INTERPRETA\u00c7\u00c3O")

  cor_ok   <- "\033[1;32m"
  cor_erro <- "\033[1;31m"
  reset    <- "\033[0m"

  if (pct_fora <= 5.0) {
    cat(sprintf("  %s[\u2713]%s Boa concord\u00e2ncia global: %.1f%% dos pontos encontram-se dentro dos limites.\n",
                cor_ok, reset, 100 - pct_fora))
  } else {
    cat(sprintf("  %s[!]%s Discrep\u00e2ncia moderada: %.1f%% dos pontos excedem os limites de 95%% (ideal \u2264 5%%).\n",
                cor_erro, reset, pct_fora))
  }

  if (p_vies < 0.05) {
    direcao_vies <- if (vies_medio > 0) paste(m1_name, "superestima", m2_name) else paste(m1_name, "subestima", m2_name)
    cat(sprintf("  %s[!]%s Vi\u00e9s sistem\u00e1tico estatisticamente significativo (%s; p = %s).\n",
                cor_erro, reset, .fmt_num(vies_medio, 2), fmt_p(p_vies)))
    cat(sprintf("      Em m\u00e9dia, %s por %s unidades.\n", direcao_vies, .fmt_num(abs(vies_medio), 2)))
  } else {
    cat(sprintf("  %s[\u2713]%s Aus\u00eancia de vi\u00e9s sistem\u00e1tico relevante (d = %s; p = %s).\n",
                cor_ok, reset, .fmt_num(vies_medio, 2), fmt_p(p_vies)))
  }

  .print_rodape()

  # ── Gráfico Padronizado ────────────────────────────────────────────────────
  p_plot <- NULL
  if (isTRUE(grafico)) {
    df_plot <- data.frame(
      Media = media_par,
      Diferenca = diferenca
    )

    p_plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Media, y = Diferenca)) +
      # Faixas de confiança dos limites
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = ic_lim_sup_baixo, ymax = ic_lim_sup_alto,
                        fill = "#E31A1C", alpha = 0.12) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = ic_lim_inf_baixo, ymax = ic_lim_inf_alto,
                        fill = "#E31A1C", alpha = 0.12) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = ic_vies_inf, ymax = ic_vies_sup,
                        fill = "#377EB8", alpha = 0.15) +
      # Linhas de referência
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4, linetype = "solid") +
      ggplot2::geom_hline(yintercept = vies_medio, color = "#377EB8", linewidth = 0.9) +
      ggplot2::geom_hline(yintercept = c(lim_inf, lim_sup), color = "#E31A1C", linetype = "dashed", linewidth = 0.75) +
      # Pontos
      ggplot2::geom_point(color = "gray25", size = 2, alpha = 0.75) +
      # Rótulos nas linhas
      ggplot2::annotate("text", x = max(media_par), y = lim_sup, label = paste0("+1,96 DP: ", .fmt_num(lim_sup, 2)),
                        hjust = 1, vjust = -0.5, color = "#E31A1C", fontface = "bold", size = 3.3) +
      ggplot2::annotate("text", x = max(media_par), y = vies_medio, label = paste0("Vi\u00e9s: ", .fmt_num(vies_medio, 2)),
                        hjust = 1, vjust = -0.5, color = "#377EB8", fontface = "bold", size = 3.3) +
      ggplot2::annotate("text", x = max(media_par), y = lim_inf, label = paste0("-1,96 DP: ", .fmt_num(lim_inf, 2)),
                        hjust = 1, vjust = 1.3, color = "#E31A1C", fontface = "bold", size = 3.3) +
      ggplot2::labs(
        title = paste("Gr\u00e1fico de Bland-Altman:", m1_name, "vs", m2_name),
        subtitle = sprintf("Vi\u00e9s M\u00e9dio: %s | Limites 95%%: [%s ; %s]", .fmt_num(vies_medio, 2), .fmt_num(lim_inf, 2), .fmt_num(lim_sup, 2)),
        x = paste("M\u00e9dia entre os M\u00e9todos (", m1_name, "+", m2_name, ") / 2"),
        y = paste("Diferen\u00e7a entre os M\u00e9todos (", m1_name, "-", m2_name, ")")
      ) +
      tema_estatR(estilo = 2, grade = "dupla")

    print(p_plot)
  }

  invisible(list(
    vies = vies_medio,
    dp = dp_dif,
    limites = c(inferior = lim_inf, superior = lim_sup),
    ic_vies = c(inferior = ic_vies_inf, superior = ic_vies_sup),
    pontos_fora = c(total = fora_limites, pct = pct_fora),
    grafico = p_plot
  ))
}
