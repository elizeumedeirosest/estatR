# ==============================================================================
# MÓDULO: SÉRIES TEMPORAIS
# Análise, Decomposição, Estacionariedade, Autocorrelação, Previsão e Seleção
# ==============================================================================

# ── 1. HELPER INTERNO: EXTRAÇÃO E VALIDAÇÃO DA SÉRIE ─────────────────────────
#' @noRd
.extrair_serie_temporal <- function(dados, y = NULL, data = NULL, frequencia = NULL) {
  # Caso 1: Objeto 'ts' direto
  if (stats::is.ts(dados)) {
    freq <- if (!is.null(frequencia)) frequencia else stats::frequency(dados)
    return(stats::ts(as.numeric(dados), start = stats::start(dados), frequency = freq))
  }

  # Caso 2: Vetor numérico direto
  if (is.numeric(dados) && is.null(y)) {
    freq <- if (!is.null(frequencia) && frequencia >= 1) frequencia else 12
    return(stats::ts(as.numeric(dados), frequency = freq))
  }

  # Caso 3: Data frame fornecido
  if (!is.data.frame(dados)) {
    stop("O argumento 'dados' deve ser um data frame, vetor num\u00e9rico ou objeto 'ts'.")
  }

  # Avalia nome de y
  y_col <- tryCatch({
    y_sub <- substitute(y)
    if (is.symbol(y_sub) || is.character(y_sub)) as.character(y_sub) else NULL
  }, error = function(e) NULL)

  if (is.null(y_col) || !y_col %in% names(dados)) {
    # Tenta como string direta se y foi passado como caractere
    if (is.character(y) && length(y) == 1 && y %in% names(dados)) {
      y_col <- y
    } else {
      # Procura a primeira coluna numérica
      nums <- names(dados)[sapply(dados, is.numeric)]
      if (length(nums) == 0) stop("Nenhuma vari\u00e1vel num\u00e9rica encontrada no conjunto de dados.")
      y_col <- nums[1]
    }
  }

  df_ord <- dados

  # Ordenação por coluna de data, se especificada
  data_col <- tryCatch({
    d_sub <- substitute(data)
    if (is.symbol(d_sub) || is.character(d_sub)) as.character(d_sub) else NULL
  }, error = function(e) NULL)

  if (!is.null(data_col) && data_col %in% names(df_ord)) {
    df_ord <- df_ord[order(df_ord[[data_col]]), , drop = FALSE]
  }

  vals <- as.numeric(df_ord[[y_col]])
  vals <- vals[!is.na(vals)]

  if (length(vals) < 4) {
    stop("A s\u00e9rie temporal deve conter pelo menos 4 observa\u00e7\u00f5es v\u00e1lidas.")
  }

  freq <- if (!is.null(frequencia) && frequencia >= 1) frequencia else 12
  stats::ts(vals, frequency = freq)
}


# ── 2. DIAGNÓSTICO E DECOMPOSIÇÃO DE SÉRIE TEMPORAL ─────────────────────────

#' @title Diagnóstico e Decomposição de Série Temporal
#' @description Realiza a análise estatística descritiva e a decomposição temporal
#'   (tendência, sazonalidade e componente aleatório/ruído) de uma série.
#'   Gera um painel gráfico com 4 gráficos alinhados com o padrão visual do estatR.
#'
#' @param dados Data frame contendo os dados, vetor numérico ou objeto de classe \code{ts}.
#' @param y Variável numérica resposta (coluna da série). Aceita símbolo ou string.
#' @param data Coluna de data opcional para ordenação cronológica.
#' @param frequencia Periodicidade da série (ex: 12 para mensal, 4 para trimestral, 7 para semanal). Padrão: 12.
#' @param tipo Tipo de decomposição: \code{"aditiva"} (padrão) ou \code{"multiplicativa"}.
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o painel visual com os 4 componentes.
#'
#' @return Objeto invisível da classe \code{estatR_analise_serie} contendo as estatísticas descritivas,
#'   a decomposição e o gráfico gerado.
#' @export
analise_serie <- function(dados, y = NULL, data = NULL, frequencia = 12,
                          tipo = c("aditiva", "multiplicativa"), grafico = TRUE) {

  tipo <- match.arg(tipo)
  y_expr <- substitute(y)
  data_expr <- substitute(data)

  y_str <- if (is.symbol(y_expr)) as.character(y_expr) else if (is.character(y)) y else NULL
  data_str <- if (is.symbol(data_expr)) as.character(data_expr) else if (is.character(data)) data else NULL

  # Extrai série ts
  ts_serie <- if (is.data.frame(dados)) {
    .extrair_serie_temporal(dados, y = y_str, data = data_str, frequencia = frequencia)
  } else {
    .extrair_serie_temporal(dados, frequencia = frequencia)
  }

  n_obs <- length(ts_serie)
  freq_efetiva <- stats::frequency(ts_serie)

  if (tipo == "multiplicativa" && any(ts_serie <= 0)) {
    stop("A decomposi\u00e7\u00e3o multiplicativa exige valores estritamente positivos (y > 0). Use tipo = 'aditiva'.")
  }

  # Estatísticas descritivas
  media_val <- mean(ts_serie)
  mediana_val <- stats::median(ts_serie)
  dp_val <- stats::sd(ts_serie)
  cv_val <- (dp_val / media_val) * 100
  min_val <- min(ts_serie)
  max_val <- max(ts_serie)
  assimetria_val <- if (exists(".calcular_assimetria", mode = "function")) .calcular_assimetria(ts_serie) else NA
  curtose_val <- if (exists(".calcular_curtose", mode = "function")) .calcular_curtose(ts_serie) else NA

  # Decomposição
  tipo_stats <- if (tipo == "multiplicativa") "multiplicative" else "additive"
  decomposicao <- tryCatch({
    if (n_obs >= 2 * freq_efetiva) {
      stats::decompose(ts_serie, type = tipo_stats)
    } else {
      # Fallback se a série tiver menos de 2 ciclos completos
      NULL
    }
  }, error = function(e) NULL)

  forca_tendencia <- NA_real_
  forca_sazonal <- NA_real_

  if (!is.null(decomposicao)) {
    tt <- as.numeric(stats::na.omit(decomposicao$trend))
    ss <- as.numeric(stats::na.omit(decomposicao$seasonal))
    rr <- as.numeric(stats::na.omit(decomposicao$random))

    # Variâncias para forças de Granger & Hyndman
    var_r <- stats::var(rr, na.rm = TRUE)
    var_tr <- stats::var(decomposicao$trend + decomposicao$random, na.rm = TRUE)
    var_sr <- stats::var(decomposicao$seasonal + decomposicao$random, na.rm = TRUE)

    forca_tendencia <- if (!is.na(var_tr) && var_tr > 0) max(0, 1 - (var_r / var_tr)) else 0
    forca_sazonal   <- if (!is.na(var_sr) && var_sr > 0) max(0, 1 - (var_r / var_sr)) else 0
  }

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("DIAGN\u00d3STICO E DECOMPOSI\u00c7\u00c3O DE S\u00c9RIE TEMPORAL")

  inicio_txt <- paste(stats::start(ts_serie), collapse = ".")
  fim_txt    <- paste(stats::end(ts_serie), collapse = ".")

  cat(sprintf("  Observa\u00e7\u00f5es:  %d\n", n_obs))
  cat(sprintf("  Frequ\u00eancia:   %d (%s)\n", freq_efetiva,
              if (freq_efetiva == 12) "Mensal" else if (freq_efetiva == 4) "Trimestral" else if (freq_efetiva == 7) "Semanal" else "Customizada"))
  cat(sprintf("  Per\u00edodo:      %s a %s\n", inicio_txt, fim_txt))
  cat(sprintf("  Modelo:       Decomposi\u00e7\u00e3o %s\n\n", if (tipo == "aditiva") "Aditiva (Y = T + S + R)" else "Multiplicativa (Y = T \u00d7 S \u00d7 R)"))

  .print_topico("ESTAT\u00cdSTICAS DESCRITIVAS")

  df_desc <- data.frame(
    "M\u00e9trica" = c("M\u00e9dia", "Mediana", "Desvio Padr\u00e3o", "CV (%)", "M\u00ednimo", "M\u00e1ximo", "Assimetria", "Curtose"),
    Valor   = c(
      .fmt_num(media_val, 2),
      .fmt_num(mediana_val, 2),
      .fmt_num(dp_val, 2),
      paste0(.fmt_num(cv_val, 1), "%"),
      .fmt_num(min_val, 2),
      .fmt_num(max_val, 2),
      if (is.na(assimetria_val)) "-" else .fmt_num(assimetria_val, 2),
      if (is.na(curtose_val)) "-" else .fmt_num(curtose_val, 2)
    ),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_desc, align = c("left", "center"))
  } else {
    print(df_desc, row.names = FALSE)
  }
  cat("\n")

  if (!is.null(decomposicao)) {
    .print_topico("FOR\u00c7A DOS COMPONENTES (HYNDMAN)")
    df_decomp_tab <- data.frame(
      Componente = c("For\u00e7a da Tend\u00eancia (Ft)", "For\u00e7a da Sazonalidade (Fs)"),
      Valor      = c(.fmt_num(forca_tendencia, 3), .fmt_num(forca_sazonal, 3)),
      "Avalia\u00e7\u00e3o"  = c(
        if (forca_tendencia > 0.6) "Forte" else if (forca_tendencia > 0.3) "Moderada" else "Fraca",
        if (forca_sazonal > 0.6) "Forte" else if (forca_sazonal > 0.3) "Moderada" else "Fraca"
      ),
      stringsAsFactors = FALSE, check.names = FALSE
    )

    if (exists(".print_tabela_estatR", mode = "function")) {
      .print_tabela_estatR(df_decomp_tab, align = c("left", "center", "center"))
    } else {
      print(df_decomp_tab, row.names = FALSE)
    }
  } else {
    cat("  * Observa\u00e7\u00e3o: S\u00e9rie com menos de 2 ciclos completos para c\u00e1lculo formal da decomposi\u00e7\u00e3o.\n")
  }

  .print_rodape()

  # ── Painel Gráfico ─────────────────────────────────────────────────────────
  p_final <- NULL
  if (isTRUE(grafico)) {
    tempos <- as.numeric(stats::time(ts_serie))

    if (!is.null(decomposicao)) {
      df_plot <- data.frame(
        Tempo       = rep(tempos, 4),
        Valor       = c(as.numeric(ts_serie), as.numeric(decomposicao$trend), as.numeric(decomposicao$seasonal), as.numeric(decomposicao$random)),
        Componente  = factor(rep(c("1. S\u00e9rie Original", "2. Tend\u00eancia", "3. Sazonalidade", "4. Ru\u00eddo (Aleat\u00f3rio)"), each = n_obs),
                             levels = c("1. S\u00e9rie Original", "2. Tend\u00eancia", "3. Sazonalidade", "4. Ru\u00eddo (Aleat\u00f3rio)"))
      )

      p1 <- ggplot2::ggplot(data.frame(Tempo = tempos, Valor = as.numeric(ts_serie)), ggplot2::aes(x = Tempo, y = Valor)) +
        ggplot2::geom_line(color = "#377EB8", linewidth = 0.85) +
        ggplot2::labs(title = "S\u00e9rie Original", x = NULL, y = "Valor") +
        tema_estatR(estilo = 2, grade = "dupla")

      p2 <- ggplot2::ggplot(data.frame(Tempo = tempos, Valor = as.numeric(decomposicao$trend)), ggplot2::aes(x = Tempo, y = Valor)) +
        ggplot2::geom_line(color = "#E31A1C", linewidth = 0.85, na.rm = TRUE) +
        ggplot2::labs(title = "Tend\u00eancia", x = NULL, y = "Tend\u00eancia") +
        tema_estatR(estilo = 2, grade = "dupla")

      p3 <- ggplot2::ggplot(data.frame(Tempo = tempos, Valor = as.numeric(decomposicao$seasonal)), ggplot2::aes(x = Tempo, y = Valor)) +
        ggplot2::geom_line(color = "seagreen4", linewidth = 0.85) +
        ggplot2::labs(title = "Sazonalidade", x = NULL, y = "Sazonal") +
        tema_estatR(estilo = 2, grade = "dupla")

      p4 <- ggplot2::ggplot(data.frame(Tempo = tempos, Valor = as.numeric(decomposicao$random)), ggplot2::aes(x = Tempo, y = Valor)) +
        ggplot2::geom_hline(yintercept = if (tipo == "aditiva") 0 else 1, linetype = "dashed", color = "gray50") +
        ggplot2::geom_line(color = "gray30", linewidth = 0.75, na.rm = TRUE) +
        ggplot2::labs(title = "Ru\u00eddo / Aleat\u00f3rio", x = "Tempo", y = "Res\u00edduo") +
        tema_estatR(estilo = 2, grade = "dupla")

      p_final <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 1)
    } else {
      p_final <- ggplot2::ggplot(data.frame(Tempo = tempos, Valor = as.numeric(ts_serie)), ggplot2::aes(x = Tempo, y = Valor)) +
        ggplot2::geom_line(color = "#377EB8", linewidth = 0.9) +
        ggplot2::labs(title = "S\u00e9rie Temporal", x = "Tempo", y = "Valor") +
        tema_estatR(estilo = 2, grade = "dupla")
    }

    print(p_final)
  }

  invisible(structure(
    list(
      estatisticas = df_desc,
      decomposicao = decomposicao,
      forca_tendencia = forca_tendencia,
      forca_sazonal = forca_sazonal,
      grafico = p_final,
      serie = ts_serie
    ),
    class = "estatR_analise_serie"
  ))
}


# ── 3. TESTES DE ESTACIONARIEDADE ───────────────────────────────────────────

#' @title Testes de Estacionariedade (ADF e KPSS)
#' @description Aplica conjuntamente os testes de Dickey-Fuller Aumentado (ADF)
#'   e KPSS para avaliar a estacionariedade da série e a presença de raiz unitária.
#'   Apresenta diagnóstico com veredito consensual padronizado com \code{[✓]} e \code{[!]}.
#'
#' @param serie Vetor numérico, objeto \code{ts} ou data frame.
#' @param frequencia Periodicidade da série (opcional).
#' @param diferencas Ordem de diferenciação a ser aplicada antes do teste (padrão: 0, avalia a série em nível).
#'
#' @return Retorna invisivelmente um data frame com os resultados estatísticos de cada teste.
#' @export
teste_estacionariedade <- function(serie, frequencia = NULL, diferencas = 0) {
  ts_serie <- .extrair_serie_temporal(serie, frequencia = frequencia)

  if (diferencas > 0) {
    ts_serie <- base::diff(ts_serie, differences = as.integer(diferencas))
  }

  vals <- as.numeric(ts_serie)
  vals <- vals[!is.na(vals)]
  n_obs <- length(vals)

  # Teste ADF (H0: Raiz unitária / Não estacionária)
  res_adf <- tryCatch({
    suppressWarnings(tseries::adf.test(vals))
  }, error = function(e) NULL)

  # Teste KPSS (H0: Estacionária no nível)
  res_kpss <- tryCatch({
    suppressWarnings(tseries::kpss.test(vals, null = "Level"))
  }, error = function(e) NULL)

  stat_adf <- if (!is.null(res_adf)) res_adf$statistic else NA_real_
  p_adf    <- if (!is.null(res_adf)) res_adf$p.value else NA_real_

  stat_kpss <- if (!is.null(res_kpss)) res_kpss$statistic else NA_real_
  p_kpss    <- if (!is.null(res_kpss)) res_kpss$p.value else NA_real_

  # Diagnóstico individual
  # ADF: rejeita H0 se p < 0.05 -> Estacionária
  adf_estac <- !is.na(p_adf) && (p_adf < 0.05)
  # KPSS: não rejeita H0 se p >= 0.05 -> Estacionária
  kpss_estac <- !is.na(p_kpss) && (p_kpss >= 0.05)

  # Formatação dos p-valores
  fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p <= 0.01) return("\u2264 0,010")
    if (p >= 0.10) return("\u2265 0,100")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }

  df_tab <- data.frame(
    Teste = c("ADF (Dickey-Fuller)", "KPSS (N\u00edvel)"),
    "Hip\u00f3tese Nula (H\u2080)" = c("Possui raiz unit\u00e1ria (n\u00e3o estacion\u00e1ria)", "Estacion\u00e1ria em torno do n\u00edvel"),
    "Estat\u00edstica" = c(
      if (is.na(stat_adf)) "-" else .fmt_num(stat_adf, 3),
      if (is.na(stat_kpss)) "-" else .fmt_num(stat_kpss, 3)
    ),
    "p-valor" = c(fmt_p(p_adf), fmt_p(p_kpss)),
    "Diagn\u00f3stico" = c(
      if (adf_estac) "Estacion\u00e1ria" else "N\u00e3o Estacion\u00e1ria",
      if (kpss_estac) "Estacion\u00e1ria" else "N\u00e3o Estacion\u00e1ria"
    ),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("TESTES DE ESTACIONARIEDADE (ADF & KPSS)")

  cat(sprintf("  S\u00e9rie avaliada:  %d observa\u00e7\u00f5es\n", n_obs))
  cat(sprintf("  Diferencia\u00e7\u00e3o:    d = %d %s\n\n",
              diferencas,
              if (diferencas == 0) "(S\u00e9rie em n\u00edvel original)" else sprintf("(Ap\u00f3s %d diferen\u00e7a(s))", diferencas)))

  .print_topico("RESULTADOS DOS TESTES")

  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_tab, align = c("left", "left", "center", "center", "center"))
  } else {
    print(df_tab, row.names = FALSE)
  }
  cat("\n")

  .print_topico("INTERPRETA\u00c7\u00c3O")

  cor_ok   <- "\033[1;32m"
  cor_erro <- "\033[1;31m"
  cor_warn <- "\033[1;33m"
  reset    <- "\033[0m"

  if (adf_estac && kpss_estac) {
    cat(sprintf("  %s[\u2713]%s A s\u00e9rie \u00e9 ESTACION\u00c1RIA ao n\u00edvel de signific\u00e2ncia de 5%%.\n", cor_ok, reset))
    cat("      Ambos os testes concordam: n\u00e3o h\u00e1 evid\u00eancia de raiz unit\u00e1ria e a vari\u00e2ncia/m\u00e9dia s\u00e3o est\u00e1veis.\n")
    cat("      Pronta para modelagem direta (ex: ARMA/ARIMA com d = 0).\n")
  } else if (!adf_estac && !kpss_estac) {
    cat(sprintf("  %s[!]%s A s\u00e9rie N\u00c3O \u00c9 ESTACION\u00c1RIA ao n\u00edvel de signific\u00e2ncia de 5%%.\n", cor_erro, reset))
    cat("      Ambos os testes detectam presen\u00e7a de raiz unit\u00e1ria ou instabilidade no n\u00edvel.\n")
    cat("      Recomenda-se aplicar diferencia\u00e7\u00e3o (ex: diferencas = 1) ou transforma\u00e7\u00e3o Box-Cox/logar\u00edtmica.\n")
  } else if (adf_estac && !kpss_estac) {
    cat(sprintf("  %s[!]%s Resultados DIVERGENTES entre os testes.\n", cor_warn, reset))
    cat("      O teste ADF rejeitou a raiz unit\u00e1ria, mas o KPSS rejeitou a estacionariedade no n\u00edvel.\n")
    cat("      Indica prov\u00e1vel estacionariedade em torno de uma tend\u00eancia determin\u00edstica (Trend-Stationary).\n")
  } else {
    cat(sprintf("  %s[!]%s Resultados DIVERGENTES entre os testes.\n", cor_warn, reset))
    cat("      O teste KPSS n\u00e3o rejeitou a estacionariedade, mas o ADF n\u00e3o rejeitou a presen\u00e7a de raiz unit\u00e1ria.\n")
    cat("      S\u00e9rie com comportamento de persist\u00eancia moderada. Avalie o gr\u00e1fico de autocorrela\u00e7\u00e3o (FAC).\n")
  }

  .print_rodape()

  invisible(df_tab)
}


# ── 4. AUTOCORRELAÇÃO (FAC E FACP COM LJUNG-BOX) ────────────────────────────

#' @title Gráfico de Autocorrelação (FAC e FACP)
#' @description Plota as funções de autocorrelação simples (FAC) e parcial (FACP)
#'   com limites de significância de Bartlett (95%) no estilo pirulito (lollipop)
#'   e executa opcionalmente o teste de Ljung-Box para diagnóstico de ruído branco.
#'
#' @param serie Vetor numérico, objeto \code{ts} ou data frame.
#' @param defasagens Número máximo de defasagens (lags) a exibir (padrão: 24).
#' @param posicao Disposição dos painéis: \code{"empilhado"} (vertical) ou \code{"lado_a_lado"}.
#' @param teste_ruido Lógico. Se \code{TRUE} (padrão), exibe o teste de Ljung-Box no console.
#'
#' @return Objeto ggplot combinado via patchwork.
#' @export
grafico_autocorrelacao <- function(serie, defasagens = 24,
                                   posicao = c("empilhado", "lado_a_lado"),
                                   teste_ruido = TRUE) {

  posicao <- match.arg(posicao)
  ts_serie <- .extrair_serie_temporal(serie)
  vals <- as.numeric(ts_serie)
  n_obs <- length(vals)

  max_lag <- min(defasagens, n_obs - 2)

  # Cálculo FAC e FACP
  acf_res  <- stats::acf(vals, lag.max = max_lag, plot = FALSE)
  pacf_res <- stats::pacf(vals, lag.max = max_lag, plot = FALSE)

  ic_val <- 1.96 / sqrt(n_obs)

  # Data frames
  df_acf <- data.frame(
    Lag   = as.numeric(acf_res$lag)[-1],
    Valor = as.numeric(acf_res$acf)[-1],
    Tipo  = "FAC (Autocorrela\u00e7\u00e3o Simples)"
  )

  df_pacf <- data.frame(
    Lag   = as.numeric(pacf_res$lag),
    Valor = as.numeric(pacf_res$acf),
    Tipo  = "FACP (Autocorrela\u00e7\u00e3o Parcial)"
  )

  # Helper do gráfico
  .plot_acf_panel <- function(df, titulo) {
    ggplot2::ggplot(df, ggplot2::aes(x = Lag, y = Valor)) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = c(-ic_val, ic_val), linetype = "dashed", color = "#E31A1C", linewidth = 0.6) +
      ggplot2::geom_segment(ggplot2::aes(x = Lag, xend = Lag, y = 0, yend = Valor),
                           color = "#377EB8", linewidth = 0.8) +
      ggplot2::geom_point(color = "#377EB8", size = 2) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      ggplot2::scale_y_continuous(limits = c(min(-ic_val * 1.2, min(df$Valor, na.rm = TRUE) * 1.1),
                                             max(ic_val * 1.2, max(df$Valor, na.rm = TRUE) * 1.1))) +
      ggplot2::labs(title = titulo, x = "Defasagem (Lag)", y = "Autocorrela\u00e7\u00e3o") +
      tema_estatR(estilo = 2, grade = "dupla")
  }

  p_acf  <- .plot_acf_panel(df_acf, "Fun\u00e7\u00e3o de Autocorrela\u00e7\u00e3o (FAC)")
  p_pacf <- .plot_acf_panel(df_pacf, "Fun\u00e7\u00e3o de Autocorrela\u00e7\u00e3o Parcial (FACP)")

  p_final <- if (posicao == "empilhado") {
    patchwork::wrap_plots(p_acf, p_pacf, ncol = 1)
  } else {
    patchwork::wrap_plots(p_acf, p_pacf, ncol = 2)
  }

  # Teste de Ljung-Box no console
  if (isTRUE(teste_ruido)) {
    box_res <- stats::Box.test(vals, lag = min(12, max_lag), type = "Ljung-Box")

    .print_titulo("TESTE DE RU\u00cdDO BRANCO (LJUNG-BOX)")
    cat(sprintf("  Hip\u00f3tese Nula (H\u2080): Os res\u00edduos s\u00e3o ru\u00eddo branco (aus\u00eancia de autocorrela\u00e7\u00e3o).\n\n"))

    p_box <- box_res$p.value
    fmt_pbox <- if (p_box < 0.001) "< 0,001" else formatC(p_box, format = "f", digits = 3, decimal.mark = ",")

    df_box <- data.frame(
      "Defasagens (Lags)" = as.character(box_res$parameter),
      "Estat\u00edstica Q"   = .fmt_num(box_res$statistic, 2),
      "p-valor"         = fmt_pbox,
      "Conclus\u00e3o"        = if (p_box >= 0.05) "Ru\u00eddo Branco" else "Autocorrelacionada",
      stringsAsFactors  = FALSE, check.names = FALSE
    )

    .print_topico("RESULTADO")
    if (exists(".print_tabela_estatR", mode = "function")) {
      .print_tabela_estatR(df_box, align = c("center", "center", "center", "center"))
    } else {
      print(df_box, row.names = FALSE)
    }
    cat("\n")

    .print_topico("INTERPRETA\u00c7\u00c3O")
    cor_ok   <- "\033[1;32m"
    cor_erro <- "\033[1;31m"
    reset    <- "\033[0m"

    if (p_box >= 0.05) {
      cat(sprintf("  %s[\u2713]%s A s\u00e9rie se comporta como RU\u00cdDO BRANCO (p \u2265 0,05).\n", cor_ok, reset))
      cat("      N\u00e3o h\u00e1 evid\u00eancia estat\u00edstica de autocorrela\u00e7\u00e3o linear nas defasagens testadas.\n")
    } else {
      cat(sprintf("  %s[!]%s Presen\u00e7a de AUTOCORRELA\u00c7\u00c3O SIGNIFICATIVA (p < 0,05).\n", cor_erro, reset))
      cat("      A s\u00e9rie cont\u00e9m padr\u00e3o temporal dependente pass\u00edvel de modelagem (AR/MA/SARIMA).\n")
    }

    .print_rodape()
  }

  print(p_final)
  invisible(p_final)
}

#' @rdname grafico_autocorrelacao
#' @export
grafico_acf <- grafico_autocorrelacao


# ── 5. MOTOR DE MÉTRICAS DE SÉRIES TEMPORAIS ────────────────────────────────

#' @title Métricas de Acurácia para Séries Temporais
#' @description Calcula as principais métricas de erro preditivo e critérios de informação
#'   para séries temporais: RMSE, MAE, MAPE, sMAPE, MASE, UMBRAE, AIC, AICc e BIC.
#'
#' @param real Vetor numérico com os valores reais observados (teste) ou modelo ajustado.
#' @param previsto Vetor numérico com os valores preditos pelo modelo.
#' @param treino Vetor numérico opcional com a série histórica de treino (necessário para escala do MASE).
#' @param frequencia Periodicidade da série para o cálculo sazonal do MASE. Padrão: 12.
#'
#' @return Data frame com as métricas calculadas.
#' @export
metricas_serie <- function(real, previsto = NULL, treino = NULL, frequencia = 12) {
  # Se o usuário passou um objeto de modelo de previsão
  if (inherits(real, c("forecast", "previsao_serie", "Arima", "ets", "HoltWinters", "stlm"))) {
    mod <- real

    # previsao_serie: contém $previsao (objeto forecast) e $serie (ts original)
    if (inherits(mod, "previsao_serie")) {
      fc_obj    <- mod$previsao
      real_vals <- as.numeric(fc_obj$x)
      prev_vals <- as.numeric(fc_obj$fitted)
      treino    <- real_vals
      freq_ef   <- stats::frequency(mod$serie)
      return(.calcular_metricas_ts(real_vals, prev_vals, treino, freq_ef, modelo_obj = mod$modelo))
    }

    # Objeto forecast direto (forecast::forecast, ses, holt, hw, etc.)
    if (inherits(mod, "forecast")) {
      real_vals <- as.numeric(mod$x)
      prev_vals <- as.numeric(mod$fitted)
      treino    <- real_vals
      freq_ef   <- if (!is.null(mod$x) && stats::is.ts(mod$x)) stats::frequency(mod$x) else frequencia
      return(.calcular_metricas_ts(real_vals, prev_vals, treino, freq_ef, modelo_obj = mod$model))
    }

    # Objetos de modelo direto (ets, Arima, HoltWinters, stlm)
    real_vals <- tryCatch(as.numeric(stats::fitted(mod)) + as.numeric(stats::residuals(mod)), error = function(e) NULL)
    prev_vals <- tryCatch(as.numeric(stats::fitted(mod)), error = function(e) NULL)
    if (is.null(real_vals) || is.null(prev_vals)) {
      stop("N\u00e3o foi poss\u00edvel extrair os valores ajustados do objeto fornecido.")
    }
    treino  <- real_vals
    freq_ef <- if (stats::is.ts(prev_vals)) stats::frequency(prev_vals) else frequencia
    return(.calcular_metricas_ts(real_vals, prev_vals, treino, freq_ef, modelo_obj = mod))
  }

  if (is.null(previsto)) {
    stop("O argumento 'previsto' deve ser fornecido quando 'real' for um vetor num\u00e9rico.")
  }

  .calcular_metricas_ts(as.numeric(real), as.numeric(previsto), treino, frequencia)
}

#' @noRd
.calcular_metricas_ts <- function(y_real, y_hat, y_train = NULL, m = 12, modelo_obj = NULL) {
  # Limpeza de NAs emparelhados
  valid <- !is.na(y_real) & !is.na(y_hat)
  y_real <- y_real[valid]
  y_hat  <- y_hat[valid]
  h <- length(y_real)

  if (h == 0) stop("Nenhum par de valores v\u00e1lidos entre 'real' e 'previsto'.")

  # 1. RMSE & MAE
  erros <- y_real - y_hat
  rmse  <- sqrt(mean(erros^2))
  mae   <- mean(abs(erros))

  # 2. MAPE & sMAPE
  # Proteção contra divisão por zero no MAPE
  idx_nz <- y_real != 0
  mape   <- if (any(idx_nz)) mean(abs(erros[idx_nz] / y_real[idx_nz])) * 100 else NA_real_
  denom_smape <- abs(y_real) + abs(y_hat)
  smape  <- if (any(denom_smape > 0)) mean(200 * abs(erros) / denom_smape) else 0

  # 3. MASE (Mean Absolute Scaled Error)
  mase <- NA_real_
  if (!is.null(y_train) && length(y_train) > m) {
    # Erro de escala ingênuo sazonal no treino
    diff_sazonal <- abs(base::diff(y_train, lag = m))
    escala_treino <- mean(diff_sazonal, na.rm = TRUE)
    if (!is.na(escala_treino) && escala_treino > 0) {
      mase <- mae / escala_treino
    }
  }

  # 4. UMBRAE (Unscaled Bounded Relative Absolute Error)
  # Benchmark ingênuo: último valor de treino se fornecido, ou lag 1
  y_naive <- if (!is.null(y_train) && length(y_train) >= 1) {
    rep(tail(y_train, 1), h)
  } else {
    c(y_real[1], y_real[-h])
  }

  e_naive <- abs(y_real - y_naive)
  brae <- abs(erros) / (abs(erros) + e_naive)
  # Trata 0/0 quando ambos acertam perfeitamente
  brae[is.na(brae)] <- 0.5
  mbrae <- mean(brae)
  umbrae <- if (!is.na(mbrae) && mbrae < 1) mbrae / (1 - mbrae) else Inf

  # 5. Critérios de Parcimônia (AIC, AICc, BIC)
  aic_val  <- NA_real_
  aicc_val <- NA_real_
  bic_val  <- NA_real_

  if (!is.null(modelo_obj)) {
    aic_val <- tryCatch(stats::AIC(modelo_obj), error = function(e) NA_real_)
    bic_val <- tryCatch(stats::BIC(modelo_obj), error = function(e) NA_real_)
    if (!is.null(modelo_obj$aicc)) {
      aicc_val <- modelo_obj$aicc
    } else if (!is.na(aic_val)) {
      k_pars <- tryCatch(length(stats::coef(modelo_obj)), error = function(e) 2)
      n_tot  <- length(y_real)
      if (n_tot - k_pars - 1 > 0) {
        aicc_val <- aic_val + (2 * k_pars * (k_pars + 1)) / (n_tot - k_pars - 1)
      }
    }
  }

  df_res <- data.frame(
    "M\u00e9trica" = c("RMSE", "MAE", "MAPE (%)", "sMAPE (%)", "MASE", "UMBRAE", "AIC", "AICc", "BIC"),
    Valor   = c(
      .fmt_num(rmse, 2),
      .fmt_num(mae, 2),
      if (is.na(mape)) "-" else paste0(.fmt_num(mape, 1), "%"),
      paste0(.fmt_num(smape, 1), "%"),
      if (is.na(mase)) "-" else .fmt_num(mase, 3),
      if (is.infinite(umbrae) || is.na(umbrae)) "-" else .fmt_num(umbrae, 3),
      if (is.na(aic_val)) "-" else .fmt_num(aic_val, 1),
      if (is.na(aicc_val)) "-" else .fmt_num(aicc_val, 1),
      if (is.na(bic_val)) "-" else .fmt_num(bic_val, 1)
    ),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  .print_titulo("M\u00c9TRICAS DE ACUR\u00c1CIA DA S\u00c9RIE")
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_res, align = c("left", "center"))
  } else {
    print(df_res, row.names = FALSE)
  }
  cat("  * UMBRAE < 1 indica desempenho superior ao benchmark ing\u00eano (Naive).\n")
  cat("  * MASE < 1 indica erro m\u00e9dio menor que o passeio aleat\u00f3rio sazonal em treino.\n\n")
  .print_rodape()

  invisible(list(
    RMSE = rmse, MAE = mae, MAPE = mape, sMAPE = smape,
    MASE = mase, UMBRAE = umbrae, AIC = aic_val, AICc = aicc_val, BIC = bic_val
  ))
}


# ── 6. PREVISÃO E PROJEÇÃO FUTURA (EXPONENTIAL SMOOTHING / ETS) ──────────────

#' @title Previsão de Série Temporal (Suavização Exponencial e ETS)
#' @description Modela a série temporal histórica por Suavização Exponencial
#'   (SES, Holt Linear, Holt Damped ou Holt-Winters) e projeta previsões
#'   futuras com intervalos de confiança em formato de leque (*Fan Chart*).
#'
#' @param serie Vetor numérico, objeto \code{ts} ou data frame.
#' @param horizonte Número de períodos a projetar para o futuro (padrão: 12).
#' @param frequencia Periodicidade da série. Padrão: 12.
#' @param metodo Método de suavização: \code{"auto"} (seleção ótima via menor AICc),
#'   \code{"holt_winters"}, \code{"holt"}, \code{"ses"} ou \code{"damped"}.
#' @param tipo Decomposição sazonal: \code{"multiplicativa"} ou \code{"aditiva"}.
#' @param confianca Nível do intervalo de confiança para as projeções (padrão: 0.95).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o gráfico das previsões.
#'
#' @return Objeto invisível da classe \code{previsao_serie} contendo o modelo e as projeções.
#' @export
previsao_serie <- function(serie, horizonte = 12, frequencia = 12,
                           metodo = c("auto", "holt_winters", "holt", "ses", "damped"),
                           tipo = c("multiplicativa", "aditiva"),
                           confianca = 0.95, grafico = TRUE) {

  metodo <- match.arg(metodo)
  tipo   <- match.arg(tipo)

  ts_serie <- .extrair_serie_temporal(serie, frequencia = frequencia)
  n_obs <- length(ts_serie)
  freq_ef <- stats::frequency(ts_serie)

  # Ajuste do modelo
  modelo <- tryCatch({
    if (metodo == "auto") {
      forecast::ets(ts_serie)
    } else if (metodo == "ses") {
      forecast::ses(ts_serie, h = horizonte, level = c(80, confianca * 100))
    } else if (metodo == "holt") {
      forecast::holt(ts_serie, h = horizonte, level = c(80, confianca * 100))
    } else if (metodo == "damped") {
      forecast::holt(ts_serie, damped = TRUE, h = horizonte, level = c(80, confianca * 100))
    } else if (metodo == "holt_winters") {
      saz_tipo <- if (tipo == "multiplicativa") "multiplicative" else "additive"
      if (any(ts_serie <= 0) && saz_tipo == "multiplicative") saz_tipo <- "additive"
      forecast::hw(ts_serie, seasonal = saz_tipo, h = horizonte, level = c(80, confianca * 100))
    }
  }, error = function(e) {
    # Fallback para ETS automático em caso de convergência difícil
    forecast::ets(ts_serie)
  })

  # Se o modelo foi gerado via ets puro, gera o forecast
  fc <- if (inherits(modelo, "forecast")) {
    modelo
  } else {
    forecast::forecast(modelo, h = horizonte, level = c(80, confianca * 100))
  }

  # Extração de projeções
  pred_media <- as.numeric(fc$mean)
  pred_inf80 <- as.numeric(fc$lower[, 1])
  pred_sup80 <- as.numeric(fc$upper[, 1])
  pred_inf95 <- as.numeric(fc$lower[, 2])
  pred_sup95 <- as.numeric(fc$upper[, 2])

  # Construção de tempos futuros
  tempo_hist <- as.numeric(stats::time(ts_serie))
  delta_t <- 1 / freq_ef
  ultimo_t <- tail(tempo_hist, 1)
  tempo_fut  <- ultimo_t + seq(1, horizonte) * delta_t

  # Rótulos para a tabela de projeção
  rotulos_futuros <- paste0("t+", seq_len(horizonte))

  df_proj_tab <- data.frame(
    "Per\u00edodo"        = rotulos_futuros,
    "Previs\u00e3o"       = sapply(pred_media, .fmt_num, decimais = 2),
    "Limite Inferior" = sapply(pred_inf95, .fmt_num, decimais = 2),
    "Limite Superior" = sapply(pred_sup95, .fmt_num, decimais = 2),
    stringsAsFactors  = FALSE, check.names = FALSE
  )

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("PREVIS\u00c3O DE S\u00c9RIE TEMPORAL")

  nome_metodo <- if (inherits(fc$model, "ets")) {
    paste0("ETS (", fc$model$method, ")")
  } else {
    fc$method
  }

  cat(sprintf("  M\u00e9todo de ajuste:   %s\n", nome_metodo))
  cat(sprintf("  Horizonte de proj: %d per\u00edodos\n", horizonte))
  cat(sprintf("  Intervalo de Conf: %.0f%%\n\n", confianca * 100))

  .print_topico(sprintf("PROJE\u00c7\u00d5ES FUTURAS (%d PER\u00cdODOS)", min(horizonte, 12)))

  tab_exibir <- if (nrow(df_proj_tab) > 12) df_proj_tab[1:12, ] else df_proj_tab
  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(tab_exibir, align = c("center", "center", "center", "center"))
  } else {
    print(tab_exibir, row.names = FALSE)
  }

  if (nrow(df_proj_tab) > 12) {
    cat(sprintf("  * Mostrando as primeiras 12 de %d proje\u00e7\u00f5es geradas.\n", horizonte))
  }
  cat("\n")
  .print_rodape()

  # ── Painel Gráfico (Fan Chart) ─────────────────────────────────────────────
  p_graf <- NULL
  if (isTRUE(grafico)) {
    # Data frames para o gráfico
    df_hist <- data.frame(Tempo = tempo_hist, Valor = as.numeric(ts_serie))

    df_fut <- data.frame(
      Tempo = tempo_fut,
      Media = pred_media,
      Inf80 = pred_inf80,
      Sup80 = pred_sup80,
      Inf95 = pred_inf95,
      Sup95 = pred_sup95
    )

    # Conectar histórico com futuro suavemente
    ponto_conexao <- data.frame(
      Tempo = ultimo_t,
      Media = tail(as.numeric(ts_serie), 1),
      Inf80 = tail(as.numeric(ts_serie), 1),
      Sup80 = tail(as.numeric(ts_serie), 1),
      Inf95 = tail(as.numeric(ts_serie), 1),
      Sup95 = tail(as.numeric(ts_serie), 1)
    )
    df_fut_con <- rbind(ponto_conexao, df_fut)

    p_graf <- ggplot2::ggplot() +
      # Faixas de confiança em camadas (Fan Chart)
      ggplot2::geom_ribbon(data = df_fut_con, ggplot2::aes(x = Tempo, ymin = Inf95, ymax = Sup95),
                           fill = "#377EB8", alpha = 0.20) +
      ggplot2::geom_ribbon(data = df_fut_con, ggplot2::aes(x = Tempo, ymin = Inf80, ymax = Sup80),
                           fill = "#377EB8", alpha = 0.35) +
      # Linha histórica
      ggplot2::geom_line(data = df_hist, ggplot2::aes(x = Tempo, y = Valor),
                         color = "gray25", linewidth = 0.8) +
      # Linha de previsão pontual
      ggplot2::geom_line(data = df_fut_con, ggplot2::aes(x = Tempo, y = Media),
                         color = "#377EB8", linewidth = 1.05) +
      ggplot2::geom_point(data = df_fut, ggplot2::aes(x = Tempo, y = Media),
                          color = "#377EB8", size = 1.8) +
      ggplot2::labs(
        title = paste("Previs\u00e3o da S\u00e9rie Temporal -", nome_metodo),
        subtitle = sprintf("Horizonte: %d per\u00edodos | Faixas de Confian\u00e7a: 80%% e %.0f%%", horizonte, confianca * 100),
        x = "Tempo", y = "Valor"
      ) +
      tema_estatR(estilo = 2, grade = "dupla")

    print(p_graf)
  }

  res_obj <- structure(
    list(
      modelo     = fc$model,
      previsao   = fc,
      projecoes  = df_proj_tab,
      grafico    = p_graf,
      serie      = ts_serie
    ),
    class = c("previsao_serie", "forecast")
  )

  invisible(res_obj)
}


# ── 7. SELEÇÃO E COMPARAÇÃO DE MODELOS (OUT-OF-SAMPLE RANKING) ──────────────

#' @title Seleção de Modelos de Séries Temporais
#' @description Compara exaustivamente modelos preditivos organizados por famílias
#'   (\code{"suavizacao"}, \code{"estocastica"}, \code{"decomposicao"} ou \code{"todas"})
#'   utilizando divisão real de validação fora da amostra (out-of-sample) de tamanho \code{h}.
#'   Retorna uma tabela formatada com o ranking dos melhores modelos segundo métricas de acurácia.
#'
#' @param serie Vetor numérico, objeto \code{ts} ou data frame contendo a série temporal.
#' @param h Tamanho da janela de validação fora da amostra (horizonte de teste out-of-sample). Padrão: 12.
#' @param frequencia Periodicidade da série (ex: 12 para mensal). Padrão: 12.
#' @param familia Família de modelos a avaliar: \code{"todas"} (padrão), \code{"suavizacao"},
#'   \code{"estocastica"} ou \code{"decomposicao"}.
#' @param top Número de modelos a exibir no ranking (padrão: 5).
#' @param ordenar_por Critério de ordenação do ranking: \code{"RMSE"} (padrão), \code{"sMAPE"},
#'   \code{"UMBRAE"}, \code{"AICc"} ou \code{"BIC"}.
#'
#' @return Retorna invisivelmente um data frame com o ranking completo e as métricas calculadas.
#' @export
selecao_modelos_serie <- function(serie, h = 12, frequencia = 12,
                                  familia = c("todas", "suavizacao", "estocastica", "decomposicao"),
                                  top = 5,
                                  ordenar_por = c("RMSE", "sMAPE", "UMBRAE", "AICc", "BIC")) {

  familia     <- match.arg(familia)
  ordenar_por <- match.arg(ordenar_por)

  ts_serie <- .extrair_serie_temporal(serie, frequencia = frequencia)
  n_total  <- length(ts_serie)
  freq_ef  <- stats::frequency(ts_serie)

  if (n_total <= h + 4) {
    stop(sprintf("A s\u00e9rie cont\u00e9m %d observa\u00e7\u00f5es, o que \u00e9 insuficiente para uma divis\u00e3o de teste h = %d (m\u00ednimo: %d).",
                 n_total, h, h + 5))
  }

  # Divisão Treino e Teste (Holdout out-of-sample)
  n_treino <- n_total - h
  y_treino <- stats::ts(as.numeric(ts_serie)[1:n_treino],
                        start = stats::start(ts_serie), frequency = freq_ef)
  y_teste  <- as.numeric(ts_serie)[(n_treino + 1):n_total]

  # Benchmark Naive sazonal para escala e UMBRAE
  y_naive_bench <- if (n_treino >= freq_ef) {
    as.numeric(tail(y_treino, freq_ef))[((seq_len(h) - 1) %% freq_ef) + 1]
  } else {
    rep(tail(y_treino, 1), h)
  }

  # Lista de modelos avaliados
  modelos_lista <- list()

  # ── Família 1: Suavização Exponencial ─────────────────────────────────────
  if (familia %in% c("todas", "suavizacao")) {
    # 1. SES
    tryCatch({
      m_ses <- forecast::ses(y_treino, h = h)
      modelos_lista[["SES"]] <- list(nome = "SES (Suav. Simples)", fam = "Suaviza\u00e7\u00e3o", fc = m_ses, mod = m_ses$model)
    }, error = function(e) NULL)

    # 2. Holt Linear
    tryCatch({
      m_holt <- forecast::holt(y_treino, h = h)
      modelos_lista[["Holt"]] <- list(nome = "Holt Linear", fam = "Suaviza\u00e7\u00e3o", fc = m_holt, mod = m_holt$model)
    }, error = function(e) NULL)

    # 3. Holt Damped
    tryCatch({
      m_damped <- forecast::holt(y_treino, damped = TRUE, h = h)
      modelos_lista[["Holt_Damped"]] <- list(nome = "Holt Amortecido (Damped)", fam = "Suaviza\u00e7\u00e3o", fc = m_damped, mod = m_damped$model)
    }, error = function(e) NULL)

    # 4. Holt-Winters Aditivo
    if (n_treino >= 2 * freq_ef) {
      tryCatch({
        m_hwa <- forecast::hw(y_treino, seasonal = "additive", h = h)
        modelos_lista[["HW_Aditivo"]] <- list(nome = "Holt-Winters Aditivo", fam = "Suaviza\u00e7\u00e3o", fc = m_hwa, mod = m_hwa$model)
      }, error = function(e) NULL)

      # 5. Holt-Winters Multiplicativo
      if (all(y_treino > 0)) {
        tryCatch({
          m_hwm <- forecast::hw(y_treino, seasonal = "multiplicative", h = h)
          modelos_lista[["HW_Multiplicativo"]] <- list(nome = "Holt-Winters Multiplicativo", fam = "Suaviza\u00e7\u00e3o", fc = m_hwm, mod = m_hwm$model)
        }, error = function(e) NULL)
      }
    }

    # 6. ETS Automático
    tryCatch({
      m_ets <- forecast::ets(y_treino)
      fc_ets <- forecast::forecast(m_ets, h = h)
      modelos_lista[["ETS_Auto"]] <- list(nome = paste0("ETS (", m_ets$method, ")"), fam = "Suaviza\u00e7\u00e3o", fc = fc_ets, mod = m_ets)
    }, error = function(e) NULL)
  }

  # ── Família 2: Estocástica (ARIMA / SARIMA) ───────────────────────────────
  if (familia %in% c("todas", "estocastica")) {
    # 1. Auto ARIMA / SARIMA
    tryCatch({
      m_arima <- forecast::auto.arima(y_treino, seasonal = (freq_ef > 1), stepwise = TRUE, approximation = FALSE)
      fc_arima <- forecast::forecast(m_arima, h = h)
      nome_ar <- paste0("ARIMA", gsub("ARIMA", "", forecast::arima.string(m_arima)))
      modelos_lista[["Auto_ARIMA"]] <- list(nome = nome_ar, fam = "Estoc\u00e1stica", fc = fc_arima, mod = m_arima)
    }, error = function(e) NULL)

    # 2. ARIMA(0,1,1)
    tryCatch({
      m_ima <- stats::arima(y_treino, order = c(0, 1, 1))
      fc_ima <- forecast::forecast(m_ima, h = h)
      modelos_lista[["IMA_011"]] <- list(nome = "ARIMA(0,1,1)", fam = "Estoc\u00e1stica", fc = fc_ima, mod = m_ima)
    }, error = function(e) NULL)

    # 3. ARIMA(1,1,0)
    tryCatch({
      m_ari <- stats::arima(y_treino, order = c(1, 1, 0))
      fc_ari <- forecast::forecast(m_ari, h = h)
      modelos_lista[["ARI_110"]] <- list(nome = "ARIMA(1,1,0)", fam = "Estoc\u00e1stica", fc = fc_ari, mod = m_ari)
    }, error = function(e) NULL)

    # 4. ARIMA(1,1,1)
    tryCatch({
      m_arima111 <- stats::arima(y_treino, order = c(1, 1, 1))
      fc_arima111 <- forecast::forecast(m_arima111, h = h)
      modelos_lista[["ARIMA_111"]] <- list(nome = "ARIMA(1,1,1)", fam = "Estoc\u00e1stica", fc = fc_arima111, mod = m_arima111)
    }, error = function(e) NULL)
  }

  # ── Família 3: Decomposição (STL + Modelos) ───────────────────────────────
  if (familia %in% c("todas", "decomposicao")) {
    # 1. Naive (Passeio Aleatório)
    tryCatch({
      m_naive <- forecast::naive(y_treino, h = h)
      modelos_lista[["Naive"]] <- list(nome = "Passeio Aleat\u00f3rio (Naive)", fam = "Decomposi\u00e7\u00e3o", fc = m_naive, mod = NULL)
    }, error = function(e) NULL)

    # 2. SNaive (Naive Sazonal)
    if (freq_ef > 1 && n_treino >= freq_ef) {
      tryCatch({
        m_snaive <- forecast::snaive(y_treino, h = h)
        modelos_lista[["SNaive"]] <- list(nome = "Naive Sazonal (SNaive)", fam = "Decomposi\u00e7\u00e3o", fc = m_snaive, mod = NULL)
      }, error = function(e) NULL)
    }

    # 3. STL + ETS
    if (freq_ef > 1 && n_treino >= 2 * freq_ef) {
      tryCatch({
        m_stlets <- forecast::stlm(y_treino, method = "ets")
        fc_stlets <- forecast::forecast(m_stlets, h = h)
        modelos_lista[["STL_ETS"]] <- list(nome = "STL + ETS", fam = "Decomposi\u00e7\u00e3o", fc = fc_stlets, mod = m_stlets$model)
      }, error = function(e) NULL)

      # 4. STL + ARIMA
      tryCatch({
        m_stlarima <- forecast::stlm(y_treino, method = "arima")
        fc_stlarima <- forecast::forecast(m_stlarima, h = h)
        modelos_lista[["STL_ARIMA"]] <- list(nome = "STL + ARIMA", fam = "Decomposi\u00e7\u00e3o", fc = fc_stlarima, mod = m_stlarima$model)
      }, error = function(e) NULL)
    }
  }

  if (length(modelos_lista) == 0) {
    stop("Nenhum modelo p\u00f4de ser ajustado com sucesso. Verifique os dados e par\u00e2metros.")
  }

  # ── Cálculo das Métricas de Validação Out-of-Sample ─────────────────────────
  nomes_vet    <- character(length(modelos_lista))
  fams_vet     <- character(length(modelos_lista))
  rmse_vet     <- numeric(length(modelos_lista))
  smape_vet    <- numeric(length(modelos_lista))
  umbrae_vet   <- numeric(length(modelos_lista))
  aicc_vet     <- numeric(length(modelos_lista))
  bic_vet      <- numeric(length(modelos_lista))

  for (i in seq_along(modelos_lista)) {
    item <- modelos_lista[[i]]
    y_prev <- as.numeric(item$fc$mean)

    # 1. RMSE
    erros <- y_teste - y_prev
    rmse_vet[i] <- sqrt(mean(erros^2, na.rm = TRUE))

    # 2. sMAPE (%)
    denom <- abs(y_teste) + abs(y_prev)
    smape_vet[i] <- mean(200 * abs(erros) / denom, na.rm = TRUE)

    # 3. UMBRAE
    e_mod   <- abs(erros)
    e_bench <- abs(y_teste - y_naive_bench)
    brae    <- e_mod / (e_mod + e_bench)
    brae[is.na(brae)] <- 0.5
    mbrae   <- mean(brae)
    umbrae_vet[i] <- if (!is.na(mbrae) && mbrae < 1) mbrae / (1 - mbrae) else 1.0

    # 4. AICc e BIC
    mod_obj <- item$mod
    aic_val <- if (!is.null(mod_obj)) tryCatch(stats::AIC(mod_obj), error = function(e) NA_real_) else NA_real_
    bic_vet[i] <- if (!is.null(mod_obj)) tryCatch(stats::BIC(mod_obj), error = function(e) NA_real_) else NA_real_

    if (!is.null(mod_obj) && !is.null(mod_obj$aicc)) {
      aicc_vet[i] <- mod_obj$aicc
    } else if (!is.na(aic_val)) {
      k_p <- tryCatch(length(stats::coef(mod_obj)), error = function(e) 2)
      if (n_treino - k_p - 1 > 0) {
        aicc_vet[i] <- aic_val + (2 * k_p * (k_p + 1)) / (n_treino - k_p - 1)
      } else {
        aicc_vet[i] <- aic_val
      }
    } else {
      aicc_vet[i] <- NA_real_
    }

    nomes_vet[i] <- item$nome
    fams_vet[i]  <- item$fam
  }

  df_rank <- data.frame(
    Modelo  = nomes_vet,
    Familia = fams_vet,
    RMSE    = rmse_vet,
    sMAPE   = smape_vet,
    UMBRAE  = umbrae_vet,
    AICc    = aicc_vet,
    BIC     = bic_vet,
    stringsAsFactors = FALSE
  )

  # Ordenação
  ord_idx <- order(df_rank[[ordenar_por]], decreasing = FALSE, na.last = TRUE)
  df_rank <- df_rank[ord_idx, ]

  # Slicing TOP
  df_top <- if (nrow(df_rank) > top) df_rank[1:top, ] else df_rank

  # Tabela formatada para exibição
  df_fmt <- data.frame(
    Rank        = paste0(seq_len(nrow(df_top)), "\u00ba"),
    Modelo      = df_top$Modelo,
    "Fam\u00edlia"     = df_top$Familia,
    RMSE        = sapply(df_top$RMSE, .fmt_num, decimais = 2),
    "sMAPE (%)" = sapply(df_top$sMAPE, function(v) paste0(.fmt_num(v, 1), "%")),
    UMBRAE      = sapply(df_top$UMBRAE, .fmt_num, decimais = 3),
    AICc        = sapply(df_top$AICc, function(v) if (is.na(v)) "-" else .fmt_num(v, 1)),
    BIC         = sapply(df_top$BIC, function(v) if (is.na(v)) "-" else .fmt_num(v, 1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("SELE\u00c7\u00c3O DE MODELOS DE S\u00c9RIES TEMPORAIS")

  cat(sprintf("  Amostra Total:    %d observa\u00e7\u00f5es\n", n_total))
  cat(sprintf("  Treino (In-Sample):   %d per\u00edodos\n", n_treino))
  cat(sprintf("  Valida\u00e7\u00e3o (Holdout):  %d per\u00edodos (h)\n", h))
  cat(sprintf("  Fam\u00edlia Avaliada:   %s\n", if (familia == "todas") "Todas as Fam\u00edlias" else df_top$Familia[1]))
  cat(sprintf("  Crit\u00e9rio de Ordem:   %s (Menor \u00e9 melhor)\n\n", ordenar_por))

  .print_topico(sprintf("TOP %d MODELOS (Ordenados por %s)", nrow(df_fmt), ordenar_por))

  if (exists(".print_tabela_estatR", mode = "function")) {
    # Alinhamento estritamente centralizado com os cabeçalhos
    .print_tabela_estatR(df_fmt, align = c("center", "left", "center", "center", "center", "center", "center", "center"))
  } else {
    print(df_fmt, row.names = FALSE)
  }

  cat("  * UMBRAE < 1 indica desempenho superior ao benchmark ing\u00eano (Naive).\n")
  cat("  * Avalia\u00e7\u00e3o baseada em previs\u00e3o real fora da amostra (out-of-sample).\n\n")

  .print_rodape()

  invisible(df_rank)
}


# ── 8. ANÁLISE RESIDUAL DE SÉRIES TEMPORAIS ─────────────────────────────────

#' @title Análise Residual para Séries Temporais
#' @description Realiza o diagnóstico formal dos resíduos de modelos de séries temporais
#'   (\code{previsao_serie}, \code{forecast}, \code{Arima}, \code{ets}, \code{stlm}).
#'   Avalia normalidade (Shapiro-Wilk ou Kolmogorov-Smirnov), ausência de viés (teste t para média zero),
#'   independência e ruído branco (teste de Ljung-Box em defasagens múltiplas) e
#'   gera um painel gráfico com 4 diagnósticos visuais.
#'
#' @param modelo Objeto ajustado da classe \code{previsao_serie}, \code{forecast},
#'   \code{Arima}, \code{ets}, \code{HoltWinters} ou \code{stlm}.
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o painel gráfico de diagnóstico.
#'
#' @return Retorna invisivelmente uma lista com os objetos dos testes estatísticos e resíduos.
#' @export
analise_residual_serie <- function(modelo, grafico = TRUE) {

  # Extração dos resíduos conforme a classe do objeto
  res <- if (inherits(modelo, "previsao_serie")) {
    if (!is.null(modelo$previsao$residuals)) {
      as.numeric(modelo$previsao$residuals)
    } else {
      as.numeric(stats::residuals(modelo$modelo))
    }
  } else if (inherits(modelo, "forecast")) {
    as.numeric(modelo$residuals)
  } else {
    as.numeric(stats::residuals(modelo))
  }

  res <- res[!is.na(res)]
  n   <- length(res)

  if (n < 4) {
    stop("Quantidade insuficiente de res\u00edduos v\u00e1lidos para realizar a an\u00e1lise (m\u00ednimo 4).")
  }

  nome_modelo <- if (inherits(modelo, "previsao_serie")) {
    if (inherits(modelo$previsao$model, "ets")) {
      paste0("ETS (", modelo$previsao$model$method, ")")
    } else {
      modelo$previsao$method
    }
  } else if (inherits(modelo, "forecast")) {
    modelo$method
  } else if (inherits(modelo, "Arima")) {
    a <- modelo$arma
    if (!is.null(a) && length(a) >= 7) {
      if (a[7] > 0) {
        sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[%d]", a[1], a[6], a[2], a[3], a[7], a[4], a[5])
      } else {
        sprintf("ARIMA(%d,%d,%d)", a[1], a[6], a[2])
      }
    } else {
      "ARIMA"
    }
  } else if (inherits(modelo, "ets")) {
    paste0("ETS (", modelo$method, ")")
  } else {
    class(modelo)[1]
  }

  # 1. Média Zero / Ausência de Viés (Teste t para H0: mu = 0)
  t_test <- stats::t.test(res, mu = 0)
  stat_t <- t_test$statistic
  p_t    <- t_test$p.value

  # 2. Normalidade (Shapiro-Wilk para n <= 5000; Kolmogorov-Smirnov para n > 5000)
  if (n <= 5000) {
    nome_norm <- "Normalidade (Shapiro-Wilk)"
    st_norm   <- stats::shapiro.test(res)
    stat_norm <- st_norm$statistic
    p_norm    <- st_norm$p.value
    txt_stat_norm <- sprintf("W = %.3f", stat_norm)
  } else {
    nome_norm <- "Normalidade (Kolmogorov-S)"
    res_pad   <- (res - mean(res)) / stats::sd(res)
    st_norm   <- suppressWarnings(stats::ks.test(res_pad, "pnorm"))
    stat_norm <- st_norm$statistic
    p_norm    <- st_norm$p.value
    txt_stat_norm <- sprintf("D = %.3f", stat_norm)
  }

  # 3. Independência e Ruído Branco (Ljung-Box)
  # Seleciona defasagens recomendadas pela literatura: min(10 ou 2m, n-2)
  lag_box <- min(12, n - 2)
  st_box  <- stats::Box.test(res, lag = lag_box, type = "Ljung-Box")
  stat_box <- st_box$statistic
  p_box   <- st_box$p.value

  # Formatação dos p-valores
  fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }

  df_tab <- data.frame(
    Pressuposto = c("M\u00e9dia Zero (Sem vi\u00e9s)", nome_norm, sprintf("Ru\u00eddo Branco (Ljung-Box, lag %d)", lag_box)),
    "Estat\u00edstica" = c(sprintf("t = %.3f", stat_t), txt_stat_norm, sprintf("Q = %.2f", stat_box)),
    "p-valor"   = c(fmt_p(p_t), fmt_p(p_norm), fmt_p(p_box)),
    "Diagn\u00f3stico" = c(
      if (p_t >= 0.05) "Sem Vi\u00e9s" else "Vi\u00e9s Detectado",
      if (p_norm >= 0.05) "Normal" else "N\u00e3o Normal",
      if (p_box >= 0.05) "Ru\u00eddo Branco" else "Autocorrelacionado"
    ),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # ── Saída no Console ───────────────────────────────────────────────────────
  .print_titulo("DIAGN\u00d3STICO DOS RES\u00cdDUOS DE S\u00c9RIE TEMPORAL")

  cat(sprintf("  Modelo avaliado:  %s\n", nome_modelo))
  cat(sprintf("  Res\u00edduos totais:  %d observa\u00e7\u00f5es\n", n))
  cat(sprintf("  M\u00e9dia residual:   %s (Desvio: %s)\n\n",
              .fmt_num(mean(res), 4), .fmt_num(stats::sd(res), 4)))

  .print_topico("RESULTADOS DOS TESTES DE DIAGN\u00d3STICO")

  if (exists(".print_tabela_estatR", mode = "function")) {
    .print_tabela_estatR(df_tab, align = c("left", "center", "center", "center"))
  } else {
    print(df_tab, row.names = FALSE)
  }
  cat("\n")

  .print_topico("INTERPRETA\u00c7\u00c3O")

  cor_ok   <- "\033[1;32m"
  cor_erro <- "\033[1;31m"
  reset    <- "\033[0m"

  # Interpretação da Média Zero
  if (p_t >= 0.05) {
    cat(sprintf("  %s[\u2713]%s A m\u00e9dia dos res\u00edduos n\u00e3o difere de zero (p = %s). O modelo n\u00e3o apresenta vi\u00e9s sistem\u00e1tico.\n",
                cor_ok, reset, fmt_p(p_t)))
  } else {
    cat(sprintf("  %s[!]%s Vi\u00e9s detectado nos res\u00edduos (p = %s). As previs\u00f5es tendem a subestimar ou superestimar a s\u00e9rie.\n",
                cor_erro, reset, fmt_p(p_t)))
  }

  # Interpretação da Normalidade
  if (p_norm >= 0.05) {
    cat(sprintf("  %s[\u2713]%s Os res\u00edduos seguem distribui\u00e7\u00e3o normal (p = %s), validando os intervalos de confian\u00e7a anal\u00edticos.\n",
                cor_ok, reset, fmt_p(p_norm)))
  } else {
    cat(sprintf("  %s[!]%s Os res\u00edduos n\u00e3o s\u00e3o estritamente normais (p = %s). Recomenda-se cautela ou intervalos bootstrap.\n",
                cor_erro, reset, fmt_p(p_norm)))
  }

  # Interpretação do Ruído Branco
  if (p_box >= 0.05) {
    cat(sprintf("  %s[\u2713]%s Res\u00edduos comportam-se como RU\u00cdDO BRANCO (p = %s, lag %d). Toda a depend\u00eancia temporal foi capturada.\n",
                cor_ok, reset, fmt_p(p_box), lag_box))
  } else {
    cat(sprintf("  %s[!]%s Autocorrela\u00e7\u00e3o residual remanescente detectada (p = %s, lag %d). Considere ajustar ordens adicionais de AR/MA.\n",
                cor_erro, reset, fmt_p(p_box), lag_box))
  }

  .print_rodape()

  # ── Painel Gráfico de 4 Resíduos ───────────────────────────────────────────
  p_painel <- NULL
  if (isTRUE(grafico)) {
    df_res_plot <- data.frame(
      Index = seq_along(res),
      Res   = res
    )

    # 1. Resíduos no Tempo
    p1 <- ggplot2::ggplot(df_res_plot, ggplot2::aes(x = Index, y = Res)) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = c(-2 * stats::sd(res), 2 * stats::sd(res)),
                          linetype = "dashed", color = "#E31A1C", linewidth = 0.5) +
      ggplot2::geom_line(color = "#377EB8", linewidth = 0.75) +
      ggplot2::geom_point(color = "#377EB8", size = 1.3) +
      ggplot2::labs(title = "Res\u00edduos ao Longo do Tempo", x = "\u00cdndice Temporal", y = "Res\u00edduo") +
      tema_estatR(estilo = 2, grade = "dupla")

    # 2. Histograma com Densidade
    p2 <- ggplot2::ggplot(df_res_plot, ggplot2::aes(x = Res)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                              bins = min(20, max(8, floor(n / 6))),
                              fill = "gray85", color = "gray40") +
      ggplot2::stat_function(fun = stats::dnorm,
                             args = list(mean = mean(res), sd = stats::sd(res)),
                             color = "#E31A1C", linewidth = 0.85) +
      ggplot2::labs(title = "Histograma & Densidade Normal", x = "Res\u00edduo", y = "Densidade") +
      tema_estatR(estilo = 2, grade = "dupla")

    # 3. QQ-Plot Normal
    p3 <- ggplot2::ggplot(df_res_plot, ggplot2::aes(sample = Res)) +
      ggplot2::stat_qq(color = "#377EB8", size = 1.5, alpha = 0.8) +
      ggplot2::stat_qq_line(color = "#E31A1C", linewidth = 0.75, linetype = "dashed") +
      ggplot2::labs(title = "QQ-Plot Normal dos Res\u00edduos", x = "Quantis Te\u00f3ricos", y = "Quantis Amostrais") +
      tema_estatR(estilo = 2, grade = "dupla")

    # 4. ACF dos Resíduos (Autocorrelação com limites Bartlett)
    max_lag_acf <- min(24, n - 2)
    acf_obj <- stats::acf(res, lag.max = max_lag_acf, plot = FALSE)
    df_acf_res <- data.frame(
      Lag   = as.numeric(acf_obj$lag)[-1],
      Valor = as.numeric(acf_obj$acf)[-1]
    )
    ic_bart <- 1.96 / sqrt(n)

    p4 <- ggplot2::ggplot(df_acf_res, ggplot2::aes(x = Lag, y = Valor)) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = c(-ic_bart, ic_bart), linetype = "dashed", color = "#E31A1C", linewidth = 0.6) +
      ggplot2::geom_segment(ggplot2::aes(x = Lag, xend = Lag, y = 0, yend = Valor),
                            color = "#377EB8", linewidth = 0.75) +
      ggplot2::geom_point(color = "#377EB8", size = 1.6) +
      ggplot2::labs(title = "FAC dos Res\u00edduos (Ru\u00eddo Branco)", x = "Defasagem (Lag)", y = "Autocorrela\u00e7\u00e3o") +
      tema_estatR(estilo = 2, grade = "dupla")

    p_painel <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)
    print(p_painel)
  }

  invisible(list(
    residuos = res,
    t_test = t_test,
    normalidade = st_norm,
    ljung_box = st_box,
    grafico = p_painel
  ))
}
