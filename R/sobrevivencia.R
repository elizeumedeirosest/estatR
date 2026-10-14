# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: ANÁLISE DE SOBREVIVÊNCIA (BIOESTATÍSTICA E CONFIABILIDADE)
# ─────────────────────────────────────────────────────────────────────────────

#' @keywords internal
.extrair_vetor_sob <- function(expr_sub, dados, arg_nome) {
  if (is.null(dados)) {
    val <- tryCatch(eval(expr_sub, envir = parent.frame(2)), error = function(e) NULL)
    if (is.null(val)) stop(sprintf("O argumento '%s' não foi encontrado.", arg_nome), call. = FALSE)
    return(list(valores = val, nome = deparse(expr_sub)))
  }

  if (is.character(expr_sub) && length(expr_sub) == 1) {
    if (!expr_sub %in% names(dados)) {
      stop(sprintf("A coluna '%s' não foi encontrada em 'dados'.", expr_sub), call. = FALSE)
    }
    return(list(valores = dados[[expr_sub]], nome = expr_sub))
  }

  col_str <- as.character(expr_sub)
  if (length(col_str) == 1 && col_str %in% names(dados)) {
    return(list(valores = dados[[col_str]], nome = col_str))
  }

  val <- tryCatch(eval(expr_sub, envir = dados), error = function(e) NULL)
  if (!is.null(val)) {
    return(list(valores = val, nome = deparse(expr_sub)))
  }

  val <- tryCatch(eval(expr_sub, envir = parent.frame(2)), error = function(e) NULL)
  if (!is.null(val)) {
    return(list(valores = val, nome = deparse(expr_sub)))
  }

  stop(sprintf("Não foi possível extrair a variável '%s'.", arg_nome), call. = FALSE)
}

#' @keywords internal
.normalizar_evento <- function(ev) {
  # Retorna vetor binario (1 = evento/falha, 0 = censura)
  if (is.logical(ev)) {
    return(as.integer(ev))
  }
  if (is.numeric(ev)) {
    u_ev <- sort(unique(stats::na.omit(ev)))
    # Caso classico de datasets como 'lung' do survival: 1 = censura, 2 = morte
    if (length(u_ev) == 2 && all(u_ev == c(1, 2))) {
      return(as.integer(ev == 2))
    }
    # Caso 0 e 1
    if (all(u_ev %in% c(0, 1))) {
      return(as.integer(ev))
    }
    # Caso continuo ou contagens: > 0 eh evento
    return(as.integer(ev > 0))
  }
  # Se for caracter ou fator
  ev_fac <- as.factor(ev)
  levs   <- levels(ev_fac)
  if (length(levs) == 2) {
    # Procura termos comuns de censura
    termo_cens <- grep("cens|vivo|control|nao|no|0", levs, ignore.case = TRUE, value = TRUE)
    if (length(termo_cens) > 0) {
      nivel_evento <- setdiff(levs, termo_cens[1])[1]
      return(as.integer(ev_fac == nivel_evento))
    } else {
      # Segundo nivel por padrao e o evento
      return(as.integer(ev_fac == levs[2]))
    }
  }
  as.integer(ev_fac != levs[1])
}


# ─────────────────────────────────────────────────────────────────────────────
# 1. ESTIMADOR DE KAPLAN-MEIER
# ─────────────────────────────────────────────────────────────────────────────

#' @title Análise de Sobrevivência com Estimador de Kaplan-Meier
#' @description Ajusta e apresenta curvas de sobrevivência não-paramétricas pelo
#'   método de Kaplan-Meier, com cálculo de medianas de sobrevivência, intervalos
#'   de confiança, taxas de censura e teste de Log-Rank para comparação de grupos.
#'
#' @param tempo Variável com os tempos até o evento ou censura (com ou sem aspas).
#' @param evento Variável indicadora do status/evento (1 = evento/óbito, 0 = censura).
#'   Aceita variáveis numéricas (0/1 ou 1/2), lógicas ou fatores.
#' @param grupo Variável categórica opcional que define os grupos de comparação.
#'   Se omitida (\code{NULL}), calcula a curva de sobrevivência global.
#' @param dados Data frame opcional contendo as variáveis.
#' @param confianca Nível de confiança para os intervalos (padrão \code{0.95}).
#' @param tabela_risco Lógico. Se \code{TRUE} (padrão), inclui a tabela de indivíduos
#'   em risco (*Risk Table*) no gráfico.
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o gráfico de sobrevivência.
#'
#' @return Objeto da classe \code{"estatR_sobrevivencia"} contendo o ajuste de
#'   Kaplan-Meier, resumo dos grupos, teste de Log-Rank e gráfico.
#' @export
#' @aliases kaplan_meier km
analise_sobrevivencia <- function(tempo,
                                  evento,
                                  grupo        = NULL,
                                  dados        = NULL,
                                  confianca    = 0.95,
                                  tabela_risco = TRUE,
                                  grafico      = TRUE) {

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("O pacote 'survival' é necessário para a análise de sobrevivência.", call. = FALSE)
  }

  sub_tempo  <- substitute(tempo)
  sub_evento <- substitute(evento)
  sub_grupo  <- substitute(grupo)

  info_tempo  <- .extrair_vetor_sob(sub_tempo,  dados, "tempo")
  info_evento <- .extrair_vetor_sob(sub_evento, dados, "evento")

  t_raw <- as.numeric(info_tempo$valores)
  e_raw <- .normalizar_evento(info_evento$valores)

  tem_grupo <- !is.null(sub_grupo) && !identical(sub_grupo, as.name("NULL"))
  g_raw <- NULL
  nm_grupo <- NULL
  if (tem_grupo) {
    info_grupo <- .extrair_vetor_sob(sub_grupo, dados, "grupo")
    g_raw <- as.factor(info_grupo$valores)
    nm_grupo <- info_grupo$nome
  }

  # Filtro de observacoes completas e tempo > 0
  if (tem_grupo) {
    ok <- !is.na(t_raw) & !is.na(e_raw) & !is.na(g_raw) & (t_raw >= 0)
    t_vec <- t_raw[ok]
    e_vec <- e_raw[ok]
    g_vec <- droplevels(g_raw[ok])
    df_sob <- data.frame(Tempo = t_vec, Evento = e_vec, Grupo = g_vec)
  } else {
    ok <- !is.na(t_raw) & !is.na(e_raw) & (t_raw >= 0)
    t_vec <- t_raw[ok]
    e_vec <- e_raw[ok]
    g_vec <- NULL
    df_sob <- data.frame(Tempo = t_vec, Evento = e_vec)
  }

  n_total   <- length(t_vec)
  n_eventos <- sum(e_vec == 1)
  n_cens    <- sum(e_vec == 0)

  if (n_total < 3) {
    stop("São necessárias ao menos 3 observações válidas para a análise.", call. = FALSE)
  }

  # ── AJUSTE DE KAPLAN-MEIER ────────────────────────────────────────────────
  surv_obj <- survival::Surv(t_vec, e_vec)

  if (tem_grupo && length(levels(g_vec)) > 1) {
    fit_km  <- survival::survfit(surv_obj ~ g_vec, conf.int = confianca)
    # Teste de Log-Rank
    logrank <- survival::survdiff(surv_obj ~ g_vec)
    gl_lr   <- length(logrank$n) - 1
    chi_lr  <- logrank$chisq
    p_lr    <- 1 - stats::pchisq(chi_lr, gl_lr)
    teste_lr <- list(
      estatistica = chi_lr,
      gl          = gl_lr,
      p_valor     = p_lr,
      obs         = logrank$obs,
      exp         = logrank$exp,
      grupos      = levels(g_vec)
    )
  } else {
    fit_km   <- survival::survfit(surv_obj ~ 1, conf.int = confianca)
    teste_lr <- NULL
  }

  # ── TABELA DE RESUMO POR GRUPO ────────────────────────────────────────────
  summ_tab <- summary(fit_km)$table

  if (is.matrix(summ_tab)) {
    # Multiplos grupos
    nomes_grps <- levels(g_vec)
    n_grps     <- as.vector(summ_tab[, "records"])
    ev_grps    <- as.vector(summ_tab[, "events"])
    cens_grps  <- n_grps - ev_grps
    pct_cens   <- (cens_grps / n_grps) * 100
    med_sob    <- as.vector(summ_tab[, "median"])
    ic_inf     <- as.vector(summ_tab[, grep("0.95LCL|LCL|lower", colnames(summ_tab))])
    ic_sup     <- as.vector(summ_tab[, grep("0.95UCL|UCL|upper", colnames(summ_tab))])

    tab_resumo <- data.frame(
      Grupo       = nomes_grps,
      Total       = n_grps,
      Eventos     = ev_grps,
      Censurados  = cens_grps,
      Taxa_Censura= pct_cens,
      Mediana     = med_sob,
      IC_inf      = ic_inf,
      IC_sup      = ic_sup,
      stringsAsFactors = FALSE
    )
  } else {
    # Curva global (1 grupo)
    n_g    <- as.numeric(summ_tab["records"])
    ev_g   <- as.numeric(summ_tab["events"])
    cens_g <- n_g - ev_g
    pct_c  <- (cens_g / n_g) * 100
    med_g  <- as.numeric(summ_tab["median"])
    ic_lo_col <- grep("0.95LCL|LCL|lower", names(summ_tab))
    ic_up_col <- grep("0.95UCL|UCL|upper", names(summ_tab))
    ic_lo  <- if (length(ic_lo_col) > 0) as.numeric(summ_tab[ic_lo_col]) else NA
    ic_up  <- if (length(ic_up_col) > 0) as.numeric(summ_tab[ic_up_col]) else NA

    tab_resumo <- data.frame(
      Grupo       = "Global",
      Total       = n_g,
      Eventos     = ev_g,
      Censurados  = cens_g,
      Taxa_Censura= pct_c,
      Mediana     = med_g,
      IC_inf      = ic_lo,
      IC_sup      = ic_up,
      stringsAsFactors = FALSE
    )
  }

  resultado <- structure(
    list(
      nomes = list(
        tempo  = info_tempo$nome,
        evento = info_evento$nome,
        grupo  = nm_grupo
      ),
      tem_grupo    = tem_grupo,
      n_total      = n_total,
      n_eventos    = n_eventos,
      n_censurados = n_cens,
      confianca    = confianca,
      tabela       = tab_resumo,
      teste_logrank= teste_lr,
      fit_km       = fit_km,
      dados        = df_sob,
      tabela_risco = tabela_risco,
      grafico      = NULL
    ),
    class = "estatR_sobrevivencia"
  )

  # Gera grafico
  g_km <- tryCatch(
    grafico_sobrevivencia(resultado, tabela_risco = tabela_risco),
    error = function(e) NULL
  )
  resultado$grafico <- g_km

  if (isTRUE(grafico) && !is.null(g_km)) {
    print(g_km)
  }

  resultado
}

#' @export
kaplan_meier <- analise_sobrevivencia

#' @export
km <- analise_sobrevivencia


# ─────────────────────────────────────────────────────────────────────────────
# MÉTODO PRINT KAPLAN-MEIER
# ─────────────────────────────────────────────────────────────────────────────

#' @export
print.estatR_sobrevivencia <- function(x, ...) {

  tit <- if (x$tem_grupo) {
    sprintf("ANÁLISE DE SOBREVIVÊNCIA DE KAPLAN-MEIER (POR %s)", toupper(x$nomes$grupo))
  } else {
    "ANÁLISE DE SOBREVIVÊNCIA DE KAPLAN-MEIER (CURVA GLOBAL)"
  }
  .print_titulo(tit)

  cat(sprintf("  Variável de Tempo:   %s\n", x$nomes$tempo))
  cat(sprintf("  Variável de Evento:  %s\n", x$nomes$evento))
  if (x$tem_grupo) {
    cat(sprintf("  Fator de Grupo:      %s (%d categorias)\n",
                x$nomes$grupo, nrow(x$tabela)))
  }
  cat(sprintf("  Total de Casos:      %d  (Eventos: %d | Censuras: %d - %.1f%%)\n\n",
              x$n_total, x$n_eventos, x$n_censurados,
              (x$n_censurados / x$n_total) * 100))

  # 1. TABELA RESUMO DE SOBREVIVÊNCIA
  .print_topico(sprintf("RESUMO DE SOBREVIVÊNCIA E TEMPO MEDIANO (IC %g%%)", x$confianca * 100))
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  Grupo             Total    Eventos   Censuras (%)    Mediana       IC (95%)    \n")
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")

  for (i in 1:nrow(x$tabela)) {
    grp_s  <- x$tabela$Grupo[i]
    tot_s  <- sprintf("%5d", x$tabela$Total[i])
    ev_s   <- sprintf("%5d", x$tabela$Eventos[i])
    cens_s <- sprintf("%4d (%4.1f%%)", x$tabela$Censurados[i], x$tabela$Taxa_Censura[i])

    med_val <- x$tabela$Mediana[i]
    med_s   <- if (!is.na(med_val)) sprintf("%8s", .fmt_num(med_val, 1)) else "     NE*"

    lo_val <- x$tabela$IC_inf[i]
    up_val <- x$tabela$IC_sup[i]
    if (!is.na(lo_val) && !is.na(up_val)) {
      ic_s <- sprintf("[%s ; %s]", .fmt_num(lo_val, 1), .fmt_num(up_val, 1))
    } else if (!is.na(lo_val)) {
      ic_s <- sprintf("[%s ;   NE]", .fmt_num(lo_val, 1))
    } else {
      ic_s <- "   [ NE ; NE ]"
    }

    cat(sprintf("  %-16s %s    %s   %12s   %8s   %15s\n",
                grp_s, tot_s, ev_s, cens_s, med_s, ic_s))
  }
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  * NE: Não Estimado (a probabilidade de sobrevivência não atingiu 50% no período).\n\n")

  # 2. TESTE DE LOG-RANK (SE HOUVER GRUPOS)
  if (!is.null(x$teste_logrank)) {
    .print_topico("COMPARAÇÃO ENTRE GRUPOS (TESTE DE LOG-RANK)")
    lr <- x$teste_logrank
    pv_str <- if (lr$p_valor < 0.001) "< 0,001 ***" else sprintf("%.3f", lr$p_valor)

    cat("  ────────────────────────────────────────────────────────────── \n")
    cat("  Grupo               Observado      Esperado    (O - E)²/E\n")
    cat("  ────────────────────────────────────────────────────────────── \n")
    for (i in seq_along(lr$grupos)) {
      o_val <- lr$obs[i]
      e_val <- lr$exp[i]
      v_val <- ((o_val - e_val)^2) / e_val
      cat(sprintf("  %-18s  %7d     %9.1f       %8.2f\n",
                  lr$grupos[i], o_val, e_val, v_val))
    }
    cat("  ────────────────────────────────────────────────────────────── \n")
    cat(sprintf("  Qui-quadrado: %s  |  GL: %d  |  p-valor: %s\n",
                .fmt_num(lr$estatistica, 2), lr$gl, pv_str))
    cat("  ────────────────────────────────────────────────────────────── \n\n")

    # 3. INTERPRETAÇÃO DO LOG-RANK
    .print_topico("INTERPRETAÇÃO ESTATÍSTICA")
    if (lr$p_valor < 0.05) {
      cat(sprintf("  O teste de Log-Rank indicou diferença estatisticamente significativa entre\n  as curvas de sobrevivência dos grupos de %s (Chi2(%d) = %s; p %s).\n\n",
                  x$nomes$grupo, lr$gl, .fmt_num(lr$estatistica, 2),
                  if (lr$p_valor < 0.001) "< 0,001" else paste0("= ", .fmt_num(lr$p_valor, 3))))

      # Identifica o grupo de melhor prognostico
      meds <- stats::setNames(x$tabela$Mediana, x$tabela$Grupo)
      meds_validas <- meds[!is.na(meds)]
      if (length(meds_validas) > 0) {
        melhor <- names(which.max(meds_validas))
        cat(sprintf("  O grupo '%s' apresentou o maior tempo mediano de sobrevivência observado\n  (%s unidades de tempo).\n",
                    melhor, .fmt_num(max(meds_validas), 1)))
      }
    } else {
      cat(sprintf("  O teste de Log-Rank não encontrou evidências de diferença estatística entre\n  as curvas de sobrevivência (Chi2(%d) = %s; p = %s),\n  sugerindo prognósticos semelhantes entre os grupos avaliados.\n",
                  lr$gl, .fmt_num(lr$estatistica, 2), .fmt_num(lr$p_valor, 3)))
    }
    cat("\n")
  } else {
    .print_topico("INTERPRETAÇÃO")
    med_g <- x$tabela$Mediana[1]
    if (!is.na(med_g)) {
      cat(sprintf("  O tempo mediano estimado de sobrevivência foi de %s unidades de tempo,\n  representando o ponto em que 50%% dos indivíduos apresentaram o evento de interesse.\n\n",
                  .fmt_num(med_g, 1)))
    } else {
      cat("  O tempo mediano não pôde ser alcançado porque mais de 50% dos indivíduos\n  permaneceram censurados até o término do estudo.\n\n")
    }
  }

  .print_rodape()
  invisible(x)
}


# ─────────────────────────────────────────────────────────────────────────────
# GRÁFICO DE SOBREVIVÊNCIA
# ─────────────────────────────────────────────────────────────────────────────

#' @title Gráfico de Curva de Sobrevivência de Kaplan-Meier
#' @description Gera gráfico da curva de sobrevivência em degraus com marcas de
#'   censura, intervalo de confiança sombreado e tabela de indivíduos em risco (*Risk Table*).
#' @param modelo Objeto retornado por \code{analise_sobrevivencia()} ou \code{kaplan_meier()}.
#' @param tabela_risco Lógico. Se \code{TRUE} (padrão), inclui a tabela de indivíduos em risco.
#' @param conf_int Lógico. Se \code{TRUE} (padrão), inclui a faixa de confiança sombreada.
#' @param ... Argumentos adicionais (não utilizados).
#' @return Objeto \code{ggplot} ou composição de gráficos.
#' @export
#' @aliases plot_sobrevivencia
grafico_sobrevivencia <- function(modelo, tabela_risco = TRUE, conf_int = TRUE, ...) {
  if (!inherits(modelo, "estatR_sobrevivencia")) {
    stop("O argumento 'modelo' deve ser um objeto retornado por analise_sobrevivencia().", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("O pacote 'ggplot2' é necessário para gerar gráficos.")
    return(invisible(NULL))
  }

  fit <- modelo$fit_km
  tem_grupo <- modelo$tem_grupo

  # Extrai dados da curva do fit
  sf <- summary(fit, censored = TRUE)

  if (tem_grupo) {
    # Nomes dos estratos limpos
    raw_strata <- as.character(sf$strata)
    strata_clean <- gsub("^.*=", "", raw_strata)
    df_plot <- data.frame(
      time   = sf$time,
      surv   = sf$surv,
      lower  = sf$lower,
      upper  = sf$upper,
      n.risk = sf$n.risk,
      n.event= sf$n.event,
      n.censor = sf$n.censor,
      grupo  = strata_clean,
      stringsAsFactors = FALSE
    )
    # Adiciona ponto inicial t=0 para cada grupo
    grps <- unique(df_plot$grupo)
    p0_list <- lapply(grps, function(g) {
      data.frame(time = 0, surv = 1, lower = 1, upper = 1,
                 n.risk = max(df_plot$n.risk[df_plot$grupo == g]),
                 n.event = 0, n.censor = 0, grupo = g, stringsAsFactors = FALSE)
    })
    df_plot <- rbind(do.call(rbind, p0_list), df_plot)
    df_plot <- df_plot[order(df_plot$grupo, df_plot$time), ]
  } else {
    df_plot <- data.frame(
      time   = c(0, sf$time),
      surv   = c(1, sf$surv),
      lower  = c(1, sf$lower),
      upper  = c(1, sf$upper),
      n.risk = c(max(sf$n.risk), sf$n.risk),
      n.event= c(0, sf$n.event),
      n.censor = c(0, sf$n.censor),
      grupo  = "Global",
      stringsAsFactors = FALSE
    )
  }

  # Paleta e Tema
  pal <- if (exists(".paleta_estatR", mode = "function")) {
    .paleta_estatR()
  } else {
    c("#2166AC", "#D6604D", "#4DAC26", "#B2182B", "#762A83", "#E08214")
  }
  n_grps <- length(unique(df_plot$grupo))
  cores  <- pal[seq_len(n_grps)]

  tema_fn <- if (exists("tema_estatR", mode = "function")) {
    tema_estatR(estilo = 2, grade = "dupla")
  } else if (exists("meu_tema", mode = "function")) {
    meu_tema(estilo = 2, grade = "dupla")
  } else {
    ggplot2::theme_minimal()
  }

  # Subtitulo com teste de log-rank se disponivel
  sub_txt <- if (!is.null(modelo$teste_logrank)) {
    lr <- modelo$teste_logrank
    pv <- if (lr$p_valor < 0.001) "p < 0,001" else sprintf("p = %s", .fmt_num(lr$p_valor, 3))
    sprintf("Teste de Log-Rank: Chi2(%d) = %s, %s", lr$gl, .fmt_num(lr$estatistica, 2), pv)
  } else {
    sprintf("Total de indivíduos: %d | Eventos: %d | Censuras: %d",
            modelo$n_total, modelo$n_eventos, modelo$n_censurados)
  }

  # Gráfico principal da curva
  if (tem_grupo) {
    p_curva <- ggplot2::ggplot(df_plot, ggplot2::aes(x = time, y = surv, color = grupo, fill = grupo))
  } else {
    p_curva <- ggplot2::ggplot(df_plot, ggplot2::aes(x = time, y = surv))
  }

  # Faixa de confianca
  if (isTRUE(conf_int)) {
    if (tem_grupo) {
      p_curva <- p_curva +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                             alpha = 0.15, linetype = 0)
    } else {
      p_curva <- p_curva +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                             fill = "#2166AC", alpha = 0.18, linetype = 0)
    }
  }

  # Degrau da curva
  if (tem_grupo) {
    p_curva <- p_curva +
      ggplot2::geom_step(linewidth = 0.9) +
      ggplot2::scale_color_manual(values = cores) +
      ggplot2::scale_fill_manual(values = cores)
  } else {
    p_curva <- p_curva +
      ggplot2::geom_step(color = "#2166AC", linewidth = 0.95)
  }

  # Marcas de censura (+)
  df_cens <- df_plot[df_plot$n.censor > 0, ]
  if (nrow(df_cens) > 0) {
    if (tem_grupo) {
      p_curva <- p_curva +
        ggplot2::geom_point(data = df_cens, ggplot2::aes(x = time, y = surv, color = grupo),
                            shape = 3, size = 2.5, stroke = 1.1, show.legend = FALSE)
    } else {
      p_curva <- p_curva +
        ggplot2::geom_point(data = df_cens, ggplot2::aes(x = time, y = surv),
                            color = "#2166AC", shape = 3, size = 2.5, stroke = 1.1)
    }
  }

  p_curva <- p_curva +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                                limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title    = "Curva de Sobrevivência de Kaplan-Meier",
      subtitle = sub_txt,
      x        = sprintf("Tempo (%s)", modelo$nomes$tempo),
      y        = "Probabilidade de Sobrevivência",
      color    = modelo$nomes$grupo,
      fill     = modelo$nomes$grupo
    ) +
    tema_fn

  if (!isTRUE(tabela_risco) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(p_curva)
  }

  # ── TABELA DE RISCO (NUMBER AT RISK) ──────────────────────────────────────
  max_t <- max(df_plot$time, na.rm = TRUE)
  tempos_quebras <- pretty(c(0, max_t), n = 5)
  tempos_quebras <- tempos_quebras[tempos_quebras <= max_t]

  # Calcula individuos em risco em cada tempo de quebra
  grps_lista <- if (tem_grupo) unique(df_plot$grupo) else "Global"
  linhas_risco <- list()

  for (g in grps_lista) {
    sub_fit <- if (tem_grupo) df_plot[df_plot$grupo == g, ] else df_plot
    for (t_pt in tempos_quebras) {
      # Pega o ultimo n.risk registrado em tempo <= t_pt
      sub_t <- sub_fit[sub_fit$time <= t_pt, ]
      n_r <- if (nrow(sub_t) > 0) sub_t$n.risk[nrow(sub_t)] - sub_t$n.event[nrow(sub_t)] - sub_t$n.censor[nrow(sub_t)] else max(sub_fit$n.risk)
      if (t_pt == 0) n_r <- max(sub_fit$n.risk)
      linhas_risco[[length(linhas_risco) + 1]] <- data.frame(
        tempo = t_pt,
        grupo = g,
        em_risco = max(0, n_r),
        stringsAsFactors = FALSE
      )
    }
  }
  df_risco <- do.call(rbind, linhas_risco)

  p_risco <- ggplot2::ggplot(df_risco, ggplot2::aes(x = tempo, y = grupo, label = em_risco)) +
    ggplot2::geom_text(size = 3.3, fontface = "bold", color = "gray20") +
    ggplot2::scale_x_continuous(limits = ggplot2::layer_scales(p_curva)$x$range$range,
                                breaks = tempos_quebras) +
    ggplot2::labs(
      x = sprintf("Tempo (%s)", modelo$nomes$tempo),
      y = ""
    ) +
    tema_fn +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(face = "bold", color = "gray20"),
      axis.ticks.y     = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(t = 0, r = 10, b = 10, l = 10)
    )

  # Remove eixo X do grafico superior para composicao perfeita
  p_curva_sem_x <- p_curva +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(t = 10, r = 10, b = 0, l = 10)
    )

  painel_final <- patchwork::wrap_plots(p_curva_sem_x, p_risco, ncol = 1, heights = c(3.8, 1))
  painel_final
}

#' @export
plot_sobrevivencia <- grafico_sobrevivencia


# ─────────────────────────────────────────────────────────────────────────────
# 2. MODELO DE RISCOS PROPORCIONAIS DE COX
# ─────────────────────────────────────────────────────────────────────────────

#' @title Regressão de Cox (Modelo de Riscos Proporcionais)
#' @description Ajusta o modelo semi-paramétrico de Cox para análise de tempo até
#'   o evento, fornecendo Hazard Ratios (HR), intervalos de confiança de 95%,
#'   estatísticas globais de qualidade de ajuste e diagnóstico formal de proporcionalidade
#'   dos riscos baseado nos resíduos de Schoenfeld.
#'
#' @param formula Fórmula do modelo. Aceita \code{Surv(tempo, evento) ~ x1 + x2}
#'   ou sintaxe direta com variáveis contidas em \code{dados}.
#' @param dados Data frame contendo as variáveis do modelo.
#' @param confianca Nível de confiança para os Hazard Ratios (padrão \code{0.95}).
#' @param diagnostico Lógico. Se \code{TRUE} (padrão), calcula o teste de Schoenfeld
#'   para checar a proporcionalidade dos riscos de cada covariável.
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o Forest Plot dos Hazard Ratios.
#'
#' @return Objeto da classe \code{"estatR_cox"} contendo estimativas, Hazard Ratios,
#'   diagnóstico de resíduos e modelo ajustado.
#' @export
#' @aliases cox
regressao_cox <- function(formula,
                          dados,
                          confianca   = 0.95,
                          diagnostico = TRUE,
                          grafico     = TRUE) {

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("O pacote 'survival' é necessário para a regressão de Cox.", call. = FALSE)
  }

  if (missing(dados) || is.null(dados)) {
    stop("O argumento 'dados' (data frame) deve ser informado.", call. = FALSE)
  }

  # Ajuste via survival::coxph
  mod_cox <- tryCatch(
    survival::coxph(formula, data = dados, model = TRUE, x = TRUE, y = TRUE),
    error = function(e) {
      stop(sprintf("Erro ao ajustar o modelo de Cox: %s", e$message), call. = FALSE)
    }
  )

  sum_cox <- summary(mod_cox)
  alpha   <- 1 - confianca
  z_crit  <- stats::qnorm(1 - alpha / 2)

  # Tabela de coeficientes e Hazard Ratios
  coef_mat <- sum_cox$coefficients
  conf_mat <- sum_cox$conf.int

  covariaveis <- rownames(coef_mat)
  beta_vals   <- coef_mat[, "coef"]
  ep_vals     <- coef_mat[, "se(coef)"]
  z_vals      <- coef_mat[, "z"]
  p_vals      <- coef_mat[, "Pr(>|z|)"]

  hr_vals     <- exp(beta_vals)
  hr_lo       <- exp(beta_vals - z_crit * ep_vals)
  hr_hi       <- exp(beta_vals + z_crit * ep_vals)

  tab_cox <- data.frame(
    Covariavel  = covariaveis,
    Beta        = beta_vals,
    EP          = ep_vals,
    HR          = hr_vals,
    IC_inf      = hr_lo,
    IC_sup      = hr_hi,
    Z_calc      = z_vals,
    p_valor     = p_vals,
    stringsAsFactors = FALSE
  )

  # Medidas globais
  c_index <- sum_cox$concordance[1]
  se_c    <- sum_cox$concordance[2]

  lrt_stat <- sum_cox$logtest["test"]
  lrt_df   <- sum_cox$logtest["df"]
  lrt_p    <- sum_cox$logtest["pvalue"]

  wald_stat <- sum_cox$waldtest["test"]
  wald_df   <- sum_cox$waldtest["df"]
  wald_p    <- sum_cox$waldtest["pvalue"]

  score_stat <- sum_cox$sctest["test"]
  score_df   <- sum_cox$sctest["df"]
  score_p    <- sum_cox$sctest["pvalue"]

  aic_val <- stats::AIC(mod_cox)
  bic_val <- stats::BIC(mod_cox)

  # Diagnostico de Riscos Proporcionais (Schoenfeld)
  diag_zph <- NULL
  if (isTRUE(diagnostico)) {
    diag_zph <- tryCatch(survival::cox.zph(mod_cox), error = function(e) NULL)
  }

  resultado <- structure(
    list(
      tabela      = tab_cox,
      n_total     = mod_cox$n,
      n_eventos   = mod_cox$nevent,
      confianca   = confianca,
      c_index     = list(valor = c_index, ep = se_c),
      lrt         = list(stat = lrt_stat, df = lrt_df, p = lrt_p),
      wald        = list(stat = wald_stat, df = wald_df, p = wald_p),
      score       = list(stat = score_stat, df = score_df, p = score_p),
      aic         = aic_val,
      bic         = bic_val,
      diagnostico = diag_zph,
      modelo_cox  = mod_cox,
      dados       = dados,
      formula     = formula,
      grafico     = NULL
    ),
    class = "estatR_cox"
  )

  # Forest plot dos Hazard Ratios
  g_forest <- tryCatch(grafico_cox(resultado), error = function(e) NULL)
  resultado$grafico <- g_forest

  if (isTRUE(grafico) && !is.null(g_forest)) {
    print(g_forest)
  }

  resultado
}

#' @export
cox <- regressao_cox


# ─────────────────────────────────────────────────────────────────────────────
# MÉTODO PRINT COX
# ─────────────────────────────────────────────────────────────────────────────

#' @export
print.estatR_cox <- function(x, ...) {

  .print_titulo("MODELO DE RISCOS PROPORCIONAIS DE COX")

  cat(sprintf("  Fórmula do Modelo:   %s\n", deparse(x$formula)))
  cat(sprintf("  Total de Casos:      %d  (Eventos: %d | Censuras: %d - %.1f%%)\n\n",
              x$n_total, x$n_eventos, x$n_total - x$n_eventos,
              ((x$n_total - x$n_eventos) / x$n_total) * 100))

  # 1. TABELA DE COEFICIENTES E HAZARD RATIOS
  .print_topico(sprintf("ESTIMATIVAS E HAZARD RATIOS (HR — IC %g%%)", x$confianca * 100))
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  Covariável              Beta       EP       HR      IC (95%)         z      p-valor \n")
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")

  for (i in 1:nrow(x$tabela)) {
    cov_s <- x$tabela$Covariavel[i]
    b_s   <- sprintf("%+8s", .fmt_num(x$tabela$Beta[i], 3))
    ep_s  <- sprintf("%6s",  .fmt_num(x$tabela$EP[i], 3))
    hr_s  <- sprintf("%6s",  .fmt_num(x$tabela$HR[i], 3))
    ic_s  <- sprintf("[%s ; %s]", .fmt_num(x$tabela$IC_inf[i], 2), .fmt_num(x$tabela$IC_sup[i], 2))
    z_s   <- sprintf("%+6s", .fmt_num(x$tabela$Z_calc[i], 2))

    pv    <- x$tabela$p_valor[i]
    pv_s  <- if (pv < 0.001) "< 0,001 ***"
             else if (pv < 0.01)  sprintf("%s ** ", .fmt_num(pv, 3))
             else if (pv < 0.05)  sprintf("%s *  ", .fmt_num(pv, 3))
             else if (pv < 0.10)  sprintf("%s .  ", .fmt_num(pv, 3))
             else sprintf("%s    ", .fmt_num(pv, 3))

    cat(sprintf("  %-20s %s   %s   %s   %-15s %s  %s\n",
                cov_s, b_s, ep_s, hr_s, ic_s, z_s, pv_s))
  }
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  * Significância: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")

  # 2. QUALIDADE DE AJUSTE GLOBAL
  .print_topico("AVALIAÇÃO GLOBAL DO MODELO")
  cat("  ────────────────────────────────────────────────────────────── \n")
  cat(sprintf("  Índice C (Concordância):     %s  (EP: %s)\n",
              .fmt_num(x$c_index$valor, 3), .fmt_num(x$c_index$ep, 3)))
  cat(sprintf("  Razão de Verossimilhança:    Chi2(%d) = %s  (p %s)\n",
              x$lrt$df, .fmt_num(x$lrt$stat, 2),
              if (x$lrt$p < 0.001) "< 0,001" else paste0("= ", .fmt_num(x$lrt$p, 3))))
  cat(sprintf("  Teste de Wald:               Chi2(%d) = %s  (p %s)\n",
              x$wald$df, .fmt_num(x$wald$stat, 2),
              if (x$wald$p < 0.001) "< 0,001" else paste0("= ", .fmt_num(x$wald$p, 3))))
  cat(sprintf("  Critérios de Informação:     AIC = %s  |  BIC = %s\n",
              .fmt_num(x$aic, 1), .fmt_num(x$bic, 1)))
  cat("  * Concordância (C-index): 0.50 (aleatório), > 0.70 (bom), > 0.80 (excelente discriminação).\n\n")

  # 3. INTERPRETAÇÃO DOS HAZARD RATIOS
  .print_topico("INTERPRETAÇÃO DOS HAZARD RATIOS")
  for (i in 1:nrow(x$tabela)) {
    cov_i <- x$tabela$Covariavel[i]
    hr_i  <- x$tabela$HR[i]
    pv_i  <- x$tabela$p_valor[i]
    if (pv_i < 0.05) {
      if (hr_i > 1) {
        pct_aumento <- (hr_i - 1) * 100
        cat(sprintf("  A covariável '%s' aumenta significativamente o risco instantâneo do evento\n  em %.1f%% (HR = %s; p = %s).\n\n",
                    cov_i, pct_aumento, .fmt_num(hr_i, 2), .fmt_num(pv_i, 3)))
      } else {
        pct_reducao <- (1 - hr_i) * 100
        cat(sprintf("  A covariável '%s' atua como fator de proteção, reduzindo o risco instantâneo\n  do evento em %.1f%% (HR = %s; p = %s).\n\n",
                    cov_i, pct_reducao, .fmt_num(hr_i, 2), .fmt_num(pv_i, 3)))
      }
    } else {
      cat(sprintf("  A covariável '%s' não apresentou efeito estatisticamente significativo sobre\n  o risco do evento (HR = %s; p = %s).\n\n",
                  cov_i, .fmt_num(hr_i, 2), .fmt_num(pv_i, 3)))
    }
  }

  .print_rodape()
  invisible(x)
}


# ─────────────────────────────────────────────────────────────────────────────
# FOREST PLOT DO MODELO DE COX
# ─────────────────────────────────────────────────────────────────────────────

#' @title Forest Plot dos Hazard Ratios do Modelo de Cox
#' @description Gera gráfico de Forest Plot com os Hazard Ratios pontuais e
#'   intervalos de confiança de 95%, com linha de referência em HR = 1.
#' @param modelo Objeto retornado por \code{regressao_cox()} ou \code{cox()}.
#' @param escala_log Lógico. Se \code{TRUE} (padrão), utiliza escala logarítmica
#'   no eixo dos Hazard Ratios.
#' @param ... Argumentos adicionais (não utilizados).
#' @return Objeto \code{ggplot}.
#' @export
#' @aliases plot_cox
grafico_cox <- function(modelo, escala_log = TRUE, ...) {
  if (!inherits(modelo, "estatR_cox")) {
    stop("O argumento 'modelo' deve ser um objeto retornado por regressao_cox().", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("O pacote 'ggplot2' é necessário para gerar gráficos.")
    return(invisible(NULL))
  }

  df_f <- modelo$tabela
  # Reordena para exibir de cima para baixo
  df_f$Covariavel <- factor(df_f$Covariavel, levels = rev(df_f$Covariavel))
  df_f$Significativo <- ifelse(df_f$p_valor < 0.05, "Significativo (p < 0.05)", "Não significativo")

  cores <- c(
    "Significativo (p < 0.05)" = "#2166AC",
    "Não significativo"         = "gray50"
  )

  tema_fn <- if (exists("tema_estatR", mode = "function")) {
    tema_estatR(estilo = 2, grade = "dupla")
  } else if (exists("meu_tema", mode = "function")) {
    meu_tema(estilo = 2, grade = "dupla")
  } else {
    ggplot2::theme_minimal()
  }

  p <- ggplot2::ggplot(df_f, ggplot2::aes(x = HR, y = Covariavel, color = Significativo)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = IC_inf, xmax = IC_sup), width = 0.25, linewidth = 0.85) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_color_manual(values = cores) +
    ggplot2::labs(
      title    = "Forest Plot — Hazard Ratios do Modelo de Cox",
      subtitle = sprintf("Concordância: %.2f | Eventos: %d/%d | Linha vermelha: HR = 1 (sem efeito)",
                         modelo$c_index$valor, modelo$n_eventos, modelo$n_total),
      x        = "Hazard Ratio (HR) e Intervalo de Confiança 95%",
      y        = "Covariável",
      color    = ""
    ) +
    tema_fn +
    ggplot2::theme(legend.position = "bottom")

  if (isTRUE(escala_log)) {
    p <- p + ggplot2::scale_x_log10()
  }

  p
}

#' @export
plot_cox <- grafico_cox


# ─────────────────────────────────────────────────────────────────────────────
# 3. ANÁLISE RESIDUAL DO MODELO DE COX
# ─────────────────────────────────────────────────────────────────────────────

#' @title Análise de Resíduos do Modelo de Cox
#' @description Executa diagnóstico completo de resíduos para modelos de riscos
#'   proporcionais de Cox, incluindo avaliação de proporcionalidade (Schoenfeld),
#'   detecção de outliers (Deviance) e forma funcional/linearidade (Martingale).
#'
#' @param modelo Objeto do modelo de Cox (\code{"estatR_cox"} ou \code{"coxph"}).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o painel de diagnóstico gráfico.
#' @return Lista invisível contendo os testes e dados de resíduos.
#' @export
analise_residual_cox <- function(modelo, grafico = TRUE) {

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("O pacote 'survival' é necessário para a análise de resíduos de Cox.", call. = FALSE)
  }

  mod_cox <- if (inherits(modelo, "estatR_cox")) modelo$modelo_cox else modelo
  if (!inherits(mod_cox, "coxph")) {
    stop("O argumento 'modelo' deve ser um modelo de Cox ajustado.", call. = FALSE)
  }

  .print_titulo("DIAGNÓSTICO DE RESÍDUOS — MODELO DE COX")

  # 1. Teste de Schoenfeld (Riscos Proporcionais)
  zph <- if (inherits(modelo, "estatR_cox") && !is.null(modelo$diagnostico)) {
    modelo$diagnostico
  } else {
    survival::cox.zph(mod_cox)
  }
  tab_zph <- zph$table

  .print_topico("1. PRESSUPOSTO DE RISCOS PROPORCIONAIS (SCHOENFELD)")
  cat("  ─────────────────────────────────────────────────────── \n")
  cat("  Covariável              rho      Chi-quadrado   p-valor \n")
  cat("  ─────────────────────────────────────────────────────── \n")
  for (i in 1:nrow(tab_zph)) {
    nm_r <- rownames(tab_zph)[i]
    rho_v <- if (!is.na(tab_zph[i, "chisq"]) && "rho" %in% colnames(tab_zph)) sprintf("%+7s", .fmt_num(tab_zph[i, "rho"], 3)) else "      -"
    chi_v <- sprintf("%8s", .fmt_num(tab_zph[i, "chisq"], 2))
    pv_v  <- tab_zph[i, "p"]
    pv_s  <- if (pv_v < 0.001) "< 0,001" else sprintf("%7s", .fmt_num(pv_v, 3))
    status_sig <- if (pv_v < 0.05) " [Violação]" else " [Atendido]"
    cat(sprintf("  %-20s %s      %s    %s %s\n",
                nm_r, rho_v, chi_v, pv_s, status_sig))
  }
  cat("  ─────────────────────────────────────────────────────── \n")
  cat("  * Critério: p >= 0.05 indica que a premissa de proporcionalidade dos riscos é válida.\n\n")

  # 2. Residuos Deviance (Outliers)
  res_dev <- stats::residuals(mod_cox, type = "deviance")
  n_outliers <- sum(abs(res_dev) > 2.5)

  .print_topico("2. DETECÇÃO DE CASOS ATÍPICOS (RESÍDUOS DEVIANCE)")
  cat(sprintf("  Média dos resíduos Deviance: %s  (DP: %s)\n",
              .fmt_num(mean(res_dev), 3), .fmt_num(stats::sd(res_dev), 3)))
  cat(sprintf("  Casos atípicos (|Deviance| > 2,5): %d de %d (%.1f%%)\n",
              n_outliers, length(res_dev), (n_outliers / length(res_dev)) * 100))
  if (n_outliers > 0) {
    idx_out <- which(abs(res_dev) > 2.5)
    cat(sprintf("  Índice dos casos atípicos: %s\n", paste(head(idx_out, 10), collapse = ", ")))
  }
  cat("\n")

  # 3. Residuos Martingale
  res_mart <- stats::residuals(mod_cox, type = "martingale")

  .print_topico("INTERPRETAÇÃO DO DIAGNÓSTICO")
  p_global <- tab_zph["GLOBAL", "p"]
  if (!is.na(p_global) && p_global >= 0.05) {
    cat("  O teste global de Schoenfeld confirmou o pressuposto de Riscos Proporcionais\n  (p >= 0,05), validando o uso do modelo de Cox para o conjunto de dados.\n\n")
  } else {
    cat("  Atenção: O teste global de Schoenfeld indicou violação do pressuposto de\n  proporcionalidade para uma ou mais covariáveis (p < 0,05).\n  Considere estratificação ou inclusão de termos dependentes do tempo.\n\n")
  }

  .print_rodape()

  # ── PAINEL GRÁFICO ────────────────────────────────────────────────────────
  p_painel <- NULL
  if (isTRUE(grafico) && requireNamespace("ggplot2", quietly = TRUE)) {
    tema_fn <- if (exists("tema_estatR", mode = "function")) {
      tema_estatR(estilo = 2, grade = "dupla")
    } else {
      ggplot2::theme_minimal()
    }

    # Gráfico 1: Deviance Residuals vs Índice
    df_dev <- data.frame(
      Indice = seq_along(res_dev),
      Deviance = res_dev,
      Outlier = ifelse(abs(res_dev) > 2.5, "Outlier (|d| > 2.5)", "Normal")
    )
    p1 <- ggplot2::ggplot(df_dev, ggplot2::aes(x = Indice, y = Deviance, color = Outlier)) +
      ggplot2::geom_point(alpha = 0.7, size = 2) +
      ggplot2::geom_hline(yintercept = c(-2.5, 0, 2.5),
                          linetype = c("dashed", "solid", "dashed"),
                          color = c("red", "gray40", "red")) +
      ggplot2::scale_color_manual(values = c("Normal" = "#2166AC", "Outlier (|d| > 2.5)" = "red")) +
      ggplot2::labs(title = "Resíduos Deviance (Outliers)", x = "Índice do Indivíduo", y = "Resíduo Deviance") +
      tema_fn +
      ggplot2::theme(legend.position = "none")

    # Gráfico 2: Martingale Residuals vs Preditor Linear
    lp_vals <- mod_cox$linear.predictors
    df_mart <- data.frame(LP = lp_vals, Martingale = res_mart)
    p2 <- ggplot2::ggplot(df_mart, ggplot2::aes(x = LP, y = Martingale)) +
      ggplot2::geom_point(color = "#4DAC26", alpha = 0.65, size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::geom_smooth(method = "loess", color = "black", linewidth = 0.8, se = FALSE) +
      ggplot2::labs(title = "Resíduos Martingale (Linearidade)", x = "Preditor Linear", y = "Resíduo Martingale") +
      tema_fn

    if (requireNamespace("patchwork", quietly = TRUE)) {
      p_painel <- patchwork::wrap_plots(p1, p2, ncol = 2)
      print(p_painel)
    } else {
      print(p1)
      p_painel <- p1
    }
  }

  invisible(list(
    schoenfeld = zph,
    deviance   = res_dev,
    martingale = res_mart,
    grafico    = p_painel
  ))
}


# ─────────────────────────────────────────────────────────────────────────────
# 4. MÉTRICAS DO MODELO DE COX
# ─────────────────────────────────────────────────────────────────────────────

#' @title Métricas de Avaliação do Modelo de Cox
#' @description Exibe e retorna as principais medidas de qualidade de ajuste,
#'   discriminação e critérios de informação para o modelo de riscos proporcionais de Cox.
#' @param modelo Objeto retornado por \code{regressao_cox()} ou \code{coxph}.
#' @return Data frame invisível contendo as métricas calculadas.
#' @export
metricas_cox <- function(modelo) {

  mod_cox <- if (inherits(modelo, "estatR_cox")) modelo$modelo_cox else modelo
  if (!inherits(mod_cox, "coxph")) {
    stop("O argumento 'modelo' deve ser um modelo de Cox (coxph ou estatR_cox).", call. = FALSE)
  }

  sum_c <- summary(mod_cox)
  n_obs <- mod_cox$n
  n_ev  <- mod_cox$nevent

  # Concordancia
  c_idx <- sum_c$concordance[1]
  se_c  <- sum_c$concordance[2]

  # Log-Verossimilhanca e Pseudo-R2
  loglik_fit  <- mod_cox$loglik[2]
  loglik_null <- mod_cox$loglik[1]
  lrt_stat    <- 2 * (loglik_fit - loglik_null)

  # Cox & Snell R2 e Nagelkerke R2
  r2_cs  <- 1 - exp(-lrt_stat / n_obs)
  r2_max <- 1 - exp(2 * loglik_null / n_obs)
  r2_nag <- if (r2_max > 0) r2_cs / r2_max else NA

  # Criterios de Informacao
  aic_val <- stats::AIC(mod_cox)
  bic_val <- stats::BIC(mod_cox)

  # Testes de significancia global
  lrt_p   <- sum_c$logtest["pvalue"]
  wald_p  <- sum_c$waldtest["pvalue"]
  score_p <- sum_c$sctest["pvalue"]

  lw <- 26
  vw <- 14
  linha <- paste0(strrep("\u2500", lw + vw + 3), "\n")

  .print_titulo("MÉTRICAS DO MODELO DE SOBREVIVÊNCIA (COX)")

  cat(sprintf("  %-*s %*s\n", lw, "Métrica", vw, "Valor"))
  cat("  ", linha, sep = "")
  cat(sprintf("  %-*s %*s\n", lw, "Concordância (C-index)", vw, .fmt_num(c_idx, 3)))
  cat(sprintf("  %-*s %*s\n", lw, "Erro Padrão do C-index",  vw, .fmt_num(se_c, 3)))
  cat(sprintf("  %-*s %*s\n", lw, "R² de Cox & Snell",       vw, .fmt_num(r2_cs, 3)))
  cat(sprintf("  %-*s %*s\n", lw, "R² de Nagelkerke",        vw, .fmt_num(r2_nag, 3)))
  cat(sprintf("  %-*s %*s\n", lw, "Log-Verossimilhança",     vw, .fmt_num(loglik_fit, 2)))
  cat(sprintf("  %-*s %*s\n", lw, "AIC",                     vw, .fmt_num(aic_val, 1)))
  cat(sprintf("  %-*s %*s\n", lw, "BIC",                     vw, .fmt_num(bic_val, 1)))
  cat(sprintf("  %-*s %*s\n", lw, "Razão Verossimilhança (p)", vw, if (lrt_p < 0.001) "< 0,001" else .fmt_num(lrt_p, 3)))
  cat(sprintf("  %-*s %*s\n", lw, "Teste de Wald (p)",       vw, if (wald_p < 0.001) "< 0,001" else .fmt_num(wald_p, 3)))
  cat(sprintf("  %-*s %*s\n", lw, "Teste Score (p)",         vw, if (score_p < 0.001) "< 0,001" else .fmt_num(score_p, 3)))
  .print_rodape()

  df_res <- data.frame(
    Metrica = c("C_index", "EP_C_index", "R2_Cox_Snell", "R2_Nagelkerke",
                "LogLik", "AIC", "BIC", "p_LRT", "p_Wald", "p_Score"),
    Valor   = c(c_idx, se_c, r2_cs, r2_nag, loglik_fit, aic_val, bic_val, lrt_p, wald_p, score_p),
    stringsAsFactors = FALSE
  )
  invisible(df_res)
}


# ─────────────────────────────────────────────────────────────────────────────
# 5. SELEÇÃO DE MODELOS DE COX (BEST SUBSETS)
# ─────────────────────────────────────────────────────────────────────────────

#' @title Seleção de Modelos de Riscos Proporcionais de Cox
#' @description Avalia todas as combinações de covariáveis para o modelo de Cox,
#'   ranqueando os melhores modelos pelos critérios de AIC, Delta AIC, pesos de Akaike
#'   e Índice C de Concordância.
#'
#' @param formula Fórmula do modelo inicial (ex: \code{Surv(tempo, evento) ~ .} ou \code{Surv(tempo, evento) ~ x1 + x2 + x3}).
#' @param dados Data frame contendo os dados.
#' @param top Número de modelos a serem exibidos no ranking (padrão \code{5}).
#' @return Data frame invisível com o ranking dos melhores modelos avaliados.
#' @export
selecao_modelos_cox <- function(formula, dados, top = 5) {

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("O pacote 'survival' é necessário para a seleção de modelos de Cox.", call. = FALSE)
  }

  if (!is.data.frame(dados)) {
    stop("O argumento 'dados' deve ser um data frame.", call. = FALSE)
  }

  # Extrai termos do lado direito da formula
  termos_todos <- attr(stats::terms(formula, data = dados), "term.labels")
  if (length(termos_todos) == 0) {
    stop("A fórmula deve conter ao menos uma covariável preditora.", call. = FALSE)
  }

  # Formula basica de resposta
  resp_call <- formula[[2]]
  resp_str  <- deparse(resp_call)

  # Limita a 10 termos para evitar explosao combinatoria
  if (length(termos_todos) > 10) {
    message("Mais de 10 variáveis detectadas. Selecionando os 10 preditores mais promissores para avaliação.")
    termos_todos <- termos_todos[1:10]
  }

  # Gera todas as combinacoes possiveis
  p_tot <- length(termos_todos)
  todas_combs <- list()
  for (k in 1:p_tot) {
    cb <- utils::combn(termos_todos, k, simplify = FALSE)
    todas_combs <- c(todas_combs, cb)
  }

  resultados <- list()

  for (i in seq_along(todas_combs)) {
    vars_i <- todas_combs[[i]]
    f_txt  <- sprintf("%s ~ %s", resp_str, paste(vars_i, collapse = " + "))
    f_obj  <- stats::as.formula(f_txt)

    fit_i <- tryCatch(
      survival::coxph(f_obj, data = dados),
      error = function(e) NULL
    )

    if (!is.null(fit_i)) {
      aic_i <- stats::AIC(fit_i)
      c_i   <- summary(fit_i)$concordance[1]
      resultados[[length(resultados) + 1]] <- data.frame(
        Modelo     = paste(vars_i, collapse = " + "),
        k          = length(vars_i),
        AIC        = aic_i,
        C_index    = c_i,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(resultados) == 0) {
    stop("Nenhum modelo de Cox convergiu satisfatoriamente.", call. = FALSE)
  }

  df_rank <- do.call(rbind, resultados)
  df_rank <- df_rank[order(df_rank$AIC), ]

  # Delta AIC e Pesos de Akaike
  min_aic <- min(df_rank$AIC)
  df_rank$Delta_AIC <- df_rank$AIC - min_aic
  pesos_raw <- exp(-0.5 * df_rank$Delta_AIC)
  df_rank$Peso_Akaike <- pesos_raw / sum(pesos_raw)
  rownames(df_rank) <- NULL

  # Exibicao do top
  n_show <- min(top, nrow(df_rank))
  top_df <- df_rank[1:n_show, ]

  .print_titulo("SELEÇÃO DE MODELOS DE COX (BEST SUBSETS)")
  cat(sprintf("  Critério de Ranqueamento: Menor AIC  |  Modelos testados: %d\n\n", nrow(df_rank)))

  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  Rank  Preditores Incluídos               k     AIC       Delta AIC   Peso   C-index \n")
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")

  for (i in 1:n_show) {
    rk_s   <- sprintf("%2d", i)
    mod_s  <- top_df$Modelo[i]
    if (nchar(mod_s) > 30) mod_s <- paste0(substr(mod_s, 1, 27), "...")
    k_s    <- sprintf("%2d", top_df$k[i])
    aic_s  <- sprintf("%8.2f", top_df$AIC[i])
    daic_s <- sprintf("%7.2f", top_df$Delta_AIC[i])
    w_s    <- sprintf("%6.3f", top_df$Peso_Akaike[i])
    c_s    <- sprintf("%5.3f", top_df$C_index[i])

    cat(sprintf("  %-4s  %-32s %s  %s    %s   %s   %s\n",
                rk_s, mod_s, k_s, aic_s, daic_s, w_s, c_s))
  }
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  * Delta AIC <= 2 indica modelos substancialmente equivalentes pelo critério de Burnham & Anderson.\n\n")

  # Modelo recomendado
  melhor_mod <- df_rank$Modelo[1]
  cat(sprintf("  ▶ Modelo Recomendado (Menor AIC): %s ~ %s\n", resp_str, melhor_mod))
  cat(sprintf("    AIC: %.2f  |  Concordância (C-index): %.3f\n\n",
              df_rank$AIC[1], df_rank$C_index[1]))

  .print_rodape()
  invisible(df_rank)
}

