# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# M\u00d3DULO: AN\u00c1LISE DE MEDIDAS REPETIDAS
# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500

#' @title An\u00e1lise de Medidas Repetidas (ANOVA Intrasujeito)
#' @description Executa ANOVA de medidas repetidas com teste de esfericidade de
#'   Mauchly, corre\u00e7\u00f5es de Greenhouse-Geisser e Huynh-Feldt quando necess\u00e1rio,
#'   c\u00e1lculo de eta\u00b2 e gr\u00e1fico de perfis longitudinais com bandas de erro.
#' @param formula F\u00f3rmula no formato \code{resposta ~ tempo} ou
#'   \code{resposta ~ tempo | grupo} para incluir fator entre-sujeitos.
#' @param sujeito Nome da coluna identificadora do sujeito (string ou s\u00edmbolo).
#' @param dados Data frame contendo os dados no formato longo.
#' @param confianca N\u00edvel de confian\u00e7a para intervalos (padr\u00e3o \code{0.95}).
#' @param grafico L\u00f3gico. Se \code{TRUE}, exibe gr\u00e1fico de perfis longitudinais.
#' @return Retorna invisivelmente uma lista com resultados da ANOVA, teste de
#'   Mauchly e epsilons de corre\u00e7\u00e3o.
#' @export
analise_medidas_repetidas <- function(formula,
                                      sujeito,
                                      dados,
                                      confianca = 0.95,
                                      grafico   = TRUE) {

  # ── Validar entradas ────────────────────────────────────────────────────────
  if (!inherits(formula, "formula"))
    stop("'formula' deve ser um objeto de f\u00f3rmula.", call. = FALSE)
  if (!is.data.frame(dados))
    stop("'dados' deve ser um data.frame.", call. = FALSE)

  # Resolver nome do sujeito
  suj_expr <- substitute(sujeito)
  suj_nome <- if (is.character(sujeito)) sujeito else deparse(suj_expr)
  if (!suj_nome %in% names(dados))
    stop(sprintf("A coluna '%s' n\u00e3o foi encontrada em 'dados'.", suj_nome), call. = FALSE)

  # Parsear f\u00f3rmula: resposta ~ tempo | grupo (grupo \u00e9 opcional)
  formula_str <- deparse(formula)
  tem_grupo   <- grepl("\\|", formula_str)

  if (tem_grupo) {
    partes      <- strsplit(formula_str, "\\|")[[1]]
    f_base      <- trimws(partes[1])         # "resposta ~ tempo"
    grupo_nome  <- trimws(partes[2])
    f_intra     <- as.formula(f_base)
  } else {
    f_intra     <- formula
    grupo_nome  <- NULL
  }

  resp_nome  <- as.character(f_intra[[2]])
  tempo_nome <- as.character(f_intra[[3]])

  for (col in c(resp_nome, tempo_nome, suj_nome)) {
    if (!col %in% names(dados))
      stop(sprintf("Coluna '%s' n\u00e3o encontrada em 'dados'.", col), call. = FALSE)
  }
  if (!is.null(grupo_nome) && !grupo_nome %in% names(dados))
    stop(sprintf("Coluna de grupo '%s' n\u00e3o encontrada em 'dados'.", grupo_nome), call. = FALSE)

  dados[[tempo_nome]]  <- as.factor(dados[[tempo_nome]])
  dados[[suj_nome]]    <- as.factor(dados[[suj_nome]])
  if (!is.null(grupo_nome))
    dados[[grupo_nome]] <- as.factor(dados[[grupo_nome]])

  niveis_tempo <- levels(dados[[tempo_nome]])
  k            <- length(niveis_tempo)
  n_suj        <- nlevels(dados[[suj_nome]])
  alpha        <- 1 - confianca

  # ── Ajustar ANOVA de medidas repetidas ──────────────────────────────────────
  if (!is.null(grupo_nome)) {
    f_aov <- as.formula(sprintf(
      "%s ~ %s * %s + Error(%s / %s)",
      resp_nome, grupo_nome, tempo_nome, suj_nome, tempo_nome
    ))
  } else {
    f_aov <- as.formula(sprintf(
      "%s ~ %s + Error(%s / %s)",
      resp_nome, tempo_nome, suj_nome, tempo_nome
    ))
  }

  aov_fit <- tryCatch(
    stats::aov(f_aov, data = dados),
    error = function(e) stop("Falha ao ajustar a ANOVA: ", e$message, call. = FALSE)
  )

  # ── Teste de Esfericidade de Mauchly ────────────────────────────────────────
  mauchly_res <- NULL
  W_mauchly   <- NA_real_
  p_mauchly   <- NA_real_
  chi_mauchly <- NA_real_
  gl_mauchly  <- NA_real_
  eps_gg      <- 1
  eps_hf      <- 1

  if (k >= 3) {
    # Construir matriz wide para mlm
    ids_unicos <- levels(dados[[suj_nome]])
    wide_mat   <- do.call(cbind, lapply(niveis_tempo, function(lv) {
      sub_d <- dados[dados[[tempo_nome]] == lv, ]
      sub_d <- sub_d[order(sub_d[[suj_nome]]), ]
      vals  <- sub_d[[resp_nome]]
      if (length(vals) != n_suj)
        stop("Dados desbalanceados: cada sujeito deve ter exatamente uma observa\u00e7\u00e3o por tempo.", call. = FALSE)
      vals
    }))
    colnames(wide_mat) <- niveis_tempo

    mlm_fit <- tryCatch({
      df_wide <- as.data.frame(wide_mat)
      lm(as.matrix(df_wide) ~ 1)
    }, error = function(e) NULL)

    if (!is.null(mlm_fit)) {
      mauchly_res <- tryCatch(
        stats::mauchly.test(mlm_fit, X = ~1),
        error = function(e) NULL
      )
      if (!is.null(mauchly_res)) {
        W_mauchly   <- unname(mauchly_res$statistic)
        p_mauchly   <- mauchly_res$p.value
        chi_mauchly <- -((n_suj - 1) - (2 * k^2 - 3 * k + 3) / (6 * (k - 1))) * log(W_mauchly)
        gl_mauchly  <- k * (k - 1) / 2 - 1
      }

      # Calcular epsilons manualmente
      tryCatch({
        S       <- cov(wide_mat)
        ones    <- rep(1, k)
        S_cent  <- S - rowMeans(S) %o% ones - ones %o% rowMeans(S) + mean(S)
        tr_S    <- sum(diag(S_cent))
        tr_S2   <- sum(S_cent^2)
        eps_gg  <- tr_S^2 / ((k - 1) * tr_S2)
        eps_gg  <- min(max(eps_gg, 1 / (k - 1)), 1)
        eps_hf_raw <- (n_suj * (k - 1) * eps_gg - 2) /
                      ((k - 1) * (n_suj - 1 - (k - 1) * eps_gg))
        eps_hf  <- min(1, max(eps_hf_raw, 1 / (k - 1)))
      }, error = function(e) {
        eps_gg <<- 1
        eps_hf <<- 1
      })
    }
  }

  # ── Extrair resultados da ANOVA ──────────────────────────────────────────────
  aov_summ <- summary(aov_fit)

  # Obter tabela do componente de tempo dentro de sujeitos
  nome_strato <- paste0("Error: ", suj_nome, ":", tempo_nome)
  # Fallback para outros nomes de estrato
  nomes_disponiveis <- names(aov_summ)
  idx_tempo <- which(nomes_disponiveis == nome_strato)
  if (length(idx_tempo) == 0) {
    # Tentar encontrar estrato com tempo_nome no nome
    idx_tempo <- which(grepl(tempo_nome, nomes_disponiveis, fixed = TRUE) &
                         grepl(suj_nome,   nomes_disponiveis, fixed = TRUE))
  }
  if (length(idx_tempo) == 0) idx_tempo <- length(nomes_disponiveis)

  tab_intra <- aov_summ[[idx_tempo[[1]]]][[1]]

  # Localizar linha do tempo e do res\u00edduo
  nomes_linhas <- rownames(tab_intra)
  idx_f_tempo  <- which(trimws(nomes_linhas) == tempo_nome)
  if (length(idx_f_tempo) == 0) idx_f_tempo <- 1
  idx_f_resid  <- which(grepl("Residuals", nomes_linhas))
  if (length(idx_f_resid) == 0) idx_f_resid <- nrow(tab_intra)

  sq_tempo  <- tab_intra[idx_f_tempo, "Sum Sq"]
  gl_tempo  <- tab_intra[idx_f_tempo, "Df"]
  qm_tempo  <- tab_intra[idx_f_tempo, "Mean Sq"]
  F_tempo   <- tab_intra[idx_f_tempo, "F value"]
  p_tempo   <- tab_intra[idx_f_tempo, "Pr(>F)"]

  sq_resid  <- tab_intra[idx_f_resid, "Sum Sq"]
  gl_resid  <- tab_intra[idx_f_resid, "Df"]
  qm_resid  <- sq_resid / gl_resid

  # eta\u00b2 parcial = SQ_tempo / (SQ_tempo + SQ_resid)
  eta2      <- sq_tempo / (sq_tempo + sq_resid)

  # p-valor corrigido (GG)
  if (!is.na(p_mauchly) && p_mauchly < 0.05 && k >= 3) {
    p_gg  <- pf(F_tempo,
                df1 = gl_tempo * eps_gg,
                df2 = gl_resid * eps_gg,
                lower.tail = FALSE)
    p_hf  <- pf(F_tempo,
                df1 = gl_tempo * eps_hf,
                df2 = gl_resid * eps_hf,
                lower.tail = FALSE)
  } else {
    p_gg  <- p_tempo
    p_hf  <- p_tempo
  }

  # ── Imprimir resultados ──────────────────────────────────────────────────────
  w_sep <- 75

  .print_titulo("AN\u00c1LISE DE MEDIDAS REPETIDAS")

  cat(sprintf(
    "  Vari\u00e1vel Resposta: %s   |   Fator Tempo: %s   |   Sujeitos: %s\n\n",
    resp_nome, tempo_nome, suj_nome
  ))
  if (!is.null(grupo_nome))
    cat(sprintf("  Fator entre-sujeitos: %s\n\n", grupo_nome))

  # ── Bloco: Esfericidade de Mauchly ──────────────────────────────────────────
  .print_topico("ESFERICIDADE DE MAUCHLY")
  if (k < 3) {
    cat("  Teste de Mauchly n\u00e3o aplic\u00e1vel (k < 3 n\u00edveis de tempo).\n\n")
  } else {
    wc_m <- c(18, 14, 10, 14, 16)
    hdr_m <- paste0(
      "  ",
      .th_pad("W de Mauchly",  wc_m[1], "center"),
      .th_pad("\u03c7\u00b2",              wc_m[2], "center"),
      .th_pad("GL",             wc_m[3], "center"),
      .th_pad("p-valor",        wc_m[4], "center"),
      .th_pad("Decis\u00e3o",        wc_m[5], "center")
    )
    sep_m <- paste0("  ", strrep("\u2500", nchar(hdr_m) - 2))
    cat(sep_m, "\n")
    cat(hdr_m, "\n")
    cat(sep_m, "\n")

    decisao_m <- if (!is.na(p_mauchly) && p_mauchly < 0.05)
      "Violada \u26a0" else "Atendida \u2713"

    cat(
      "  ",
      .th_pad(formatC(W_mauchly,   format = "f", digits = 3, decimal.mark = ","), wc_m[1], "center"),
      .th_pad(formatC(chi_mauchly, format = "f", digits = 3, decimal.mark = ","), wc_m[2], "center"),
      .th_pad(formatC(gl_mauchly,  format = "f", digits = 0),                     wc_m[3], "center"),
      .th_pad(.th_fmt_p(p_mauchly),                                                wc_m[4], "center"),
      .th_pad(decisao_m,                                                           wc_m[5], "center"),
      "\n", sep = ""
    )
    cat(sep_m, "\n\n")

    if (!is.na(p_mauchly) && p_mauchly < 0.05) {
      cat(sprintf(
        "  Epsilon GG: %.3f   |   Epsilon HF: %.3f\n\n",
        eps_gg, eps_hf
      ))
    }
  }

  # ── Bloco: Tabela ANOVA ──────────────────────────────────────────────────────
  .print_topico("ANOVA \u2014 EFEITO DO TEMPO DENTRO DOS SUJEITOS")
  wc_a <- c(22, 10, 8, 10, 12, 12, 8)
  hdr_a <- paste0(
    "  ",
    .th_pad("Fonte de Varia\u00e7\u00e3o", wc_a[1], "left"),
    .th_pad("SQ",           wc_a[2], "center"),
    .th_pad("GL",           wc_a[3], "center"),
    .th_pad("QM",           wc_a[4], "center"),
    .th_pad("F calc",       wc_a[5], "center"),
    .th_pad("p-valor",      wc_a[6], "center"),
    .th_pad("\u03b7\u00b2",           wc_a[7], "center")
  )
  sep_a <- paste0("  ", strrep("\u2500", nchar(hdr_a) - 2))
  cat(sep_a, "\n")
  cat(hdr_a, "\n")
  cat(sep_a, "\n")

  p_uso     <- if (!is.na(p_mauchly) && p_mauchly < 0.05) p_gg else p_tempo
  p_str_a   <- paste(.th_fmt_p(p_uso), .th_asterisk(p_uso))

  cat(
    "  ",
    .th_pad(tempo_nome,                                                           wc_a[1], "left"),
    .th_pad(formatC(sq_tempo, format = "f", digits = 2, decimal.mark = ","),     wc_a[2], "center"),
    .th_pad(formatC(gl_tempo, format = "f", digits = 0),                         wc_a[3], "center"),
    .th_pad(formatC(qm_tempo, format = "f", digits = 2, decimal.mark = ","),     wc_a[4], "center"),
    .th_pad(formatC(F_tempo,  format = "f", digits = 2, decimal.mark = ","),     wc_a[5], "center"),
    .th_pad(p_str_a,                                                              wc_a[6], "center"),
    .th_pad(formatC(eta2,     format = "f", digits = 2, decimal.mark = ","),     wc_a[7], "center"),
    "\n", sep = ""
  )
  cat(
    "  ",
    .th_pad(paste0("Res\u00edduo (", suj_nome, ")"),                             wc_a[1], "left"),
    .th_pad(formatC(sq_resid, format = "f", digits = 2, decimal.mark = ","),     wc_a[2], "center"),
    .th_pad(formatC(gl_resid, format = "f", digits = 0),                         wc_a[3], "center"),
    .th_pad(formatC(qm_resid, format = "f", digits = 2, decimal.mark = ","),     wc_a[4], "center"),
    "\n", sep = ""
  )
  cat(sep_a, "\n")
  .th_nota_sig(p_uso)
  cat("\n")

  # ── Bloco: Interpreta\u00e7\u00e3o ────────────────────────────────────────────────────
  .print_topico("INTERPRETA\u00c7\u00c3O")

  sig_str <- if (p_uso < alpha) {
    sprintf("SIGNIFICATIVO (%s) ao n\u00edvel de %.0f%%.", .th_fmt_p_interp(p_uso), alpha * 100)
  } else {
    sprintf("n\u00e3o significativo (%s) ao n\u00edvel de %.0f%%.", .th_fmt_p_interp(p_uso), alpha * 100)
  }
  cat(sprintf("  O efeito de '%s' foi %s\n", tempo_nome, sig_str))

  if (k >= 3 && !is.na(p_mauchly)) {
    if (p_mauchly < 0.05) {
      cat(sprintf(
        "  A esfericidade foi VIOLADA (Mauchly W = %.3f, %s).\n",
        W_mauchly, .th_fmt_p_interp(p_mauchly)
      ))
      cat(sprintf(
        "  Corre\u00e7\u00e3o de Greenhouse-Geisser aplicada (\u03b5 = %.3f); p-valor ajustado: %s\n",
        eps_gg, .th_fmt_p_interp(p_gg)
      ))
    } else {
      cat(sprintf(
        "  A esfericidade foi atendida (Mauchly W = %.3f, %s).\n",
        W_mauchly, .th_fmt_p_interp(p_mauchly)
      ))
    }
  }

  tamanho_str <- if (eta2 >= 0.14) "grande" else if (eta2 >= 0.06) "m\u00e9dio" else "pequeno"
  cat(sprintf("  Tamanho de efeito: \u03b7\u00b2 = %.3f (%s).\n\n", eta2, tamanho_str))

  .print_rodape()

  # ── Gr\u00e1fico de perfis longitudinais ────────────────────────────────────────
  if (grafico) {
    tryCatch({
      z_alpha <- qnorm(1 - alpha / 2)

      if (!is.null(grupo_nome)) {
        summ_d <- do.call(rbind, lapply(
          split(dados, list(dados[[tempo_nome]], dados[[grupo_nome]])),
          function(sub) {
            m <- mean(sub[[resp_nome]], na.rm = TRUE)
            s <- sd(sub[[resp_nome]],   na.rm = TRUE)
            n <- sum(!is.na(sub[[resp_nome]]))
            ep <- s / sqrt(n)
            data.frame(
              tempo  = as.character(sub[[tempo_nome]][1]),
              grupo  = as.character(sub[[grupo_nome]][1]),
              media  = m,
              ic_inf = m - z_alpha * ep,
              ic_sup = m + z_alpha * ep,
              stringsAsFactors = FALSE
            )
          }
        ))
        summ_d$tempo <- factor(summ_d$tempo, levels = niveis_tempo)

        p_graf <- ggplot2::ggplot(
          summ_d,
          ggplot2::aes(x = tempo, y = media, group = grupo, colour = grupo)
        ) +
          ggplot2::geom_ribbon(
            ggplot2::aes(ymin = ic_inf, ymax = ic_sup, fill = grupo),
            alpha = 0.15, colour = NA
          ) +
          ggplot2::geom_line(linewidth = 0.9) +
          ggplot2::geom_point(size = 3) +
          ggplot2::labs(
            title    = "Perfis Longitudinais por Grupo",
            subtitle = sprintf("IC %.0f%% da m\u00e9dia", confianca * 100),
            x        = tempo_nome,
            y        = paste("M\u00e9dia de", resp_nome),
            colour   = grupo_nome,
            fill     = grupo_nome
          )
      } else {
        summ_d <- do.call(rbind, lapply(
          split(dados, dados[[tempo_nome]]),
          function(sub) {
            m <- mean(sub[[resp_nome]], na.rm = TRUE)
            s <- sd(sub[[resp_nome]],   na.rm = TRUE)
            n <- sum(!is.na(sub[[resp_nome]]))
            ep <- s / sqrt(n)
            data.frame(
              tempo  = as.character(sub[[tempo_nome]][1]),
              media  = m,
              ic_inf = m - z_alpha * ep,
              ic_sup = m + z_alpha * ep,
              stringsAsFactors = FALSE
            )
          }
        ))
        summ_d$tempo <- factor(summ_d$tempo, levels = niveis_tempo)

        p_graf <- ggplot2::ggplot(
          summ_d,
          ggplot2::aes(x = tempo, y = media, group = 1)
        ) +
          ggplot2::geom_ribbon(
            ggplot2::aes(ymin = ic_inf, ymax = ic_sup),
            fill = "#377EB8", alpha = 0.2
          ) +
          ggplot2::geom_line(colour = "#377EB8", linewidth = 0.9) +
          ggplot2::geom_point(colour = "#377EB8", size = 3) +
          ggplot2::labs(
            title    = "Perfis Longitudinais",
            subtitle = sprintf("IC %.0f%% da m\u00e9dia", confianca * 100),
            x        = tempo_nome,
            y        = paste("M\u00e9dia de", resp_nome)
          )
      }

      tema_fn <- if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }
      p_graf <- p_graf + tema_fn
      print(p_graf)
    }, error = function(e) {
      message("  [aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico: ", e$message)
    })
  }

  # ── Retorno invis\u00edvel ────────────────────────────────────────────────────────
  invisible(list(
    anova        = tab_intra,
    mauchly_W    = W_mauchly,
    mauchly_p    = p_mauchly,
    epsilon_gg   = eps_gg,
    epsilon_hf   = eps_hf,
    eta2         = eta2,
    p_valor      = p_uso,
    p_gg         = p_gg,
    p_hf         = p_hf,
    aov_fit      = aov_fit
  ))
}
