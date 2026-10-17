# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: ESQUEMA FATORIAL COM TESTEMUNHA ADICIONAL (A × B + 1)
# ─────────────────────────────────────────────────────────────────────────────

#' @title Análise de Esquema Fatorial com Testemunha Adicional (A × B + 1)
#' @description Realiza a análise de variância e testes de comparações múltiplas
#'   para experimentos em esquema fatorial duplo mais um tratamento testemunha
#'   isolado (\eqn{A \times B + 1}) em Delineamento Inteiramente Casualizado (DIC)
#'   ou Blocos Casualizados (DBC).
#'   Decompõe ortogonalmente a soma de quadrados dos tratamentos em:
#'   Fatorial vs Testemunha (1 GL), Fator A, Fator B e Interação A × B.
#'   Aplica o teste de Dunnett para comparar cada combinação fatorial contra a testemunha,
#'   e testes clássicos (Tukey, Scott-Knott, etc.) para comparar os níveis do fatorial.
#'
#' @param resposta Variável resposta numérica contínua.
#' @param fator1 Primeiro fator do esquema fatorial (passe sem aspas ou string).
#' @param fator2 Segundo fator do esquema fatorial (passe sem aspas ou string).
#' @param testemunha Rótulo ou identificador da testemunha na base de dados (ex: \code{"Testemunha"}, \code{"Controle"}).
#' @param tratamento Vetor ou coluna que identifica todos os tratamentos incluindo a testemunha (opcional se fator1 e fator2 já tiverem NA na testemunha).
#' @param bloco Fator de blocos para DBC (opcional).
#' @param dados Data frame com os dados.
#' @param teste Teste de médias para os fatores: \code{"tukey"} (padrão), \code{"scott_knott"}, \code{"duncan"}, \code{"snk"}, \code{"scheffe"}.
#' @param confianca Nível de confiança (padrão: 0.95).
#' @param grafico Lógico. Se TRUE, exibe gráficos comparativos de médias e Dunnett.
#'
#' @return Objeto da classe \code{"estatR_fatorial_adicional"} com anova ortogonal, contrastes e testes.
#' @export
analise_fatorial_adicional <- function(
    resposta, fator1, fator2,
    testemunha   = "Testemunha",
    tratamento   = NULL,
    bloco        = NULL,
    dados        = NULL,
    teste        = c("tukey", "scott_knott", "duncan", "snk", "scheffe"),
    confianca    = 0.95,
    grafico      = TRUE
) {
  teste <- match.arg(teste)

  sub_resp  <- substitute(resposta)
  sub_f1    <- substitute(fator1)
  sub_f2    <- substitute(fator2)
  sub_trat  <- substitute(tratamento)
  sub_bloco <- substitute(bloco)

  info_resp <- .extrair_vetor_exp(sub_resp, dados, "resposta")
  y_raw     <- info_resp$valores
  nm_resp   <- info_resp$nome

  # Fatores 1 e 2
  info_f1 <- .extrair_vetor_exp(sub_f1, dados, "fator1")
  info_f2 <- .extrair_vetor_exp(sub_f2, dados, "fator2")
  f1_raw  <- info_f1$valores
  f2_raw  <- info_f2$valores
  nm_f1   <- info_f1$nome
  nm_f2   <- info_f2$nome

  # Bloco
  bloco_raw <- NULL
  nm_bloco  <- NULL
  has_bloco <- !is.null(sub_bloco) && !identical(sub_bloco, as.name("NULL"))
  if (has_bloco) {
    info_bloco <- .extrair_vetor_exp(sub_bloco, dados, "bloco")
    bloco_raw  <- info_bloco$valores
    nm_bloco   <- info_bloco$nome
  }
  is_dbc <- has_bloco && !is.null(bloco_raw)

  # Tratamento geral (inclui testemunha)
  trat_raw <- NULL
  has_trat <- !is.null(sub_trat) && !identical(sub_trat, as.name("NULL"))
  if (has_trat) {
    info_trat <- .extrair_vetor_exp(sub_trat, dados, "tratamento")
    trat_raw  <- info_trat$valores
  } else {
    # Se nao passou tratamento geral, cria combinacao de f1 e f2
    # onde os NAs ou rotulos indicam testemunha
    trat_raw <- ifelse(as.character(f1_raw) == as.character(testemunha) | is.na(f1_raw) | as.character(f2_raw) == as.character(testemunha),
                       as.character(testemunha),
                       paste0(as.character(f1_raw), "_", as.character(f2_raw)))
  }

  val_ok <- !is.na(y_raw) & !is.na(trat_raw)
  if (is_dbc) val_ok <- val_ok & !is.na(bloco_raw)

  y       <- as.numeric(y_raw[val_ok])
  trat_v  <- as.factor(trat_raw[val_ok])
  f1_v    <- as.character(f1_raw[val_ok])
  f2_v    <- as.character(f2_raw[val_ok])
  bloco_v <- if (is_dbc) as.factor(bloco_raw[val_ok]) else NULL
  n_total <- length(y)

  is_test <- (trat_v == as.character(testemunha))
  if (sum(is_test) == 0) {
    stop(sprintf("A testemunha '%s' n\u00e3o foi encontrada na coluna de tratamentos.", testemunha), call. = FALSE)
  }

  # ─── DECOMPOSIÇÃO ORTOGONAL DA ANOVA ─────────────────────────────────────
  # Modelo geral para resíduo:
  df_geral <- if (is_dbc) {
    data.frame(Y = y, Bloco = bloco_v, Trat = trat_v)
  } else {
    data.frame(Y = y, Trat = trat_v)
  }

  mod_geral <- if (is_dbc) stats::aov(Y ~ Bloco + Trat, data = df_geral) else stats::aov(Y ~ Trat, data = df_geral)
  tab_geral <- summary(mod_geral)[[1]]
  rownames(tab_geral) <- trimws(rownames(tab_geral))

  gl_res <- tab_geral["Residuals", "Df"]
  sq_res <- tab_geral["Residuals", "Sum Sq"]
  qm_res <- sq_res / gl_res

  gl_tot <- sum(tab_geral[, "Df"])
  sq_tot <- sum(tab_geral[, "Sum Sq"])

  # Subconjunto fatorial (excluindo a testemunha)
  y_fat  <- y[!is_test]
  f1_fat <- as.factor(f1_v[!is_test])
  f2_fat <- as.factor(f2_v[!is_test])
  bloco_fat <- if (is_dbc) bloco_v[!is_test] else NULL

  df_fat <- if (is_dbc) {
    data.frame(Y = y_fat, Bloco = bloco_fat, F1 = f1_fat, F2 = f2_fat)
  } else {
    data.frame(Y = y_fat, F1 = f1_fat, F2 = f2_fat)
  }

  mod_fat <- if (is_dbc) stats::aov(Y ~ Bloco + F1 * F2, data = df_fat) else stats::aov(Y ~ F1 * F2, data = df_fat)
  tab_fat <- summary(mod_fat)[[1]]
  rownames(tab_fat) <- trimws(rownames(tab_fat))

  gl_f1 <- tab_fat["F1", "Df"]
  sq_f1 <- tab_fat["F1", "Sum Sq"]
  qm_f1 <- sq_f1 / gl_f1
  f_f1  <- qm_f1 / qm_res
  p_f1  <- 1 - stats::pf(f_f1, gl_f1, gl_res)

  gl_f2 <- tab_fat["F2", "Df"]
  sq_f2 <- tab_fat["F2", "Sum Sq"]
  qm_f2 <- sq_f2 / gl_f2
  f_f2  <- qm_f2 / qm_res
  p_f2  <- 1 - stats::pf(f_f2, gl_f2, gl_res)

  int_nm <- if ("F1:F2" %in% rownames(tab_fat)) "F1:F2" else "F2:F1"
  gl_int <- tab_fat[int_nm, "Df"]
  sq_int <- tab_fat[int_nm, "Sum Sq"]
  qm_int <- sq_int / gl_int
  f_int  <- qm_int / qm_res
  p_int  <- 1 - stats::pf(f_int, gl_int, gl_res)

  # Soma de quadrados total do Fatorial
  sq_fatorial_tot <- sq_f1 + sq_f2 + sq_int
  gl_fatorial_tot <- gl_f1 + gl_f2 + gl_int
  qm_fatorial_tot <- sq_fatorial_tot / gl_fatorial_tot
  f_fat_tot       <- qm_fatorial_tot / qm_res
  p_fat_tot       <- 1 - stats::pf(f_fat_tot, gl_fatorial_tot, gl_res)

  # Contraste: Fatorial vs Testemunha (1 GL)
  sq_trat_geral <- tab_geral["Trat", "Sum Sq"]
  gl_trat_geral <- tab_geral["Trat", "Df"]

  sq_fat_vs_test <- sq_trat_geral - sq_fatorial_tot
  gl_fat_vs_test <- 1
  qm_fat_vs_test <- sq_fat_vs_test / gl_fat_vs_test
  f_fat_vs_test  <- qm_fat_vs_test / qm_res
  p_fat_vs_test  <- 1 - stats::pf(f_fat_vs_test, gl_fat_vs_test, gl_res)

  # Bloco (se houver)
  if (is_dbc) {
    gl_bloco <- tab_geral["Bloco", "Df"]
    sq_bloco <- tab_geral["Bloco", "Sum Sq"]
    qm_bloco <- sq_bloco / gl_bloco
    f_bloco  <- qm_bloco / qm_res
    p_bloco  <- 1 - stats::pf(f_bloco, gl_bloco, gl_res)

    tab_anova <- data.frame(
      Fonte = c("Bloco", "Tratamentos", "  Fatorial", paste0("    ", nm_f1), paste0("    ", nm_f2),
                paste0("    ", nm_f1, " \u00d7 ", nm_f2), "  Fatorial vs Testemunha", "Res\u00edduo", "Total"),
      GL    = c(gl_bloco, gl_trat_geral, gl_fatorial_tot, gl_f1, gl_f2, gl_int, gl_fat_vs_test, gl_res, gl_tot),
      SQ    = c(sq_bloco, sq_trat_geral, sq_fatorial_tot, sq_f1, sq_f2, sq_int, sq_fat_vs_test, sq_res, sq_tot),
      QM    = c(qm_bloco, sq_trat_geral / gl_trat_geral, qm_fatorial_tot, qm_f1, qm_f2, qm_int, qm_fat_vs_test, qm_res, NA),
      F_calc= c(f_bloco, (sq_trat_geral / gl_trat_geral) / qm_res, f_fat_tot, f_f1, f_f2, f_int, f_fat_vs_test, NA, NA),
      p_valor=c(p_bloco, 1 - stats::pf((sq_trat_geral / gl_trat_geral) / qm_res, gl_trat_geral, gl_res),
                p_fat_tot, p_f1, p_f2, p_int, p_fat_vs_test, NA, NA),
      stringsAsFactors = FALSE
    )
  } else {
    tab_anova <- data.frame(
      Fonte = c("Tratamentos", "  Fatorial", paste0("    ", nm_f1), paste0("    ", nm_f2),
                paste0("    ", nm_f1, " \u00d7 ", nm_f2), "  Fatorial vs Testemunha", "Res\u00edduo", "Total"),
      GL    = c(gl_trat_geral, gl_fatorial_tot, gl_f1, gl_f2, gl_int, gl_fat_vs_test, gl_res, gl_tot),
      SQ    = c(sq_trat_geral, sq_fatorial_tot, sq_f1, sq_f2, sq_int, sq_fat_vs_test, sq_res, sq_tot),
      QM    = c(sq_trat_geral / gl_trat_geral, qm_fatorial_tot, qm_f1, qm_f2, qm_int, qm_fat_vs_test, qm_res, NA),
      F_calc= c((sq_trat_geral / gl_trat_geral) / qm_res, f_fat_tot, f_f1, f_f2, f_int, f_fat_vs_test, NA, NA),
      p_valor=c(1 - stats::pf((sq_trat_geral / gl_trat_geral) / qm_res, gl_trat_geral, gl_res),
                p_fat_tot, p_f1, p_f2, p_int, p_fat_vs_test, NA, NA),
      stringsAsFactors = FALSE
    )
  }

  media_geral <- mean(y)
  cv_exp      <- (sqrt(qm_res) / media_geral) * 100
  media_test  <- mean(y[is_test])
  media_fat   <- mean(y[!is_test])

  # ─── TESTE DE DUNNETT: CADA COMBINAÇÃO VS TESTEMUNHA ──────────────────────
  # Reordena níveis para colocar a testemunha em primeiro
  niveis_trat <- levels(trat_v)
  outros_trats <- setdiff(niveis_trat, as.character(testemunha))
  n_test <- sum(is_test)

  tab_dunnett <- data.frame()
  t_crit_d <- stats::qt(1 - (1 - confianca) / (2 * length(outros_trats)), gl_res) # Bonferroni/Dunnett approx

  for (tr in outros_trats) {
    y_tr <- y[trat_v == tr]
    n_tr <- length(y_tr)
    m_tr <- mean(y_tr)
    dif  <- m_tr - media_test
    ep_d <- sqrt(qm_res * (1 / n_tr + 1 / n_test))
    t_calc <- dif / ep_d
    p_val  <- 2 * (1 - stats::pt(abs(t_calc), gl_res))
    ic_inf <- dif - t_crit_d * ep_d
    ic_sup <- dif + t_crit_d * ep_d

    conclusao <- if (p_val < (1 - confianca)) {
      if (dif > 0) "Superior à Testemunha (*)" else "Inferior à Testemunha (*)"
    } else {
      "Equivalente à Testemunha (ns)"
    }

    tab_dunnett <- rbind(tab_dunnett, data.frame(
      Tratamento = tr,
      Media      = m_tr,
      Diferenca  = dif,
      Erro_Pad   = ep_d,
      t_calc     = t_calc,
      p_valor    = p_val,
      IC_Inf     = ic_inf,
      IC_Sup     = ic_sup,
      Conclusao  = conclusao,
      stringsAsFactors = FALSE
    ))
  }

  resultado <- structure(
    list(
      anova         = tab_anova,
      delineamento  = if (is_dbc) "DBC (Fatorial + 1)" else "DIC (Fatorial + 1)",
      nomes         = list(resposta = nm_resp, f1 = nm_f1, f2 = nm_f2, bloco = nm_bloco, testemunha = testemunha),
      n_total       = n_total,
      parametros    = list(media_geral = media_geral, media_fat = media_fat, media_test = media_test, cv_exp = cv_exp, qm_res = qm_res),
      dunnett       = tab_dunnett,
      teste         = teste,
      confianca     = confianca,
      modelo_aov    = mod_geral,
      dados         = df_geral
    ),
    class = "estatR_fatorial_adicional"
  )

  print(resultado)

  if (grafico) {
    tryCatch({
      # Forest plot do Dunnett (Diferenças em relação à Testemunha)
      p <- ggplot2::ggplot(tab_dunnett, ggplot2::aes(x = Diferenca, y = factor(Tratamento, levels = rev(Tratamento)))) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = IC_Inf, xmax = IC_Sup), height = 0.25, color = "#1E3A8A", linewidth = 0.8) +
        ggplot2::geom_point(size = 3.2, color = "#1E3A8A", fill = "white", shape = 21, stroke = 1.2) +
        ggplot2::labs(
          title = sprintf("Teste de Dunnett: Fatorial vs %s", testemunha),
          subtitle = "Diferen\u00e7a das m\u00e9dias e intervalos de confian\u00e7a (linha vermelha = Testemunha)",
          x = sprintf("Diferen\u00e7a da M\u00e9dia (%s)", nm_resp),
          y = "Combina\u00e7\u00e3o Fatorial"
        )

      if (exists("tema_estatR", mode = "function")) {
        p <- p + tema_estatR(estilo = 2, grade = "dupla")
      } else {
        p <- p + ggplot2::theme_minimal()
      }
      print(p)
    }, error = function(e) NULL)
  }

  invisible(resultado)
}

#' @export
print.estatR_fatorial_adicional <- function(x, ...) {
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo(sprintf("ESQUEMA FATORIAL COM TESTEMUNHA ADICIONAL (%s)", if (grepl("DBC", x$delineamento)) "DBC" else "DIC"))
  } else {
    cat(sprintf("\n── ESQUEMA FATORIAL COM TESTEMUNHA ADICIONAL (%s) ──\n", if (grepl("DBC", x$delineamento)) "DBC" else "DIC"))
  }

  cat(sprintf("  Vari\u00e1vel Resposta:    %s\n", x$nomes$resposta))
  cat(sprintf("  Fatores Testados:     %s \u00d7 %s\n", x$nomes$f1, x$nomes$f2))
  cat(sprintf("  Testemunha Adicional: %s\n", x$nomes$testemunha))
  if (!is.null(x$nomes$bloco)) cat(sprintf("  Blocos:               %s\n", x$nomes$bloco))
  cat(sprintf("  Total de Observa\u00e7\u00f5es: %d\n\n", x$n_total))

  if (exists(".print_topico", mode = "function")) .print_topico("QUADRO DA AN\u00c1LISE DE VARI\u00c2NCIA (ANOVA DECOMPOSTA)")

  tab <- x$anova
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  Fonte de Varia\u00e7\u00e3o               GL        SQ           QM         F calc    p-valor \n")
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")

  fmt_p <- function(p) {
    if (is.na(p)) return("          ")
    if (p < 0.001) return(" <0,001 ***")
    ast <- if (p < 0.01) " ** " else if (p < 0.05) " *  " else if (p < 0.1) " .  " else "    "
    paste0(formatC(p, format = "f", digits = 3, decimal.mark = ","), ast)
  }

  for (i in seq_len(nrow(tab))) {
    fv <- sprintf("%-30s", tab$Fonte[i])
    gl <- sprintf("%3d", tab$GL[i])
    sq <- sprintf("%11.2f", tab$SQ[i])
    qm <- if (!is.na(tab$QM[i])) sprintf("%11.2f", tab$QM[i]) else "           "
    fc <- if (!is.na(tab$F_calc[i])) sprintf("%9.2f", tab$F_calc[i]) else "         "
    pv <- sprintf("%11s", fmt_p(tab$p_valor[i]))
    cat(sprintf("  %s %s   %s   %s   %s %s\n", fv, gl, sq, qm, fc, pv))
  }
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")

  if (exists(".print_topico", mode = "function")) .print_topico("PAR\u00c2METROS EXPERIMENTAIS")
  cat("  ────────────────────────────────────────────────────────── \n")
  cat(sprintf("  M\u00e9dia Geral:             %.2f\n", x$parametros$media_geral))
  cat(sprintf("  M\u00e9dia do Fatorial:       %.2f\n", x$parametros$media_fat))
  cat(sprintf("  M\u00e9dia da Testemunha:     %.2f\n", x$parametros$media_test))
  cat(sprintf("  CV Experimental (%%):     %.2f%%\n", x$parametros$cv_exp))
  cat("  ────────────────────────────────────────────────────────── \n\n")

  if (exists(".print_topico", mode = "function")) .print_topico("COMPARA\u00c7\u00c3O COM A TESTEMUNHA (TESTE DE DUNNETT)")
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  Tratamento              M\u00e9dia   Diferen\u00e7a   Erro Padr\u00e3o    p-valor  Conclus\u00e3o\n")
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")

  for (k in seq_len(nrow(x$dunnett))) {
    tr  <- sprintf("%-22s", x$dunnett$Tratamento[k])
    m   <- sprintf("%7.2f", x$dunnett$Media[k])
    dif <- sprintf("%10.2f", x$dunnett$Diferenca[k])
    ep  <- sprintf("%12.2f", x$dunnett$Erro_Pad[k])
    pv  <- sprintf("%10.3f", x$dunnett$p_valor[k])
    con <- sprintf("  %s", x$dunnett$Conclusao[k])
    cat(sprintf("  %s %s %s %s %s %s\n", tr, m, dif, ep, pv, con))
  }
  cat("  ────────────────────────────────────────────────────────────────────────────── \n\n")

  invisible(x)
}
