# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: ANÁLISE FATORIAL DUPLO (A × B) — DIC e DBC
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# HELPER INTERNO: comparação múltipla dentro de subgrupos
# ─────────────────────────────────────────────────────────────────────────────

#' @keywords internal
.comparar_dentro <- function(y_sub, fator_sub, qm_res_global, gl_res_global,
                              teste = "tukey", confianca = 0.95, upper = FALSE) {
  # Retorna vetor nomeado de letras para os niveis de fator_sub
  # usando qm_res_global (residuo do modelo fatorial completo)
  alpha    <- 1 - confianca
  fator_f  <- as.factor(fator_sub)
  niveis_s <- levels(fator_f)
  k        <- length(niveis_s)

  medias_s <- tapply(y_sub, fator_f, mean)
  ns_s     <- as.vector(table(fator_f))
  names(ns_s) <- niveis_s
  r_media  <- mean(ns_s)

  if (k == 1) {
    letras <- stats::setNames("a", niveis_s)
    if (upper) letras <- toupper(letras)
    return(letras)
  }

  # Matriz de p-valores
  p_mat <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis_s, niveis_s))

  if (teste == "tukey") {
    # Usando TukeyHSD via modelo local (mas residuo global)
    # Ajusta modelo local
    df_loc  <- data.frame(Y = y_sub, Trat = fator_f)
    # Para usar qm_res global, calcula estatistica manualmente
    pares_lista <- utils::combn(niveis_s, 2, simplify = FALSE)
    for (pr in pares_lista) {
      n1 <- pr[1]; n2 <- pr[2]
      diff_v <- medias_s[n1] - medias_s[n2]
      ep_v   <- sqrt(qm_res_global * (1/ns_s[n1] + 1/ns_s[n2]))
      q_obs  <- abs(diff_v) / (ep_v / sqrt(2))
      p_v    <- 1 - stats::ptukey(q_obs, k, gl_res_global)
      p_mat[n1, n2] <- p_v
      p_mat[n2, n1] <- p_v
    }
  } else if (teste == "bonferroni") {
    pares_lista <- utils::combn(niveis_s, 2, simplify = FALSE)
    m_comp <- length(pares_lista)
    for (pr in pares_lista) {
      n1 <- pr[1]; n2 <- pr[2]
      diff_v <- medias_s[n1] - medias_s[n2]
      ep_v   <- sqrt(qm_res_global * (1/ns_s[n1] + 1/ns_s[n2]))
      t_v    <- abs(diff_v) / ep_v
      p_v    <- min(1, 2 * (1 - stats::pt(t_v, gl_res_global)) * m_comp)
      p_mat[n1, n2] <- p_v
      p_mat[n2, n1] <- p_v
    }
  } else if (teste %in% c("duncan", "snk")) {
    med_ord <- sort(medias_s, decreasing = TRUE)
    nomes_ord <- names(med_ord)
    for (i in 1:(k - 1)) {
      for (j in (i + 1):k) {
        n1 <- nomes_ord[i]; n2 <- nomes_ord[j]
        dist_p <- j - i + 1
        diff_v <- med_ord[i] - med_ord[j]
        q_obs  <- (diff_v / sqrt(qm_res_global / r_media)) * sqrt(2)
        p_v    <- 1 - stats::ptukey(abs(q_obs), dist_p, gl_res_global)
        p_mat[n1, n2] <- p_v
        p_mat[n2, n1] <- p_v
      }
    }
  } else if (teste == "scheffe") {
    pares_lista <- utils::combn(niveis_s, 2, simplify = FALSE)
    for (pr in pares_lista) {
      n1 <- pr[1]; n2 <- pr[2]
      diff_v <- medias_s[n1] - medias_s[n2]
      ep_v   <- sqrt(qm_res_global * (1/ns_s[n1] + 1/ns_s[n2]))
      t_v    <- diff_v / ep_v
      f_v    <- (t_v^2) / (k - 1)
      p_v    <- 1 - stats::pf(f_v, k - 1, gl_res_global)
      p_mat[n1, n2] <- p_v
      p_mat[n2, n1] <- p_v
    }
  } else {
    # Fallback: tukey
    pares_lista <- utils::combn(niveis_s, 2, simplify = FALSE)
    for (pr in pares_lista) {
      n1 <- pr[1]; n2 <- pr[2]
      diff_v <- medias_s[n1] - medias_s[n2]
      ep_v   <- sqrt(qm_res_global * (1/ns_s[n1] + 1/ns_s[n2]))
      q_obs  <- abs(diff_v) / (ep_v / sqrt(2))
      p_v    <- 1 - stats::ptukey(q_obs, k, gl_res_global)
      p_mat[n1, n2] <- p_v
      p_mat[n2, n1] <- p_v
    }
  }

  med_ord  <- sort(medias_s, decreasing = TRUE)
  letras   <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
  # Reordenar para os niveis originais
  letras   <- letras[niveis_s]
  if (upper) letras <- toupper(letras)
  letras
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Análise Fatorial Duplo (A × B) — DIC e DBC
#' @description Executa a Análise de Variância (ANOVA) para experimentos em
#'   esquema fatorial duplo (A × B) nos delineamentos DIC ou DBC.
#'   Quando a interação é significativa, realiza o desdobramento dos fatores
#'   e apresenta uma tabela matricial combinada com letras minúsculas (comparação
#'   de A dentro de B) e letras maiúsculas (comparação de B dentro de A).
#'   Quando a interação não é significativa, os efeitos principais são avaliados
#'   de forma independente.
#'
#' @param resposta Variável resposta contínua (com ou sem aspas).
#' @param fator1 Primeiro fator (com ou sem aspas). Aparece nas colunas da tabela.
#' @param fator2 Segundo fator (com ou sem aspas). Aparece nas linhas da tabela.
#' @param bloco Fator de blocos (com ou sem aspas). Se fornecido, adota o DBC.
#' @param dados Data frame opcional contendo as variáveis.
#' @param delineamento \code{"auto"} (padrão), \code{"dic"} ou \code{"dbc"}.
#'   No modo \code{"auto"}, usa DBC se \code{bloco} for fornecido, DIC caso contrário.
#' @param teste Teste de comparações múltiplas: \code{"tukey"} (padrão),
#'   \code{"scott_knott"}, \code{"duncan"}, \code{"snk"}, \code{"scheffe"} ou
#'   \code{"bonferroni"}.
#' @param confianca Nível de confiança (padrão \code{0.95}).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o gráfico de perfil de interação.
#' @param residuos Lógico. Se \code{TRUE}, executa a análise de resíduos.
#'
#' @return Objeto da classe \code{"estatR_fatorial"} com ANOVA, desdobramentos,
#'   tabela matricial de médias e gráfico de interação.
#' @export
#' @aliases fatorial fat
analise_fatorial_duplo <- function(resposta,
                              fator1,
                              fator2,
                              bloco       = NULL,
                              dados       = NULL,
                              delineamento = c("auto", "dic", "dbc"),
                              teste       = c("tukey", "scott_knott", "duncan", "snk", "scheffe", "bonferroni"),
                              confianca   = 0.95,
                              grafico     = TRUE,
                              residuos    = FALSE) {

  delineamento <- match.arg(delineamento)
  teste        <- match.arg(teste)

  sub_resp  <- substitute(resposta)
  sub_f1    <- substitute(fator1)
  sub_f2    <- substitute(fator2)
  sub_bloco <- substitute(bloco)

  info_resp  <- .extrair_vetor_exp(sub_resp,  dados, "resposta")
  info_f1    <- .extrair_vetor_exp(sub_f1,    dados, "fator1")
  info_f2    <- .extrair_vetor_exp(sub_f2,    dados, "fator2")

  y_raw   <- info_resp$valores
  f1_raw  <- info_f1$valores
  f2_raw  <- info_f2$valores
  nm_resp <- info_resp$nome
  nm_f1   <- info_f1$nome
  nm_f2   <- info_f2$nome

  bloco_raw <- NULL
  nm_bloco  <- NULL
  has_bloco <- !is.null(sub_bloco) && !identical(sub_bloco, as.name("NULL"))
  if (has_bloco) {
    info_bloco <- .extrair_vetor_exp(sub_bloco, dados, "bloco")
    bloco_raw  <- info_bloco$valores
    nm_bloco   <- info_bloco$nome
  }

  # Delineamento efetivo
  delin_efetivo <- if (delineamento == "auto") {
    if (has_bloco && !is.null(bloco_raw)) "dbc" else "dic"
  } else {
    delineamento
  }

  # Limpeza de NAs
  if (delin_efetivo == "dbc" && !is.null(bloco_raw)) {
    val_ok <- !is.na(y_raw) & !is.na(f1_raw) & !is.na(f2_raw) & !is.na(bloco_raw)
  } else {
    val_ok <- !is.na(y_raw) & !is.na(f1_raw) & !is.na(f2_raw)
  }
  y       <- as.numeric(y_raw[val_ok])
  f1_vec  <- as.factor(f1_raw[val_ok])
  f2_vec  <- as.factor(f2_raw[val_ok])
  bloco_v <- if (!is.null(bloco_raw)) as.factor(bloco_raw[val_ok]) else NULL
  n_total <- length(y)

  niveis_f1 <- levels(f1_vec)
  niveis_f2 <- levels(f2_vec)
  k1        <- length(niveis_f1)
  k2        <- length(niveis_f2)
  alpha     <- 1 - confianca

  if (k1 < 2) stop("Fator1 deve ter ao menos 2 níveis.", call. = FALSE)
  if (k2 < 2) stop("Fator2 deve ter ao menos 2 níveis.", call. = FALSE)

  # ─── AJUSTE DO MODELO ANOVA ───────────────────────────────────────────────
  if (delin_efetivo == "dbc") {
    df_mod  <- data.frame(Y = y, F1 = f1_vec, F2 = f2_vec, Bloco = bloco_v)
    mod_aov <- stats::aov(Y ~ Bloco + F1 * F2, data = df_mod)
    tab_aov <- summary(mod_aov)[[1]]
    rownames(tab_aov) <- trimws(rownames(tab_aov))

    gl_bloco <- tab_aov["Bloco", "Df"]
    sq_bloco <- tab_aov["Bloco", "Sum Sq"]
    qm_bloco <- tab_aov["Bloco", "Mean Sq"]
    f_bloco  <- tab_aov["Bloco", "F value"]
    p_bloco  <- tab_aov["Bloco", "Pr(>F)"]
  } else {
    df_mod  <- data.frame(Y = y, F1 = f1_vec, F2 = f2_vec)
    mod_aov <- stats::aov(Y ~ F1 * F2, data = df_mod)
    tab_aov <- summary(mod_aov)[[1]]
    rownames(tab_aov) <- trimws(rownames(tab_aov))
  }

  gl_f1  <- tab_aov["F1",     "Df"]
  sq_f1  <- tab_aov["F1",     "Sum Sq"]
  qm_f1  <- tab_aov["F1",     "Mean Sq"]
  f_f1   <- tab_aov["F1",     "F value"]
  p_f1   <- tab_aov["F1",     "Pr(>F)"]

  gl_f2  <- tab_aov["F2",     "Df"]
  sq_f2  <- tab_aov["F2",     "Sum Sq"]
  qm_f2  <- tab_aov["F2",     "Mean Sq"]
  f_f2   <- tab_aov["F2",     "F value"]
  p_f2   <- tab_aov["F2",     "Pr(>F)"]

  int_row <- if ("F1:F2" %in% rownames(tab_aov)) "F1:F2" else "F2:F1"
  gl_int <- tab_aov[int_row, "Df"]
  sq_int <- tab_aov[int_row, "Sum Sq"]
  qm_int <- tab_aov[int_row, "Mean Sq"]
  f_int  <- tab_aov[int_row, "F value"]
  p_int  <- tab_aov[int_row, "Pr(>F)"]

  gl_res <- tab_aov["Residuals", "Df"]
  sq_res <- tab_aov["Residuals", "Sum Sq"]
  qm_res <- tab_aov["Residuals", "Mean Sq"]

  gl_tot <- sum(tab_aov[, "Df"])
  sq_tot <- sum(tab_aov[, "Sum Sq"])

  media_geral <- mean(y)
  cv_exp      <- (sqrt(qm_res) / media_geral) * 100

  # R² do modelo
  if (delin_efetivo == "dbc") {
    sq_model <- sq_bloco + sq_f1 + sq_f2 + sq_int
  } else {
    sq_model <- sq_f1 + sq_f2 + sq_int
  }
  r2 <- sq_model / sq_tot

  # Tabela ANOVA formatada
  nm_int <- paste0(nm_f1, " \u00d7 ", nm_f2)  # × symbol

  if (delin_efetivo == "dbc") {
    tab_anova_out <- data.frame(
      Fonte   = c("Bloco", nm_f1, nm_f2, nm_int, "Res\u00edduo", "Total"),
      GL      = c(gl_bloco, gl_f1, gl_f2, gl_int, gl_res, gl_tot),
      SQ      = c(sq_bloco, sq_f1, sq_f2, sq_int, sq_res, sq_tot),
      QM      = c(qm_bloco, qm_f1, qm_f2, qm_int, qm_res, NA),
      F_calc  = c(f_bloco,  f_f1,  f_f2,  f_int,  NA,     NA),
      p_valor = c(p_bloco,  p_f1,  p_f2,  p_int,  NA,     NA),
      stringsAsFactors = FALSE
    )
  } else {
    tab_anova_out <- data.frame(
      Fonte   = c(nm_f1, nm_f2, nm_int, "Res\u00edduo", "Total"),
      GL      = c(gl_f1, gl_f2, gl_int, gl_res, gl_tot),
      SQ      = c(sq_f1, sq_f2, sq_int, sq_res, sq_tot),
      QM      = c(qm_f1, qm_f2, qm_int, qm_res, NA),
      F_calc  = c(f_f1,  f_f2,  f_int,  NA,     NA),
      p_valor = c(p_f1,  p_f2,  p_int,  NA,     NA),
      stringsAsFactors = FALSE
    )
  }

  # ─── MÉDIAS DA COMBINAÇÃO ─────────────────────────────────────────────────
  medias_comb <- tapply(y, list(F1 = f1_vec, F2 = f2_vec), mean)
  dp_comb     <- tapply(y, list(F1 = f1_vec, F2 = f2_vec), stats::sd)
  n_comb      <- tapply(y, list(F1 = f1_vec, F2 = f2_vec), length)
  medias_f1   <- tapply(y, f1_vec, mean)
  medias_f2   <- tapply(y, f2_vec, mean)

  # ─── DESDOBRAMENTO ────────────────────────────────────────────────────────
  interacao_sig <- !is.na(p_int) && p_int < alpha

  letras_min <- NULL  # F1 dentro de cada nivel de F2 (minúsculas, por linha)
  letras_mai <- NULL  # F2 dentro de cada nivel de F1 (maiúsculas, por coluna)
  letras_f1  <- NULL  # F1 efeito principal (interação n.s.)
  letras_f2  <- NULL  # F2 efeito principal (interação n.s.)

  desdobramento_f1_em_f2 <- NULL  # data.frame: F-teste de F1 dentro de cada B
  desdobramento_f2_em_f1 <- NULL  # data.frame: F-teste de F2 dentro de cada A

  if (interacao_sig) {
    # --- F1 dentro de cada nível de F2 (letras minúsculas — linha) ---
    letras_min <- matrix("", nrow = k2, ncol = k1,
                         dimnames = list(niveis_f2, niveis_f1))
    desdf1_list <- list()

    for (b in niveis_f2) {
      idx_b  <- f2_vec == b
      y_b    <- y[idx_b]
      f1_b   <- droplevels(f1_vec[idx_b])

      # F-teste desdobrado (usando residuo global)
      medias_b <- tapply(y_b, f1_b, mean)
      ns_b     <- as.vector(table(f1_b))
      soma_sq  <- sum(ns_b * (medias_b - mean(y_b))^2)
      gl_desf  <- length(levels(f1_b)) - 1
      f_desf   <- if (gl_desf > 0) (soma_sq / gl_desf) / qm_res else NA
      p_desf   <- if (!is.na(f_desf)) 1 - stats::pf(f_desf, gl_desf, gl_res) else NA

      desdf1_list[[b]] <- data.frame(
        Nivel  = b,
        GL     = gl_desf,
        F_calc = f_desf,
        p_valor = p_desf,
        stringsAsFactors = FALSE
      )

      # Letras de comparação (apenas se F for significativo)
      if (!is.na(p_desf) && p_desf < alpha) {
        letras_lin <- .comparar_dentro(y_b, f1_b, qm_res, gl_res,
                                       teste = teste, confianca = confianca,
                                       upper = FALSE)
      } else {
        letras_lin <- stats::setNames(rep("a", length(levels(f1_b))), levels(f1_b))
      }
      letras_min[b, levels(f1_b)] <- letras_lin
    }
    desdobramento_f1_em_f2 <- do.call(rbind, desdf1_list)

    # --- F2 dentro de cada nível de F1 (letras maiúsculas — coluna) ---
    letras_mai <- matrix("", nrow = k2, ncol = k1,
                         dimnames = list(niveis_f2, niveis_f1))
    desdf2_list <- list()

    for (a in niveis_f1) {
      idx_a  <- f1_vec == a
      y_a    <- y[idx_a]
      f2_a   <- droplevels(f2_vec[idx_a])

      medias_a <- tapply(y_a, f2_a, mean)
      ns_a     <- as.vector(table(f2_a))
      soma_sq  <- sum(ns_a * (medias_a - mean(y_a))^2)
      gl_desf  <- length(levels(f2_a)) - 1
      f_desf   <- if (gl_desf > 0) (soma_sq / gl_desf) / qm_res else NA
      p_desf   <- if (!is.na(f_desf)) 1 - stats::pf(f_desf, gl_desf, gl_res) else NA

      desdf2_list[[a]] <- data.frame(
        Nivel  = a,
        GL     = gl_desf,
        F_calc = f_desf,
        p_valor = p_desf,
        stringsAsFactors = FALSE
      )

      if (!is.na(p_desf) && p_desf < alpha) {
        letras_col <- .comparar_dentro(y_a, f2_a, qm_res, gl_res,
                                       teste = teste, confianca = confianca,
                                       upper = TRUE)
      } else {
        letras_col <- stats::setNames(rep("A", length(levels(f2_a))), levels(f2_a))
      }
      letras_mai[levels(f2_a), a] <- letras_col
    }
    desdobramento_f2_em_f1 <- do.call(rbind, desdf2_list)

  } else {
    # Interação NÃO significativa — testa efeitos principais
    if (!is.na(p_f1) && p_f1 < alpha) {
      letras_f1 <- .comparar_dentro(y, f1_vec, qm_res, gl_res,
                                    teste = teste, confianca = confianca,
                                    upper = FALSE)
    } else {
      letras_f1 <- stats::setNames(rep("a", k1), niveis_f1)
    }
    if (!is.na(p_f2) && p_f2 < alpha) {
      letras_f2 <- .comparar_dentro(y, f2_vec, qm_res, gl_res,
                                    teste = teste, confianca = confianca,
                                    upper = FALSE)
    } else {
      letras_f2 <- stats::setNames(rep("a", k2), niveis_f2)
    }
  }

  # ─── OBJETO DE RETORNO ────────────────────────────────────────────────────
  resultado <- structure(
    list(
      nomes = list(
        resposta = nm_resp,
        fator1   = nm_f1,
        fator2   = nm_f2,
        bloco    = nm_bloco,
        interacao = nm_int
      ),
      delineamento       = delin_efetivo,
      niveis_f1          = niveis_f1,
      niveis_f2          = niveis_f2,
      n_total            = n_total,
      teste              = teste,
      confianca          = confianca,
      interacao_sig      = interacao_sig,
      anova              = tab_anova_out,
      parametros         = list(
        media_geral = media_geral,
        cv_exp      = cv_exp,
        r2          = r2
      ),
      medias_comb        = medias_comb,
      dp_comb            = dp_comb,
      n_comb             = n_comb,
      medias_f1          = medias_f1,
      medias_f2          = medias_f2,
      letras_min         = letras_min,
      letras_mai         = letras_mai,
      letras_f1          = letras_f1,
      letras_f2          = letras_f2,
      desdobramento_f1_em_f2 = desdobramento_f1_em_f2,
      desdobramento_f2_em_f1 = desdobramento_f2_em_f1,
      modelo_aov         = mod_aov,
      dados_int          = df_mod,
      grafico            = NULL
    ),
    class = "estatR_fatorial"
  )

  # Gráfico automático
  g_fat <- tryCatch(grafico_fatorial(resultado), error = function(e) NULL)
  resultado$grafico <- g_fat

  # Análise de resíduos
  if (isTRUE(residuos)) {
    if (exists("analise_residual", mode = "function")) {
      cat("\n")
      analise_residual(mod_aov)
    }
  }

  if (isTRUE(grafico) && !is.null(g_fat)) {
    print(g_fat)
  }

  resultado
}

#' @export
fatorial <- analise_fatorial_duplo

#' @export
fat <- analise_fatorial_duplo


# ─────────────────────────────────────────────────────────────────────────────
# MÉTODO PRINT
# ─────────────────────────────────────────────────────────────────────────────

#' @export
print.estatR_fatorial <- function(x, ...) {

  alpha     <- 1 - x$confianca
  delin_str <- if (x$delineamento == "dbc") "DELINEAMENTO EM BLOCOS CASUALIZADOS (DBC)"
               else "DELINEAMENTO INTEIRAMENTE CASUALIZADO (DIC)"

  .print_titulo(sprintf("ANÁLISE FATORIAL DUPLO (%s) — %s",
                        x$nomes$interacao, delin_str))

  cat(sprintf("  Variável Resposta:  %s\n", x$nomes$resposta))
  cat(sprintf("  Fator 1:            %s (%d níveis: %s)\n",
              x$nomes$fator1, length(x$niveis_f1),
              paste(x$niveis_f1, collapse = ", ")))
  cat(sprintf("  Fator 2:            %s (%d níveis: %s)\n",
              x$nomes$fator2, length(x$niveis_f2),
              paste(x$niveis_f2, collapse = ", ")))
  if (!is.null(x$nomes$bloco)) {
    cat(sprintf("  Blocos:             %s\n", x$nomes$bloco))
  }
  cat(sprintf("  Total de Obs.:      %d\n\n", x$n_total))

  # 1. QUADRO DE ANOVA
  .print_topico("QUADRO DA ANÁLISE DE VARIÂNCIA (ANOVA)")
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  Fonte de Variação               GL        SQ           QM         F calc    p-valor \n")
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  for (i in 1:nrow(x$anova)) {
    f_str <- x$anova$Fonte[i]
    gl_s  <- sprintf("%2d", x$anova$GL[i])
    sq_s  <- sprintf("%10s", .fmt_num_dic(x$anova$SQ[i], 2))
    qm_s  <- if (!is.na(x$anova$QM[i]))   sprintf("%10s", .fmt_num_dic(x$anova$QM[i], 2))   else "          "
    fc_s  <- if (!is.na(x$anova$F_calc[i])) sprintf("%8s", .fmt_num_dic(x$anova$F_calc[i], 2)) else "        "
    pv_s  <- sprintf("%11s", .fmt_p_dic(x$anova$p_valor[i]))
    cat(sprintf("  %-30s %4s   %10s   %10s   %8s %11s\n", f_str, gl_s, sq_s, qm_s, fc_s, pv_s))
  }
  cat("  ────────────────────────────────────────────────────────────────────────────── \n")
  cat("  * Significância: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")

  # 2. PARÂMETROS
  .print_topico("PARÂMETROS EXPERIMENTAIS")
  cat("  ────────────────────────────────────────────────────────── \n")
  cat(sprintf("  Média Geral:             %s\n", .fmt_num_dic(x$parametros$media_geral, 2)))
  cat(sprintf("  CV Experimental (%%):     %s%%\n", .fmt_num_dic(x$parametros$cv_exp, 2)))
  if (identical(x$tipo_A, "qualitativo") && identical(x$tipo_B, "qualitativo")) {
    cat(sprintf("  Eta² (Tamanho de Efeito): %s%% (%s)\n",
                .fmt_num_dic(x$parametros$r2 * 100, 2),
                .fmt_num_dic(x$parametros$r2, 3)))
  } else {
    cat(sprintf("  R² do Modelo:            %s\n", .fmt_num_dic(x$parametros$r2, 3)))
  }
  cat("  ────────────────────────────────────────────────────────── \n")
  cat("  * CV experimental: < 10% (Baixo / Alta precisão), 10-20% (Médio), 20-30% (Alto), > 30% (Muito alto - Pimentel-Gomes).\n\n")

  nm_teste_exib <- switch(x$teste,
    tukey       = "Tukey",
    scott_knott = "Scott-Knott",
    duncan      = "Duncan",
    snk         = "SNK",
    scheffe     = "Scheffé",
    bonferroni  = "Bonferroni"
  )

  if (x$interacao_sig) {
    # ── INTERAÇÃO SIGNIFICATIVA ──────────────────────────────────────────────
    .print_topico(sprintf("DESDOBRAMENTO DE %s DENTRO DE %s",
                          toupper(x$nomes$fator1), toupper(x$nomes$fator2)))
    cat("  ─────────────────────────────────────────────────────── \n")
    cat(sprintf("  %-20s  GL    F calc    p-valor\n", paste0(x$nomes$fator2, " / Nível")))
    cat("  ─────────────────────────────────────────────────────── \n")
    for (i in 1:nrow(x$desdobramento_f1_em_f2)) {
      d <- x$desdobramento_f1_em_f2[i, ]
      cat(sprintf("  %-20s  %2d  %8s  %11s\n",
                  d$Nivel, d$GL,
                  .fmt_num_dic(d$F_calc, 2),
                  .fmt_p_dic(d$p_valor)))
    }
    cat("  ─────────────────────────────────────────────────────── \n")
    cat("  * Significância: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")

    .print_topico(sprintf("DESDOBRAMENTO DE %s DENTRO DE %s",
                          toupper(x$nomes$fator2), toupper(x$nomes$fator1)))
    cat("  ─────────────────────────────────────────────────────── \n")
    cat(sprintf("  %-20s  GL    F calc    p-valor\n", paste0(x$nomes$fator1, " / Nível")))
    cat("  ─────────────────────────────────────────────────────── \n")
    for (i in 1:nrow(x$desdobramento_f2_em_f1)) {
      d <- x$desdobramento_f2_em_f1[i, ]
      cat(sprintf("  %-20s  %2d  %8s  %11s\n",
                  d$Nivel, d$GL,
                  .fmt_num_dic(d$F_calc, 2),
                  .fmt_p_dic(d$p_valor)))
    }
    cat("  ─────────────────────────────────────────────────────── \n")
    cat("  * Significância: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")

    # Tabela matricial de médias com letras combinadas
    .print_topico(sprintf("MÉDIAS POR COMBINAÇÃO — %s × %s (Teste: %s)",
                          x$nomes$fator1, x$nomes$fator2, nm_teste_exib))

    niveis_f1 <- x$niveis_f1
    niveis_f2 <- x$niveis_f2
    k1 <- length(niveis_f1)
    k2 <- length(niveis_f2)

    # Determina larguras das colunas
    w_row <- max(nchar(niveis_f2), nchar(x$nomes$fator2), 8) + 2
    # Formato de célula: "00,00 ab AB" -> ~14 chars
    w_col <- max(nchar(niveis_f1), 14) + 2

    sep_total <- w_row + k1 * w_col + 4
    sep_line  <- paste0("  ", paste(rep("\u2500", sep_total), collapse = ""), "\n")

    # Cabeçalho — Nome Fator1
    header_f1 <- sprintf("  %-*s", w_row, x$nomes$fator1)
    for (nv in niveis_f1) {
      header_f1 <- paste0(header_f1, sprintf(" %-*s", w_col - 1, nv))
    }
    cat(sep_line)
    cat(sprintf("  %-*s", w_row, x$nomes$fator2))
    for (nv in niveis_f1) {
      cat(sprintf(" %-*s", w_col - 1, nv))
    }
    cat("\n")
    cat(sep_line)

    for (b in niveis_f2) {
      cat(sprintf("  %-*s", w_row, b))
      for (a in niveis_f1) {
        med_v  <- x$medias_comb[a, b]
        l_min  <- if (!is.null(x$letras_min)) x$letras_min[b, a] else ""
        l_mai  <- if (!is.null(x$letras_mai)) x$letras_mai[b, a] else ""
        celula <- sprintf("%s %s%s", .fmt_num_dic(med_v, 2), l_min, l_mai)
        cat(sprintf(" %-*s", w_col - 1, celula))
      }
      cat("\n")
    }
    cat(sep_line)
    cat(sprintf("  * a, b - Médias seguidas pela mesma letra minúscula na linha não diferem\n"))
    cat(sprintf("    entre si pelo teste de %s ao nível de %g%% de significância.\n",
                nm_teste_exib, alpha * 100))
    cat(sprintf("  * A, B - Médias seguidas pela mesma letra maiúscula na coluna não diferem\n"))
    cat(sprintf("    entre si pelo teste de %s ao nível de %g%% de significância.\n\n",
                nm_teste_exib, alpha * 100))

  } else {
    # ── INTERAÇÃO NÃO SIGNIFICATIVA ─────────────────────────────────────────
    .print_topico(sprintf("EFEITO PRINCIPAL: %s", toupper(x$nomes$fator1)))
    if (!is.na(x$anova$p_valor[which(x$anova$Fonte == x$nomes$fator1)]) &&
        x$anova$p_valor[which(x$anova$Fonte == x$nomes$fator1)] < (1 - x$confianca)) {
      cat("  ────────────────────────────────────────────────── \n")
      cat(sprintf("  %-20s  %s\n", "Nível", "Média   Grupo"))
      cat("  ────────────────────────────────────────────────── \n")
      med_ord_f1 <- sort(x$medias_f1, decreasing = TRUE)
      for (nv in names(med_ord_f1)) {
        letra <- if (!is.null(x$letras_f1) && nv %in% names(x$letras_f1)) x$letras_f1[nv] else ""
        cat(sprintf("  %-20s  %s    %s\n", nv, .fmt_num_dic(med_ord_f1[nv], 2), letra))
      }
      cat("  ────────────────────────────────────────────────── \n")
      cat(sprintf("  * Médias seguidas pela mesma letra não diferem pelo teste de %s (5%%).\n\n",
                  nm_teste_exib))
    } else {
      cat(sprintf("  Efeito de %s não significativo (p = %s).\n\n",
                  x$nomes$fator1,
                  .fmt_num_dic(x$anova$p_valor[which(x$anova$Fonte == x$nomes$fator1)], 3)))
    }

    .print_topico(sprintf("EFEITO PRINCIPAL: %s", toupper(x$nomes$fator2)))
    if (!is.na(x$anova$p_valor[which(x$anova$Fonte == x$nomes$fator2)]) &&
        x$anova$p_valor[which(x$anova$Fonte == x$nomes$fator2)] < (1 - x$confianca)) {
      cat("  ────────────────────────────────────────────────── \n")
      cat(sprintf("  %-20s  %s\n", "Nível", "Média   Grupo"))
      cat("  ────────────────────────────────────────────────── \n")
      med_ord_f2 <- sort(x$medias_f2, decreasing = TRUE)
      for (nv in names(med_ord_f2)) {
        letra <- if (!is.null(x$letras_f2) && nv %in% names(x$letras_f2)) x$letras_f2[nv] else ""
        cat(sprintf("  %-20s  %s    %s\n", nv, .fmt_num_dic(med_ord_f2[nv], 2), letra))
      }
      cat("  ────────────────────────────────────────────────── \n")
      cat(sprintf("  * Médias seguidas pela mesma letra não diferem pelo teste de %s (5%%).\n\n",
                  nm_teste_exib))
    } else {
      cat(sprintf("  Efeito de %s não significativo (p = %s).\n\n",
                  x$nomes$fator2,
                  .fmt_num_dic(x$anova$p_valor[which(x$anova$Fonte == x$nomes$fator2)], 3)))
    }
  }

  # 4. INTERPRETAÇÃO
  .print_topico("INTERPRETAÇÃO")

  p_f1_val  <- x$anova$p_valor[which(x$anova$Fonte == x$nomes$fator1)]
  p_f2_val  <- x$anova$p_valor[which(x$anova$Fonte == x$nomes$fator2)]
  f_int_val <- x$anova$F_calc[which(x$anova$Fonte == x$nomes$interacao)]
  p_int_val <- x$anova$p_valor[which(x$anova$Fonte == x$nomes$interacao)]
  gl_int_v  <- x$anova$GL[which(x$anova$Fonte == x$nomes$interacao)]
  gl_res_v  <- x$anova$GL[which(x$anova$Fonte == "Resíduo")]

  if (x$interacao_sig) {
    cat(sprintf("  A interação %s foi estatisticamente significativa (F(%d, %d) = %s; p %s),\n  indicando que o efeito de %s sobre %s depende do nível de %s.\n\n",
                x$nomes$interacao,
                gl_int_v, gl_res_v,
                .fmt_num_dic(f_int_val, 2),
                if (!is.na(p_int_val) && p_int_val < 0.001) "< 0,001" else paste0("= ", .fmt_num_dic(p_int_val, 3)),
                x$nomes$fator1, x$nomes$resposta, x$nomes$fator2))
    cat(sprintf("  O desdobramento revelou os efeitos simples de cada fator dentro dos níveis do outro,\n  sendo as comparações múltiplas realizadas pelo teste de %s a %g%% de significância.\n\n",
                nm_teste_exib, (1 - x$confianca) * 100))
  } else {
    cat(sprintf("  A interação %s não foi significativa (F(%d, %d) = %s; p = %s),\n  indicando que os fatores atuam de forma independente sobre %s.\n\n",
                x$nomes$interacao,
                gl_int_v, gl_res_v,
                .fmt_num_dic(f_int_val, 2),
                .fmt_num_dic(p_int_val, 3),
                x$nomes$resposta))

    f1_sig <- !is.na(p_f1_val) && p_f1_val < (1 - x$confianca)
    f2_sig <- !is.na(p_f2_val) && p_f2_val < (1 - x$confianca)
    if (f1_sig && f2_sig) {
      cat(sprintf("  Ambos os fatores (%s e %s) apresentaram efeito significativo sobre %s.\n  As comparações múltiplas foram realizadas pelo teste de %s a %g%% de significância.\n\n",
                  x$nomes$fator1, x$nomes$fator2, x$nomes$resposta,
                  nm_teste_exib, (1 - x$confianca) * 100))
    } else if (f1_sig) {
      cat(sprintf("  Apenas %s apresentou efeito significativo sobre %s.\n\n",
                  x$nomes$fator1, x$nomes$resposta))
    } else if (f2_sig) {
      cat(sprintf("  Apenas %s apresentou efeito significativo sobre %s.\n\n",
                  x$nomes$fator2, x$nomes$resposta))
    } else {
      cat(sprintf("  Nenhum dos fatores avaliados apresentou efeito significativo sobre %s.\n\n",
                  x$nomes$resposta))
    }
  }

  cat(sprintf("  O Coeficiente de Variação Experimental (CV = %s%%) indica precisão experimental %s\n  segundo os critérios de Pimentel-Gomes.\n\n",
              .fmt_num_dic(x$parametros$cv_exp, 2),
              if (x$parametros$cv_exp < 10) "alta (baixo CV)"
              else if (x$parametros$cv_exp < 20) "média"
              else if (x$parametros$cv_exp < 30) "baixa (alto CV)"
              else "muito baixa (CV muito alto)"))

  .print_rodape()
  invisible(x)
}


# ─────────────────────────────────────────────────────────────────────────────
# GRÁFICO DE INTERAÇÃO (PERFIL)
# ─────────────────────────────────────────────────────────────────────────────

#' @title Gráfico de Perfil de Interação — Esquema Fatorial
#' @description Gera o gráfico de perfil de interação para um fatorial duplo,
#'   com linhas por nível de Fator2, barras de erro (IC) e anotação das letras
#'   de comparação múltipla (quando a interação é significativa).
#' @param modelo Objeto retornado por \code{analise_fatorial()} ou \code{fatorial()}.
#' @param ... Argumentos adicionais (não utilizados).
#' @return Objeto \code{ggplot}.
#' @export
#' @aliases plot_fatorial
grafico_fatorial <- function(modelo, ...) {
  if (!inherits(modelo, "estatR_fatorial")) {
    stop("O argumento 'modelo' deve ser um objeto retornado por analise_fatorial().", call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("O pacote 'ggplot2' é necessário para gerar gráficos.")
    return(invisible(NULL))
  }

  niveis_f1 <- modelo$niveis_f1
  niveis_f2 <- modelo$niveis_f2
  alpha_val  <- 1 - modelo$confianca
  qm_res     <- summary(modelo$modelo_aov)[[1]]["Residuals", "Mean Sq"]
  gl_res     <- summary(modelo$modelo_aov)[[1]]["Residuals", "Df"]

  nm_teste_exib <- switch(modelo$teste,
    tukey       = "Tukey",
    scott_knott = "Scott-Knott",
    duncan      = "Duncan",
    snk         = "SNK",
    scheffe     = "Scheffé",
    bonferroni  = "Bonferroni",
    "Tukey"
  )

  # Montar data.frame para o gráfico
  rows <- list()
  for (b in niveis_f2) {
    for (a in niveis_f1) {
      med_v <- modelo$medias_comb[a, b]
      n_v   <- modelo$n_comb[a, b]
      ep_v  <- sqrt(qm_res / n_v)
      t_crit <- stats::qt(1 - alpha_val / 2, gl_res)
      ic_lo  <- med_v - t_crit * ep_v
      ic_hi  <- med_v + t_crit * ep_v

      if (modelo$interacao_sig) {
        l_min <- if (!is.null(modelo$letras_min)) modelo$letras_min[b, a] else ""
        l_mai <- if (!is.null(modelo$letras_mai)) modelo$letras_mai[b, a] else ""
        letra <- paste0(l_min, l_mai)
      } else {
        letra <- ""
      }

      rows[[length(rows) + 1]] <- data.frame(
        Fator1 = a, Fator2 = b,
        Media = med_v, IC_lo = ic_lo, IC_hi = ic_hi,
        Letra = letra,
        stringsAsFactors = FALSE
      )
    }
  }
  df_plot <- do.call(rbind, rows)
  df_plot$Fator1 <- factor(df_plot$Fator1, levels = niveis_f1)
  df_plot$Fator2 <- factor(df_plot$Fator2, levels = niveis_f2)

  # Paleta de cores
  pal <- if (exists(".paleta_estatR", mode = "function")) {
    .paleta_estatR()
  } else {
    c("#2166AC", "#D6604D", "#4DAC26", "#B2182B", "#762A83",
      "#E08214", "#1B7837", "#542788", "#8073AC", "#F46D43")
  }
  cores <- pal[seq_along(niveis_f2)]

  # Tema
  tema_fn <- if (exists("tema_estatR", mode = "function")) {
    tema_estatR(estilo = 2, grade = "dupla")
  } else if (exists("meu_tema", mode = "function")) {
    meu_tema(estilo = 2, grade = "dupla")
  } else {
    ggplot2::theme_minimal()
  }

  g <- ggplot2::ggplot(df_plot,
         ggplot2::aes(x = Fator1, y = Media, color = Fator2,
                      group = Fator2, linetype = Fator2)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = IC_lo, ymax = IC_hi),
                           width = 0.15, linewidth = 0.6) +
    ggplot2::scale_color_manual(values = cores) +
    ggplot2::labs(
      title    = sprintf("Interação: %s × %s", modelo$nomes$fator1, modelo$nomes$fator2),
      subtitle = sprintf("Variável: %s  |  Teste: %s  |  Barras de erro: IC %g%%",
                         modelo$nomes$resposta, nm_teste_exib, modelo$confianca * 100),
      x     = modelo$nomes$fator1,
      y     = modelo$nomes$resposta,
      color = modelo$nomes$fator2,
      linetype = modelo$nomes$fator2
    ) +
    tema_fn

  # Adiciona letras se interação significativa
  if (modelo$interacao_sig && any(df_plot$Letra != "")) {
    g <- g + ggplot2::geom_text(
      ggplot2::aes(label = Letra, y = IC_hi),
      vjust = -0.5, size = 3.2,
      show.legend = FALSE,
      position = ggplot2::position_dodge(width = 0)
    )
  }

  g
}

#' @export
plot_fatorial <- grafico_fatorial


# ─────────────────────────────────────────────────────────────────────────────
# PAINEL GRÁFICO FATORIAL
# ─────────────────────────────────────────────────────────────────────────────

#' @title Painel Gráfico Completo do Fatorial
#' @description Exibe o gráfico de perfil de interação do esquema fatorial duplo.
#'   Quando a interação não é significativa, pode incluir gráficos de efeito
#'   principal lado a lado usando patchwork.
#' @param modelo Objeto retornado por \code{analise_fatorial_duplo()} ou \code{fatorial()}.
#' @param ... Argumentos adicionais passados para as funções gráficas individuais.
#' @return Objeto \code{ggplot} ou composição de gráficos.
#' @export
#' @aliases plot_fat
grafico_fat <- function(modelo, ...) {
  if (!inherits(modelo, "estatR_fatorial")) {
    stop("O argumento 'modelo' deve ser um objeto retornado por analise_fatorial_duplo().", call. = FALSE)
  }
  grafico_fatorial(modelo, ...)
}

#' @export
plot_fat <- grafico_fat

