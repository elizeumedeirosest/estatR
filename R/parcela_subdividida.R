# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: ANÁLISE DE PARCELAS SUBDIVIDIDAS (SPLIT-PLOT) — DIC e DBC
# ─────────────────────────────────────────────────────────────────────────────

#' @title Análise de Experimentos em Parcelas Subdivididas (Split-Plot)
#' @description Realiza a análise estatística completa de experimentos instalados
#'   no esquema de Parcelas Subdivididas (Split-Plot) em Delineamento Inteiramente
#'   Casualizado (DIC) ou em Blocos Casualizados (DBC).
#'   Calcula o quadro de ANOVA com dois resíduos (Erro a para parcelas e Erro b para subparcelas),
#'   dois coeficientes de variação (CVa e CVb), testes de comparações múltiplas (Tukey,
#'   Scott-Knott, Duncan, SNK, Scheffé ou Bonferroni) e desdobramento automático
#'   quando a interação Parcela × Subparcela for significativa (utilizando combinações
#'   lineares de erros conforme Satterthwaite/Banzatto & Kronka).
#'
#' @param resposta Variável resposta numérica contínua.
#' @param parcela Fator aplicado nas parcelas principais (fator primário).
#' @param subparcela Fator aplicado nas subparcelas (fator secundário).
#' @param bloco Fator de blocos (opcional). Se informado, o experimento é ajustado
#'   como DBC em parcelas subdivididas; se omitido, como DIC em parcelas subdivididas.
#' @param repeticao Vetor ou fator que identifica a repetição da parcela principal
#'   no caso de DIC (opcional). Se omitido em DIC, assume repetições balanceadas.
#' @param dados Data frame contendo os dados (opcional se vetores forem fornecidos diretamente).
#' @param teste Teste de médias a ser aplicado: \code{"tukey"} (padrão), \code{"scott_knott"},
#'   \code{"duncan"}, \code{"snk"}, \code{"scheffe"} ou \code{"bonferroni"}.
#' @param confianca Nível de confiança para os testes e intervalos (padrão: 0.95).
#' @param grafico Lógico. Se TRUE, gera gráfico de interação ou efeitos principais.
#'
#' @return Objeto da classe \code{"estatR_subdividida"} contendo a ANOVA, coeficientes de variação,
#'   tabelas de médias e desdobramentos.
#' @export
analise_parcela_subdividida <- function(
    resposta, parcela, subparcela,
    bloco       = NULL,
    repeticao   = NULL,
    dados       = NULL,
    teste       = c("tukey", "scott_knott", "duncan", "snk", "scheffe", "bonferroni"),
    confianca   = 0.95,
    grafico     = TRUE
) {
  teste <- match.arg(teste)

  sub_resp  <- substitute(resposta)
  sub_parc  <- substitute(parcela)
  sub_subp  <- substitute(subparcela)
  sub_bloco <- substitute(bloco)
  sub_rep   <- substitute(repeticao)

  info_resp <- .extrair_vetor_exp(sub_resp, dados, "resposta")
  info_parc <- .extrair_vetor_exp(sub_parc, dados, "parcela")
  info_subp <- .extrair_vetor_exp(sub_subp, dados, "subparcela")

  y_raw    <- info_resp$valores
  parc_raw <- info_parc$valores
  subp_raw <- info_subp$valores
  nm_resp  <- info_resp$nome
  nm_parc  <- info_parc$nome
  nm_subp  <- info_subp$nome

  # Identificação do bloco
  bloco_raw <- NULL
  nm_bloco  <- NULL
  has_bloco <- !is.null(sub_bloco) && !identical(sub_bloco, as.name("NULL"))
  if (has_bloco) {
    info_bloco <- .extrair_vetor_exp(sub_bloco, dados, "bloco")
    bloco_raw  <- info_bloco$valores
    nm_bloco   <- info_bloco$nome
  }

  is_dbc <- has_bloco && !is.null(bloco_raw)

  # Repetição no DIC
  rep_raw <- NULL
  has_rep <- !is.null(sub_rep) && !identical(sub_rep, as.name("NULL"))
  if (has_rep) {
    info_rep <- .extrair_vetor_exp(sub_rep, dados, "repeticao")
    rep_raw  <- info_rep$valores
  }

  # Filtro de NAs
  if (is_dbc) {
    val_ok <- !is.na(y_raw) & !is.na(parc_raw) & !is.na(subp_raw) & !is.na(bloco_raw)
  } else {
    val_ok <- !is.na(y_raw) & !is.na(parc_raw) & !is.na(subp_raw)
    if (has_rep) val_ok <- val_ok & !is.na(rep_raw)
  }

  y    <- as.numeric(y_raw[val_ok])
  parc <- as.factor(parc_raw[val_ok])
  subp <- as.factor(subp_raw[val_ok])
  bloco_v <- if (is_dbc) as.factor(bloco_raw[val_ok]) else NULL
  rep_v   <- if (has_rep) as.factor(rep_raw[val_ok]) else NULL

  n_total <- length(y)
  a_niv   <- levels(parc)
  b_niv   <- levels(subp)
  I       <- length(a_niv) # niveis da parcela
  J       <- length(b_niv) # niveis da subparcela

  if (I < 2) stop("O fator da parcela principal deve possuir ao menos 2 n\u00edveis.", call. = FALSE)
  if (J < 2) stop("O fator da subparcela deve possuir ao menos 2 n\u00edveis.", call. = FALSE)

  # ─── CÁLCULO DA ANOVA DE PARCELAS SUBDIVIDIDAS ─────────────────────────────
  # Em DBC: Y ~ Bloco + Parcela + Bloco:Parcela (Erro a) + Subparcela + Parcela:Subparcela + Erro b
  # Em DIC: Y ~ Parcela + Rep/Parcela (Erro a) + Subparcela + Parcela:Subparcela + Erro b

  if (is_dbc) {
    K <- length(levels(bloco_v)) # numero de blocos (repeticoes)
    df_mod <- data.frame(Y = y, Bloco = bloco_v, Parcela = parc, Subparcela = subp)

    # Ajuste por aov com Error stratum ou manual para precisão absoluta:
    sq_tot <- sum((y - mean(y))^2)
    gl_tot <- n_total - 1

    # Médias e totais
    t_bloco <- tapply(y, bloco_v, sum)
    t_parc  <- tapply(y, parc, sum)
    t_subp  <- tapply(y, subp, sum)
    t_parc_bloco <- tapply(y, list(parc, bloco_v), sum)
    t_parc_subp  <- tapply(y, list(parc, subp), sum)

    C <- (sum(y)^2) / n_total

    sq_bloco <- sum(t_bloco^2) / (I * J) - C
    gl_bloco <- K - 1
    qm_bloco <- sq_bloco / gl_bloco

    sq_parc <- sum(t_parc^2) / (J * K) - C
    gl_parc <- I - 1
    qm_parc <- sq_parc / gl_parc

    # Erro (a) = Bloco:Parcela
    sq_parc_bloco_tot <- sum(t_parc_bloco^2, na.rm = TRUE) / J - C
    sq_ea <- sq_parc_bloco_tot - sq_bloco - sq_parc
    gl_ea <- (I - 1) * (K - 1)
    qm_ea <- sq_ea / gl_ea

    f_bloco <- qm_bloco / qm_ea
    p_bloco <- 1 - stats::pf(f_bloco, gl_bloco, gl_ea)

    f_parc  <- qm_parc / qm_ea
    p_parc  <- 1 - stats::pf(f_parc, gl_parc, gl_ea)

    # Subparcela e Interacao
    sq_subp <- sum(t_subp^2) / (I * K) - C
    gl_subp <- J - 1
    qm_subp <- sq_subp / gl_subp

    sq_inter_tot <- sum(t_parc_subp^2, na.rm = TRUE) / K - C
    sq_int <- sq_inter_tot - sq_parc - sq_subp
    gl_int <- (I - 1) * (J - 1)
    qm_int <- sq_int / gl_int

    # Erro (b) = Residuo global menos componentes de parcela
    sq_eb <- sq_tot - sq_bloco - sq_parc - sq_ea - sq_subp - sq_int
    gl_eb <- I * (K - 1) * (J - 1)
    qm_eb <- sq_eb / gl_eb

    f_subp <- qm_subp / qm_eb
    p_subp <- 1 - stats::pf(f_subp, gl_subp, gl_eb)

    f_int  <- qm_int / qm_eb
    p_int  <- 1 - stats::pf(f_int, gl_int, gl_eb)

    tab_anova <- data.frame(
      Fonte = c("Bloco", nm_parc, "Erro (a)", nm_subp, paste0(nm_parc, " \u00d7 ", nm_subp), "Erro (b)", "Total"),
      GL    = c(gl_bloco, gl_parc, gl_ea, gl_subp, gl_int, gl_eb, gl_tot),
      SQ    = c(sq_bloco, sq_parc, sq_ea, sq_subp, sq_int, sq_eb, sq_tot),
      QM    = c(qm_bloco, qm_parc, qm_ea, qm_subp, qm_int, qm_eb, NA),
      F_calc= c(f_bloco, f_parc, NA, f_subp, f_int, NA, NA),
      p_valor=c(p_bloco, p_parc, NA, p_subp, p_int, NA, NA),
      stringsAsFactors = FALSE
    )

  } else {
    # DIC em Parcelas Subdivididas
    # Se rep_v não foi fornecido, gera repetição com base na ordem dos dados
    if (is.null(rep_v)) {
      # repete 1:r para cada combinacao
      r_est <- n_total / (I * J)
      rep_v <- factor(rep(seq_len(r_est), each = J, times = I))
    }
    K <- length(unique(rep_v)) # repeticoes

    C <- (sum(y)^2) / n_total
    sq_tot <- sum((y - mean(y))^2)
    gl_tot <- n_total - 1

    t_parc <- tapply(y, parc, sum)
    sq_parc <- sum(t_parc^2) / (J * K) - C
    gl_parc <- I - 1
    qm_parc <- sq_parc / gl_parc

    # Parcelas = Parcela + Rep/Parcela
    t_parc_rep <- tapply(y, list(parc, rep_v), sum)
    sq_parc_tot <- sum(t_parc_rep^2, na.rm = TRUE) / J - C
    sq_ea <- sq_parc_tot - sq_parc
    gl_ea <- I * (K - 1)
    qm_ea <- sq_ea / gl_ea

    f_parc <- qm_parc / qm_ea
    p_parc <- 1 - stats::pf(f_parc, gl_parc, gl_ea)

    t_subp <- tapply(y, subp, sum)
    sq_subp <- sum(t_subp^2) / (I * K) - C
    gl_subp <- J - 1
    qm_subp <- sq_subp / gl_subp

    t_parc_subp <- tapply(y, list(parc, subp), sum)
    sq_inter_tot <- sum(t_parc_subp^2, na.rm = TRUE) / K - C
    sq_int <- sq_inter_tot - sq_parc - sq_subp
    gl_int <- (I - 1) * (J - 1)
    qm_int <- sq_int / gl_int

    sq_eb <- sq_tot - sq_parc - sq_ea - sq_subp - sq_int
    gl_eb <- I * (K - 1) * (J - 1)
    qm_eb <- sq_eb / gl_eb

    f_subp <- qm_subp / qm_eb
    p_subp <- 1 - stats::pf(f_subp, gl_subp, gl_eb)

    f_int  <- qm_int / qm_eb
    p_int  <- 1 - stats::pf(f_int, gl_int, gl_eb)

    tab_anova <- data.frame(
      Fonte = c(nm_parc, "Erro (a)", nm_subp, paste0(nm_parc, " \u00d7 ", nm_subp), "Erro (b)", "Total"),
      GL    = c(gl_parc, gl_ea, gl_subp, gl_int, gl_eb, gl_tot),
      SQ    = c(sq_parc, sq_ea, sq_subp, sq_int, sq_eb, sq_tot),
      QM    = c(qm_parc, qm_ea, qm_subp, qm_int, qm_eb, NA),
      F_calc= c(f_parc, NA, f_subp, f_int, NA, NA),
      p_valor=c(p_parc, NA, p_subp, p_int, NA, NA),
      stringsAsFactors = FALSE
    )
  }

  media_geral <- mean(y)
  cv_a <- (sqrt(qm_ea) / media_geral) * 100
  cv_b <- (sqrt(qm_eb) / media_geral) * 100

  # ─── COMPARAÇÕES MÚLTIPLAS / DESDOBRAMENTO ───────────────────────────────
  # Matriz de médias combinadas
  mat_medias <- tapply(y, list(parc, subp), mean)
  mat_dp     <- tapply(y, list(parc, subp), stats::sd)

  # Se interacao for significativa (p < 0.05):
  # 1. Subparcela dentro de Parcela: usa Erro b (qm_eb, gl_eb)
  # 2. Parcela dentro de Subparcela: usa erro combinado (qm_comb, gl_comb de Satterthwaite)
  interacao_sig <- !is.na(p_int) && (p_int < 0.05)

  desdobramento_subp_em_parc <- list()
  desdobramento_parc_em_subp <- list()
  mat_letras_min <- matrix("", nrow = I, ncol = J, dimnames = list(a_niv, b_niv))
  mat_letras_mai <- matrix("", nrow = I, ncol = J, dimnames = list(a_niv, b_niv))

  # Letras de Subparcela dentro de cada Parcela (letras minúsculas na linha)
  for (i_idx in seq_along(a_niv)) {
    nv_a <- a_niv[i_idx]
    sub_df <- (parc == nv_a)
    y_sub  <- y[sub_df]
    b_sub  <- subp[sub_df]
    letras <- .comparar_dentro(y_sub, b_sub, qm_res_global = qm_eb, gl_res_global = gl_eb, teste = teste, confianca = confianca, upper = FALSE)
    mat_letras_min[nv_a, names(letras)] <- letras
  }

  # Letras de Parcela dentro de cada Subparcela (letras maiúsculas na coluna)
  # Erro combinado: qm_comb = (qm_ea + (J - 1) * qm_eb) / J
  # GL Satterthwaite: gl_comb = (qm_ea + (J-1)*qm_eb)^2 / [ (qm_ea^2 / gl_ea) + ((J-1)*qm_eb)^2 / gl_eb ]
  qm_comb <- (qm_ea + (J - 1) * qm_eb) / J
  gl_comb <- ((qm_ea + (J - 1) * qm_eb)^2) / ((qm_ea^2 / gl_ea) + (((J - 1) * qm_eb)^2 / gl_eb))

  for (j_idx in seq_along(b_niv)) {
    nv_b <- b_niv[j_idx]
    sub_df <- (subp == nv_b)
    y_sub  <- y[sub_df]
    a_sub  <- parc[sub_df]
    letras <- .comparar_dentro(y_sub, a_sub, qm_res_global = qm_comb, gl_res_global = gl_comb, teste = teste, confianca = confianca, upper = TRUE)
    mat_letras_mai[names(letras), nv_b] <- letras
  }

  # Médias independentes para efeitos principais se interacao NAO for significativa
  medias_parc <- tapply(y, parc, mean)
  medias_subp <- tapply(y, subp, mean)
  letras_parc <- .comparar_dentro(y, parc, qm_res_global = qm_ea, gl_res_global = gl_ea, teste = teste, confianca = confianca, upper = FALSE)
  letras_subp <- .comparar_dentro(y, subp, qm_res_global = qm_eb, gl_res_global = gl_eb, teste = teste, confianca = confianca, upper = FALSE)

  # Objeto com modelo aov para análise residual
  df_fit <- if (is_dbc) {
    data.frame(Y = y, Bloco = bloco_v, Parcela = parc, Subparcela = subp)
  } else {
    data.frame(Y = y, Parcela = parc, Subparcela = subp)
  }
  mod_aov <- if (is_dbc) {
    stats::aov(Y ~ Bloco + Parcela * Subparcela, data = df_fit)
  } else {
    stats::aov(Y ~ Parcela * Subparcela, data = df_fit)
  }

  resultado <- structure(
    list(
      anova         = tab_anova,
      delineamento  = if (is_dbc) "DBC (Parcelas Subdivididas)" else "DIC (Parcelas Subdivididas)",
      nomes         = list(resposta = nm_resp, parcela = nm_parc, subparcela = nm_subp, bloco = nm_bloco),
      niveis_parcela= a_niv,
      niveis_subp   = b_niv,
      n_total       = n_total,
      parametros    = list(media_geral = media_geral, cv_a = cv_a, cv_b = cv_b, qm_ea = qm_ea, qm_eb = qm_eb, qm_comb = qm_comb),
      interacao_sig = interacao_sig,
      medias_comb   = mat_medias,
      dp_comb       = mat_dp,
      letras_min    = mat_letras_min,
      letras_mai    = mat_letras_mai,
      medias_parc   = medias_parc,
      letras_parc   = letras_parc,
      medias_subp   = medias_subp,
      letras_subp   = letras_subp,
      teste         = teste,
      confianca     = confianca,
      modelo_aov    = mod_aov,
      dados         = df_fit
    ),
    class = "estatR_subdividida"
  )

  print(resultado)

  if (grafico) {
    tryCatch({
      # Gráfico de interação estilo Forest / Linhas de perfil com barras de erro
      df_g <- data.frame(
        Parcela    = rep(a_niv, each = J),
        Subparcela = rep(b_niv, times = I),
        Media      = as.vector(t(mat_medias)),
        DP         = as.vector(t(mat_dp)),
        Letras     = paste0(as.vector(t(mat_letras_min)), as.vector(t(mat_letras_mai)))
      )

      p <- ggplot2::ggplot(df_g, ggplot2::aes(x = Parcela, y = Media, color = Subparcela, group = Subparcela)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = Media - DP, ymax = Media + DP), width = 0.15, position = ggplot2::position_dodge(0.2)) +
        ggplot2::geom_line(linewidth = 0.9, position = ggplot2::position_dodge(0.2)) +
        ggplot2::geom_point(size = 3.2, position = ggplot2::position_dodge(0.2)) +
        ggplot2::geom_text(ggplot2::aes(label = Letras), vjust = -1.2, position = ggplot2::position_dodge(0.2), show.legend = FALSE, fontface = "bold") +
        ggplot2::labs(
          title = sprintf("Perfil de Intera\u00e7\u00e3o: %s \u00d7 %s", nm_parc, nm_subp),
          subtitle = sprintf("Letras min\u00fasculas comparam %s dentro de %s; Mai\u00fasculas comparam %s dentro de %s", nm_subp, nm_parc, nm_parc, nm_subp),
          x = nm_parc,
          y = nm_resp,
          color = nm_subp
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
print.estatR_subdividida <- function(x, ...) {
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo(sprintf("EXPERIMENTO EM PARCELAS SUBDIVIDIDAS (%s)", if (grepl("DBC", x$delineamento)) "DBC" else "DIC"))
  } else {
    cat(sprintf("\n── EXPERIMENTO EM PARCELAS SUBDIVIDIDAS (%s) ──\n", if (grepl("DBC", x$delineamento)) "DBC" else "DIC"))
  }

  cat(sprintf("  Vari\u00e1vel Resposta:    %s\n", x$nomes$resposta))
  cat(sprintf("  Parcela Principal:    %s (%d n\u00edveis: %s)\n", x$nomes$parcela, length(x$niveis_parcela), paste(x$niveis_parcela, collapse = ", ")))
  cat(sprintf("  Subparcela:           %s (%d n\u00edveis: %s)\n", x$nomes$subparcela, length(x$niveis_subp), paste(x$niveis_subp, collapse = ", ")))
  if (!is.null(x$nomes$bloco)) cat(sprintf("  Blocos:               %s\n", x$nomes$bloco))
  cat(sprintf("  Total de Observa\u00e7\u00f5es: %d\n\n", x$n_total))

  if (exists(".print_topico", mode = "function")) .print_topico("QUADRO DA AN\u00c1LISE DE VARI\u00c2NCIA (ANOVA)")

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
  cat(sprintf("  CV(a) - Parcela:         %.2f%%\n", x$parametros$cv_a))
  cat(sprintf("  CV(b) - Subparcela:      %.2f%%\n", x$parametros$cv_b))
  cat("  ────────────────────────────────────────────────────────── \n")
  cat("  * CV(a) avalia a precis\u00e3o nas parcelas principais; CV(b) nas subparcelas.\n\n")

  nm_teste_exib <- toupper(x$teste)

  if (x$interacao_sig) {
    if (exists(".print_topico", mode = "function")) {
      .print_topico(sprintf("DESDOBRAMENTO DA INTERA\u00c7\u00c3O (TESTE DE %s - \u03b1 = %.0f%%)", nm_teste_exib, (1 - x$confianca) * 100))
    }

    cat("  ────────────────────────────────────────────────────────────────────────────── \n")
    cat(sprintf("  Médias com as mesmas letras minúsculas na LINHA não diferem entre si (%s).\n", x$nomes$subparcela))
    cat(sprintf("  Médias com as mesmas letras maiúsculas na COLUNA não diferem entre si (%s).\n", x$nomes$parcela))
    cat("  ────────────────────────────────────────────────────────────────────────────── \n\n")

    b_niv <- x$niveis_subp
    a_niv <- x$niveis_parcela

    # Cabeçalho da tabela
    cat(sprintf("  %-20s", x$nomes$parcela))
    for (nb in b_niv) cat(sprintf("  %14s", nb))
    cat("\n  ────────────────────────────────────────────────────────────────────────────── \n")

    for (na in a_niv) {
      cat(sprintf("  %-20s", na))
      for (nb in b_niv) {
        val_m <- x$medias_comb[na, nb]
        l_min <- x$letras_min[na, nb]
        l_mai <- x$letras_mai[na, nb]
        txt <- sprintf("%.2f %s%s", val_m, l_min, l_mai)
        cat(sprintf("  %14s", txt))
      }
      cat("\n")
    }
    cat("  ────────────────────────────────────────────────────────────────────────────── \n\n")

  } else {
    if (exists(".print_topico", mode = "function")) {
      .print_topico(sprintf("EFEITOS PRINCIPAIS INDEPENDENTES (TESTE DE %s)", nm_teste_exib))
    }

    cat(sprintf("\n  FATOR PARCELA: %s\n", x$nomes$parcela))
    cat("  ──────────────────────────────────────────── \n")
    for (na in names(x$medias_parc)) {
      cat(sprintf("  %-20s  %10.2f  %s\n", na, x$medias_parc[na], x$letras_parc[na]))
    }

    cat(sprintf("\n  FATOR SUBPARCELA: %s\n", x$nomes$subparcela))
    cat("  ──────────────────────────────────────────── \n")
    for (nb in names(x$medias_subp)) {
      cat(sprintf("  %-20s  %10.2f  %s\n", nb, x$medias_subp[nb], x$letras_subp[nb]))
    }
    cat("  ──────────────────────────────────────────── \n\n")
  }

  invisible(x)
}
