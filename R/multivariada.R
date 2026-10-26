#' @title Análise de Componentes Principais (ACP / PCA)
#' @description Realiza Análise de Componentes Principais (ACP) com diagnóstico completo
#' de adequação da amostra (KMO e Bartlett), critérios de retenção de componentes
#' (Kaiser, Análise Paralela de Horn e Variância Acumulada), cargas fatoriais (loadings),
#' escores e interpretação automatizada.
#' @param dados Data frame ou matriz contendo as variáveis numéricas.
#' @param padronizar Lógico. Se TRUE (padrão), padroniza as variáveis (matriz de correlação).
#'   Se FALSE, utiliza a matriz de covariância.
#' @param n_comp Inteiro opcional. Número de componentes a exibir nas cargas e reter.
#'   Se NULL (padrão), retém automaticamente pelo Critério de Kaiser (autovalor > 1).
#' @param analise_paralela Lógico. Se TRUE (padrão), executa a Análise Paralela de Horn.
#' @param n_sim Inteiro. Número de simulações/reamostragens para a Análise Paralela (padrão: 100).
#' @param grafico Lógico. Se TRUE, exibe automaticamente o biplot ao final.
#' @return Objeto da classe \code{acp_estatR} contendo:
#'   \item{autovalores}{Data frame com Componente, Autovalor, Var_Expl_Pct e Var_Acum_Pct}
#'   \item{cargas}{Data frame com as cargas (loadings) dos componentes retidos}
#'   \item{escores}{Data frame com os escores dos componentes para cada observação}
#'   \item{kmo}{Lista com o KMO global, nível e KMO individual por variável}
#'   \item{bartlett}{Lista com estatística qui-quadrado, graus de liberdade e p-valor}
#'   \item{criterios}{Data frame com resumo dos critérios de retenção}
#'   \item{prcomp}{Objeto nativo da função stats::prcomp}
#' @export
acp <- function(dados,
                padronizar       = TRUE,
                n_comp           = NULL,
                analise_paralela = TRUE,
                n_sim            = 100,
                grafico          = TRUE) {

  # 1. Validação e seleção de variáveis numéricas
  if (!is.data.frame(dados) && !is.matrix(dados)) {
    stop("O argumento 'dados' deve ser um data frame ou matriz.")
  }

  dados_df <- as.data.frame(dados)
  cols_num <- vapply(dados_df, is.numeric, logical(1))
  if (sum(cols_num) < 2) {
    stop("São necessárias pelo menos 2 variáveis numéricas para a ACP.")
  }

  df_num <- dados_df[, cols_num, drop = FALSE]
  casos_completos <- stats::complete.cases(df_num)
  if (sum(casos_completos) < nrow(df_num)) {
    df_num <- df_num[casos_completos, , drop = FALSE]
    dados_df <- dados_df[casos_completos, , drop = FALSE]
  }

  n <- nrow(df_num)
  p <- ncol(df_num)
  var_nomes <- colnames(df_num)

  # 2. Diagnóstico de Adequação: KMO e Teste de Bartlett
  mat_cor <- stats::cor(df_num)

  # Bartlett
  det_r <- det(mat_cor)
  if (is.na(det_r) || det_r <= 0) {
    bart_chi2 <- NA; bart_df <- p * (p - 1) / 2; bart_p <- NA
  } else {
    bart_df   <- p * (p - 1) / 2
    bart_chi2 <- -( (n - 1) - (2 * p + 5) / 6 ) * log(det_r)
    bart_p    <- stats::pchisq(bart_chi2, df = bart_df, lower.tail = FALSE)
  }

  # KMO
  inv_r <- tryCatch(solve(mat_cor), error = function(e) NULL)
  if (is.null(inv_r)) {
    kmo_global <- NA
    kmo_ind <- rep(NA, p)
  } else {
    anti_img <- -inv_r / sqrt(outer(diag(inv_r), diag(inv_r)))
    diag(anti_img) <- 1

    r2 <- mat_cor^2
    a2 <- anti_img^2
    diag(r2) <- 0
    diag(a2) <- 0

    kmo_global <- sum(r2) / (sum(r2) + sum(a2))
    kmo_ind <- rowSums(r2) / (rowSums(r2) + rowSums(a2))
  }

  .classifica_kmo <- function(val) {
    if (is.na(val)) return("Indeterminado")
    if (val >= 0.90) return("Excelente")
    if (val >= 0.80) return("Meritório")
    if (val >= 0.70) return("Médio")
    if (val >= 0.60) return("Medíocre")
    if (val >= 0.50) return("Fraco")
    return("Inadmissível")
  }

  kmo_global_nivel <- .classifica_kmo(kmo_global)
  kmo_ind_nivel    <- vapply(kmo_ind, .classifica_kmo, character(1))

  # 3. Ajuste da ACP via prcomp
  mod_prcomp <- stats::prcomp(df_num, center = TRUE, scale. = padronizar)
  autovalores_obs <- mod_prcomp$sdev^2
  var_expl_pct    <- (autovalores_obs / sum(autovalores_obs)) * 100
  var_acum_pct    <- cumsum(var_expl_pct)

  df_autovalores <- data.frame(
    Componente   = paste0("PC", seq_len(p)),
    Autovalor    = autovalores_obs,
    Var_Expl_Pct = var_expl_pct,
    Var_Acum_Pct = var_acum_pct,
    stringsAsFactors = FALSE
  )

  # 4. Critérios de Retenção
  # a) Kaiser (autovalor > 1 quando padronizado; se não, autovalor > média dos autovalores)
  corte_kaiser <- if (padronizar) 1.0 else mean(autovalores_obs)
  kaiser_sugerido <- max(1, sum(autovalores_obs > corte_kaiser))

  # b) Análise Paralela de Horn
  ap_p95 <- rep(NA, p)
  ap_sugerido <- NA
  if (analise_paralela) {
    sim_ev <- matrix(0, nrow = n_sim, ncol = p)
    for (s in seq_len(n_sim)) {
      mat_sim <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
      cor_sim <- stats::cor(mat_sim)
      ev_sim  <- eigen(cor_sim, only.values = TRUE)$values
      sim_ev[s, ] <- ev_sim
    }
    ap_p95 <- apply(sim_ev, 2, stats::quantile, probs = 0.95)
    ap_sugerido <- max(1, sum(autovalores_obs > ap_p95))
  }

  # c) Variância Acumulada >= 80%
  n_80 <- which(var_acum_pct >= 80)
  n_80_sugerido <- if (length(n_80) > 0) n_80[1] else p

  # Número final de componentes retidos
  n_retidos <- if (!is.null(n_comp)) {
    min(max(1, as.integer(n_comp)), p)
  } else {
    kaiser_sugerido
  }

  # 5. Cargas dos Componentes (Loadings)
  mat_cargas <- mod_prcomp$rotation[, seq_len(n_retidos), drop = FALSE]
  df_cargas  <- data.frame(
    Variavel = var_nomes,
    round(mat_cargas, 4),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  colnames(df_cargas) <- c("Variavel", paste0("PC", seq_len(n_retidos)))

  # Escores
  df_escores <- as.data.frame(mod_prcomp$x[, seq_len(n_retidos), drop = FALSE])
  colnames(df_escores) <- paste0("PC", seq_len(n_retidos))
  rownames(df_escores) <- rownames(df_num)

  # Helpers de formatação
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = ",")
  .fmt_p <- function(x) {
    if (is.na(x)) return("-")
    if (x < 0.001) return("< 0,001")
    formatC(x, format = "f", digits = 3, decimal.mark = ",")
  }
  .ast <- function(p) {
    if (is.na(p) || p >= 0.10) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    return("*")
  }

  # ==========================================================================
  # IMPRESSÃO NO CONSOLE
  # ==========================================================================
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("ANÁLISE DE COMPONENTES PRINCIPAIS (ACP / PCA)")
  } else {
    cat("\n── ANÁLISE DE COMPONENTES PRINCIPAIS (ACP / PCA) ──\n")
  }

  # SEÇÃO: ADEQUAÇÃO DA AMOSTRA (PRESSUPOSTOS)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("ADEQUAÇÃO DA AMOSTRA (PRESSUPOSTOS)")
  } else {
    cat("\n▶ ADEQUAÇÃO DA AMOSTRA (PRESSUPOSTOS)\n")
  }

  w_crit <- c(crit = 42, val = 26)
  sep_crit <- paste0("  ", strrep("─", sum(w_crit) + 3))
  cat(sep_crit, "\n")
  cat("  ", .pad("Critério", w_crit["crit"], "left"), .pad("Valor", w_crit["val"], "center"), "\n", sep = "")
  cat(sep_crit, "\n")

  kmo_glob_txt <- if (!is.na(kmo_global)) sprintf("%s (%s)", .fmt_f(kmo_global, 3), kmo_global_nivel) else "Não aplicável"
  bart_txt <- if (!is.na(bart_chi2)) sprintf("%s (df = %d, p %s)%s", .fmt_f(bart_chi2, 2), bart_df, .fmt_p(bart_p), .ast(bart_p)) else "Não aplicável"

  cat("  ", .pad("KMO Global (Kaiser-Meyer-Olkin)", w_crit["crit"], "left"), .pad(kmo_glob_txt, w_crit["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Teste de Esfericidade de Bartlett (χ²)", w_crit["crit"], "left"), .pad(bart_txt, w_crit["val"], "center"), "\n", sep = "")
  cat(sep_crit, "\n")
  cat("  * Bartlett significativo (p < 0,05) indica correlações suficientes entre as variáveis.\n")
  cat("  * KMO Global: >= 0,90 Excelente | 0,80-0,89 Meritório | 0,70-0,79 Médio | 0,60-0,69 Medíocre | < 0,50 Inadmissível\n\n")

  # Tabela KMO Individual por Variável
  cat("  KMO Individual por Variável:\n")
  w_kmo <- c(var = 28, kmo = 12, niv = 18)
  sep_kmo <- paste0("  ", strrep("─", sum(w_kmo) + 3))
  cat(sep_kmo, "\n")
  cat("  ", .pad("Variável", w_kmo["var"], "left"), .pad("KMO", w_kmo["kmo"], "center"), .pad("Nível do KMO", w_kmo["niv"], "center"), "\n", sep = "")
  cat(sep_kmo, "\n")
  for (i in seq_along(var_nomes)) {
    cat("  ",
        .pad(var_nomes[i], w_kmo["var"], "left"),
        .pad(.fmt_f(kmo_ind[i], 3), w_kmo["kmo"], "center"),
        .pad(kmo_ind_nivel[i], w_kmo["niv"], "center"),
        "\n", sep = "")
  }
  cat(sep_kmo, "\n")
  cat("  Nota: >= 0,90 Excelente | 0,80 - 0,89 Meritório | 0,70 - 0,79 Médio\n")
  cat("        0,60 - 0,69 Medíocre  | 0,50 - 0,59 Fraco     | < 0,50 Inadmissível\n\n")

  # SEÇÃO: VARIÂNCIA EXPLICADA POR COMPONENTE
  if (exists(".print_topico", mode = "function")) {
    .print_topico("VARIÂNCIA EXPLICADA POR COMPONENTE")
  } else {
    cat("\n▶ VARIÂNCIA EXPLICADA POR COMPONENTE\n")
  }

  w_var <- c(comp = 14, auto = 15, expl = 22, acum = 22)
  sep_var <- paste0("  ", strrep("─", sum(w_var) + 3))
  cat(sep_var, "\n")
  cat("  ",
      .pad("Componente", w_var["comp"], "left"),
      .pad("Autovalor", w_var["auto"], "center"),
      .pad("Var. Explicada (%)", w_var["expl"], "center"),
      .pad("Var. Acumulada (%)", w_var["acum"], "center"),
      "\n", sep = "")
  cat(sep_var, "\n")
  for (i in seq_len(p)) {
    cat("  ",
        .pad(paste0("PC", i), w_var["comp"], "left"),
        .pad(.fmt_f(autovalores_obs[i], 2), w_var["auto"], "center"),
        .pad(.fmt_f(var_expl_pct[i], 2), w_var["expl"], "center"),
        .pad(.fmt_f(var_acum_pct[i], 2), w_var["acum"], "center"),
        "\n", sep = "")
  }
  cat(sep_var, "\n\n")

  # SEÇÃO: CRITÉRIOS DE RETENÇÃO DE COMPONENTES
  if (exists(".print_topico", mode = "function")) {
    .print_topico("CRITÉRIOS DE RETENÇÃO DE COMPONENTES")
  } else {
    cat("\n▶ CRITÉRIOS DE RETENÇÃO DE COMPONENTES\n")
  }

  w_ret <- c(crit = 22, sug = 16, just = 36)
  sep_ret <- paste0("  ", strrep("─", sum(w_ret) + 3))
  cat(sep_ret, "\n")
  cat("  ", .pad("Critério", w_ret["crit"], "left"), .pad("Sugerido", w_ret["sug"], "center"), .pad("Justificativa", w_ret["just"], "left"), "\n", sep = "")
  cat(sep_ret, "\n")

  txt_kaiser_just <- sprintf("Autovalores maiores que %s (%s)",
                             .fmt_f(corte_kaiser, 1),
                             paste0("PC", seq_len(kaiser_sugerido), collapse = " e "))
  cat("  ", .pad("Critério de Kaiser", w_ret["crit"], "left"),
      .pad(sprintf("%d componente(s)", kaiser_sugerido), w_ret["sug"], "center"),
      .pad(txt_kaiser_just, w_ret["just"], "left"), "\n", sep = "")

  if (analise_paralela) {
    cat("  ", .pad("Análise Paralela", w_ret["crit"], "left"),
        .pad(sprintf("%d componente(s)", ap_sugerido), w_ret["sug"], "center"),
        .pad("Autovalor empírico > P95 da simulação", w_ret["just"], "left"), "\n", sep = "")
  }

  txt_80_just <- sprintf("Acumulam %s%% da variância total", .fmt_f(var_acum_pct[n_80_sugerido], 2))
  cat("  ", .pad("Variância >= 80%", w_ret["crit"], "left"),
      .pad(sprintf("%d componente(s)", n_80_sugerido), w_ret["sug"], "center"),
      .pad(txt_80_just, w_ret["just"], "left"), "\n", sep = "")
  cat(sep_ret, "\n\n")

  # Tabela 2: Comparação Autovalores Obs vs Análise Paralela (Horn)
  if (analise_paralela) {
    cat("  Comparação: Autovalores Observados vs Análise Paralela (Horn):\n")
    w_ap <- c(comp = 14, obs = 18, sim = 22, acum = 20)
    sep_ap <- paste0("  ", strrep("─", sum(w_ap) + 3))
    cat(sep_ap, "\n")
    cat("  ",
        .pad("Componente", w_ap["comp"], "left"),
        .pad("Autovalor Obs.", w_ap["obs"], "center"),
        .pad("Autovalor Sim. (AP)", w_ap["sim"], "center"),
        .pad("Var. Acumulada (%)", w_ap["acum"], "center"),
        "\n", sep = "")
    cat(sep_ap, "\n")
    n_mostrar_ap <- min(p, max(kaiser_sugerido + 2, 5))
    for (i in seq_len(n_mostrar_ap)) {
      cat("  ",
          .pad(paste0("PC", i), w_ap["comp"], "left"),
          .pad(.fmt_f(autovalores_obs[i], 2), w_ap["obs"], "center"),
          .pad(.fmt_f(ap_p95[i], 2), w_ap["sim"], "center"),
          .pad(.fmt_f(var_acum_pct[i], 2), w_ap["acum"], "center"),
          "\n", sep = "")
    }
    cat(sep_ap, "\n")
    cat("  * Reter componentes onde Autovalor Obs. > Autovalor Simulado (P95).\n\n")
  }

  # SEÇÃO: CARGAS DOS COMPONENTES (LOADINGS)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("CARGAS DOS COMPONENTES (LOADINGS)")
  } else {
    cat("\n▶ CARGAS DOS COMPONENTES (LOADINGS)\n")
  }

  w_var_col <- 24
  w_pc_col  <- max(10, min(14, floor((54 - w_var_col) / n_retidos)))
  sep_load  <- paste0("  ", strrep("─", w_var_col + n_retidos * (w_pc_col + 1)))

  cat(sep_load, "\n")
  cat("  ", .pad("Variável", w_var_col, "left"))
  for (j in seq_len(n_retidos)) {
    cat(.pad(paste0("PC", j), w_pc_col, "center"))
  }
  cat("\n")
  cat(sep_load, "\n")

  for (i in seq_along(var_nomes)) {
    cat("  ", .pad(var_nomes[i], w_var_col, "left"))
    for (j in seq_len(n_retidos)) {
      val <- mat_cargas[i, j]
      sinal_txt <- if (val >= 0) sprintf("+%s", .fmt_f(val, 3)) else sprintf("-%s", .fmt_f(abs(val), 3))
      cat(.pad(sinal_txt, w_pc_col, "center"))
    }
    cat("\n")
  }
  cat(sep_load, "\n")
  cat(sprintf("  * Exibindo cargas dos %d componente(s) retido(s).\n\n", n_retidos))

  # SEÇÃO: INTERPRETAÇÃO (estilo texto fluido, sem bolinhas, parágrafos espaçados)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO")
  } else {
    cat("\n▶ INTERPRETAÇÃO\n")
  }

  padrao_txt <- if (padronizar) "padronizados (matriz de correlação)" else "não padronizados (matriz de covariância)"
  bart_interp <- if (!is.na(bart_p) && bart_p < 0.05) {
    sprintf("O teste de esfericidade de Bartlett foi estatisticamente significativo (p %s),\n  rejeitando a hipótese de que a matriz de correlação é uma identidade.", .fmt_p(bart_p))
  } else {
    "O teste de esfericidade de Bartlett não foi significativo, indicando baixas correlações entre as variáveis."
  }
  kmo_interp <- if (!is.na(kmo_global)) {
    sprintf("O KMO global de %s (%s) indicou adequação satisfatória da amostra para redução dimensional.",
            .fmt_f(kmo_global, 3), kmo_global_nivel)
  } else ""

  cat(sprintf("  A ACP foi ajustada com os dados %s.\n  %s\n  %s\n\n",
              padrao_txt, bart_interp, kmo_interp))

  cat(sprintf("  Pelo Critério de Kaiser (autovalor > %s), foram retidos %d componente(s),\n  que juntos explicam %s%% da variabilidade total contida nas %d variáveis originais.\n\n",
              .fmt_f(corte_kaiser, 1), n_retidos, .fmt_f(var_acum_pct[n_retidos], 2), p))

  # Destaque de variáveis para PC1
  cargas_pc1 <- mat_cargas[, 1]
  vars_pos1  <- names(sort(cargas_pc1[cargas_pc1 > 0.25], decreasing = TRUE))
  vars_neg1  <- names(sort(cargas_pc1[cargas_pc1 < -0.25]))

  pc1_txt <- sprintf("  O primeiro componente (PC1) é responsável por %s%% da variância explicada.", .fmt_f(var_expl_pct[1], 2))
  if (length(vars_pos1) > 0 || length(vars_neg1) > 0) {
    detalhes <- c()
    if (length(vars_pos1) > 0) detalhes <- c(detalhes, sprintf("cargas positivas elevadas em %s", paste(vars_pos1, collapse = ", ")))
    if (length(vars_neg1) > 0) detalhes <- c(detalhes, sprintf("cargas negativas em %s", paste(vars_neg1, collapse = ", ")))
    pc1_txt <- paste0(pc1_txt, sprintf("\n  Ele é caracterizado principalmente por %s.", paste(detalhes, collapse = " contrapostas a ")))
  }
  cat(pc1_txt, "\n\n")

  # Destaque para PC2 se houver
  if (n_retidos >= 2) {
    cargas_pc2 <- mat_cargas[, 2]
    vars_pos2  <- names(sort(cargas_pc2[cargas_pc2 > 0.25], decreasing = TRUE))
    vars_neg2  <- names(sort(cargas_pc2[cargas_pc2 < -0.25]))

    pc2_txt <- sprintf("  O segundo componente (PC2) responde por %s%% da variância.", .fmt_f(var_expl_pct[2], 2))
    if (length(vars_pos2) > 0 || length(vars_neg2) > 0) {
      detalhes2 <- c()
      if (length(vars_pos2) > 0) detalhes2 <- c(detalhes2, sprintf("cargas positivas em %s", paste(vars_pos2, collapse = ", ")))
      if (length(vars_neg2) > 0) detalhes2 <- c(detalhes2, sprintf("cargas negativas em %s", paste(vars_neg2, collapse = ", ")))
      pc2_txt <- paste0(pc2_txt, sprintf("\n  Destacam-se %s.", paste(detalhes2, collapse = " em oposição a ")))
    }
    cat(pc2_txt, "\n\n")
  }

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(strrep("─", sum(w_crit) + 3), "\n\n")
  }

  # Montagem do objeto de retorno
  res <- list(
    autovalores = df_autovalores,
    cargas      = df_cargas,
    escores     = df_escores,
    kmo         = list(
      global = kmo_global,
      classificacao_global = kmo_global_nivel,
      por_variavel = data.frame(
        Variavel = var_nomes,
        KMO      = kmo_ind,
        Nivel    = kmo_ind_nivel,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    ),
    bartlett    = list(
      qui_quadrado  = bart_chi2,
      gl            = bart_df,
      p_valor       = bart_p,
      significativo = (!is.na(bart_p) && bart_p < 0.05)
    ),
    criterios   = list(
      kaiser_sugerido  = kaiser_sugerido,
      ap_sugerido      = ap_sugerido,
      n_80_sugerido    = n_80_sugerido,
      n_retidos        = n_retidos,
      autovalores_simulados_p95 = ap_p95
    ),
    prcomp      = mod_prcomp,
    dados_originais = df_num
  )
  class(res) <- c("acp_estatR", "list")

  if (grafico) {
    tryCatch({
      p_g <- biplot_acp(res)
      print(p_g)
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o biplot: ", e$message)
    })
  }

  invisible(res)
}

#' @title Biplot de ACP
#' @description Gera gráfico Biplot 2D (escores das observações e vetores das cargas das variáveis).
#' @param modelo Objeto retornado pela função \code{acp()}.
#' @param grupo Vetor ou fator com grupos das observações (opcional).
#' @param elipses Lógico. Se TRUE e 'grupo' informado, traça elipses de confiança 95%.
#' @param escala Numérico. Fator de escala das setas das variáveis em relação aos escores.
#' @param rotulos Lógico. Se TRUE, exibe os nomes das variáveis nas setas.
#' @param pc_x Inteiro. Componente do eixo X (padrão: 1).
#' @param pc_y Inteiro. Componente do eixo Y (padrão: 2).
#' @export
biplot_acp <- function(modelo,
                       grupo   = NULL,
                       elipses = TRUE,
                       escala  = 0.7,
                       rotulos = TRUE,
                       pc_x    = 1,
                       pc_y    = 2) {

  if (!inherits(modelo, "acp_estatR")) {
    stop("O argumento 'modelo' deve ser um objeto gerado pela função acp().")
  }

  var_expl <- modelo$autovalores$Var_Expl_Pct
  nome_x <- paste0("PC", pc_x)
  nome_y <- paste0("PC", pc_y)

  rot <- modelo$prcomp$rotation
  if (pc_x > ncol(rot) || pc_y > ncol(rot)) {
    stop("Os componentes selecionados excedem o número de componentes disponíveis.")
  }

  scores_df <- as.data.frame(modelo$prcomp$x[, c(pc_x, pc_y), drop = FALSE])
  colnames(scores_df) <- c("EixoX", "EixoY")

  cargas_df <- as.data.frame(rot[, c(pc_x, pc_y), drop = FALSE])
  colnames(cargas_df) <- c("EixoX", "EixoY")
  cargas_df$Variavel  <- rownames(cargas_df)

  # Fator multiplicador de escala para as setas ficarem proporcionais
  max_scores <- max(abs(scores_df$EixoX), abs(scores_df$EixoY), na.rm = TRUE)
  max_cargas <- max(abs(cargas_df$EixoX), abs(cargas_df$EixoY), na.rm = TRUE)
  mult_escala <- (max_scores / max_cargas) * escala

  cargas_df$X_fim <- cargas_df$EixoX * mult_escala
  cargas_df$Y_fim <- cargas_df$EixoY * mult_escala

  lab_x <- sprintf("PC%d (%s%%)", pc_x, formatC(var_expl[pc_x], format = "f", digits = 1, decimal.mark = ","))
  lab_y <- sprintf("PC%d (%s%%)", pc_y, formatC(var_expl[pc_y], format = "f", digits = 1, decimal.mark = ","))

  p <- ggplot2::ggplot()

  if (!is.null(grupo)) {
    scores_df$Grupo <- as.factor(grupo)
    if (elipses && nlevels(scores_df$Grupo) > 1) {
      p <- p + ggplot2::stat_ellipse(
        data = scores_df,
        ggplot2::aes(x = EixoX, y = EixoY, color = Grupo, fill = Grupo),
        geom = "polygon", alpha = 0.1, linewidth = 0.5, level = 0.95
      )
    }
    p <- p + ggplot2::geom_point(data = scores_df, ggplot2::aes(x = EixoX, y = EixoY, color = Grupo), size = 2.2, alpha = 0.75)
  } else {
    p <- p + ggplot2::geom_point(data = scores_df, ggplot2::aes(x = EixoX, y = EixoY), color = "#555555", size = 2, alpha = 0.65)
  }

  # Linhas dos eixos centrais
  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4)

  # Setas das variáveis
  p <- p + ggplot2::geom_segment(
    data = cargas_df,
    ggplot2::aes(x = 0, y = 0, xend = X_fim, yend = Y_fim),
    color = "#D90429", linewidth = 0.8,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"), type = "closed")
  )

  # Rótulos
  if (rotulos) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = cargas_df,
        ggplot2::aes(x = X_fim, y = Y_fim, label = Variavel),
        color = "#1D3557", size = 4, fontface = "bold", box.padding = 0.3
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = cargas_df,
        ggplot2::aes(x = X_fim * 1.08, y = Y_fim * 1.08, label = Variavel),
        color = "#1D3557", size = 4, fontface = "bold"
      )
    }
  }

  p <- p + ggplot2::labs(
    title    = "Biplot — ACP (Análise de Componentes Principais)",
    subtitle = paste0(lab_x, " vs ", lab_y),
    x = lab_x, y = lab_y
  )

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2)
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}

#' @title Scree Plot de ACP
#' @description Gera gráfico Scree Plot com linha de Kaiser (autovalor = 1) e Análise Paralela de Horn.
#' @param modelo Objeto retornado pela função \code{acp()}.
#' @param analise_paralela Lógico. Se TRUE, exibe a linha dos autovalores simulados (P95).
#' @export
screeplot_acp <- function(modelo, analise_paralela = TRUE) {

  if (!inherits(modelo, "acp_estatR")) {
    stop("O argumento 'modelo' deve ser um objeto gerado pela função acp().")
  }

  df_scree <- data.frame(
    Componente_num = seq_len(nrow(modelo$autovalores)),
    Componente     = modelo$autovalores$Componente,
    Autovalor      = modelo$autovalores$Autovalor
  )

  tem_ap <- analise_paralela && !is.null(modelo$criterios$autovalores_simulados_p95) &&
            !any(is.na(modelo$criterios$autovalores_simulados_p95))

  cores_legenda <- c(
    "Autovalores Observados"             = "#1D3557",
    "Critério de Kaiser (Autovalor = 1)" = "#D90429"
  )
  linetypes_legenda <- c(
    "Autovalores Observados"             = "solid",
    "Critério de Kaiser (Autovalor = 1)" = "dashed"
  )

  if (tem_ap) {
    cores_legenda["Análise Paralela de Horn (P95)"]     <- "#457B9D"
    linetypes_legenda["Análise Paralela de Horn (P95)"] <- "dotted"
  }

  p <- ggplot2::ggplot(df_scree, ggplot2::aes(x = Componente_num, y = Autovalor)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1.0, color = "Critério de Kaiser (Autovalor = 1)",
                                     linetype = "Critério de Kaiser (Autovalor = 1)"),
                        linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(color = "Autovalores Observados", linetype = "Autovalores Observados"),
                       linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(color = "Autovalores Observados"), size = 3)

  if (tem_ap) {
    df_ap <- data.frame(
      Componente_num = seq_along(modelo$criterios$autovalores_simulados_p95),
      Autovalor_AP   = modelo$criterios$autovalores_simulados_p95
    )
    p <- p +
      ggplot2::geom_line(data = df_ap,
                         ggplot2::aes(x = Componente_num, y = Autovalor_AP,
                                      color = "Análise Paralela de Horn (P95)",
                                      linetype = "Análise Paralela de Horn (P95)"),
                         linewidth = 1) +
      ggplot2::geom_point(data = df_ap,
                          ggplot2::aes(x = Componente_num, y = Autovalor_AP,
                                       color = "Análise Paralela de Horn (P95)"),
                          size = 2.5, shape = 17)
  }

  p <- p +
    ggplot2::scale_x_continuous(breaks = df_scree$Componente_num, labels = df_scree$Componente) +
    ggplot2::scale_color_manual(name = "Critério", values = cores_legenda) +
    ggplot2::scale_linetype_manual(name = "Critério", values = linetypes_legenda) +
    ggplot2::labs(
      title = "Scree Plot — ACP",
      x     = "Componente Principal",
      y     = "Autovalor"
    ) +
    meu_tema(grade = "dupla") +
    ggplot2::theme(legend.position = c(0.85, 0.85), legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.7), color = "gray80"))

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}


# ==============================================================================
# ANÁLISE FATORIAL EXPLORATÓRIA (AFE)
# ==============================================================================

#' @title Análise Fatorial Exploratória (AFE)
#' @description Realiza Análise Fatorial Exploratória (AFE) completa com diagnóstico
#' de adequação amostral (KMO global e individual, Teste de Bartlett), critérios de retenção
#' de fatores (Kaiser e Análise Paralela de Horn), cargas fatoriais rotacionadas,
#' comunalidades (h²), unicidades (u²), escores fatoriais e interpretação automatizada.
#' @param dados Data frame ou matriz contendo variáveis numéricas.
#' @param n_fatores Inteiro opcional. Número de fatores a reter. Se NULL (padrão),
#'   determina automaticamente pelo Critério de Kaiser (autovalor > 1).
#' @param rotacao Caractere indicando o método de rotação: \code{"varimax"} (padrão),
#'   \code{"oblimin"}, \code{"promax"} ou \code{"none"}.
#' @param metodo Caractere indicando o método de extração fatorial: \code{"pa"} (Eixos Principais - padrão),
#'   \code{"minres"} (Mínimos Resíduos) ou \code{"ml"} (Máxima Verossimilhança).
#' @param analise_paralela Lógico. Se TRUE (padrão), calcula a Análise Paralela de Horn.
#' @param n_sim Inteiro. Número de simulações para a Análise Paralela (padrão: 100).
#' @param corte_carga Numérico. Valor de corte para considerar a carga fatorial saliente na interpretação (padrão: 0.30).
#' @param grafico Lógico. Se TRUE, exibe automaticamente o gráfico de cargas fatoriais.
#' @return Objeto da classe \code{afe_estatR} contendo:
#'   \item{cargas}{Data frame com variáveis, cargas fatoriais por fator, comunalidades (h²) e unicidades (u²)}
#'   \item{comunalidades}{Vetor com as comunalidades h² de cada variável}
#'   \item{unicidades}{Vetor com as unicidades u² de cada variável}
#'   \item{variancia}{Data frame com autovalores fatoriais, variância explicada e acumulada}
#'   \item{escores}{Data frame com os escores fatoriais das observações}
#'   \item{kmo}{Lista com KMO global, nível e data frame individual por variável}
#'   \item{bartlett}{Lista com estatística qui-quadrado, graus de liberdade e p-valor}
#'   \item{criterios}{Lista com resumo dos critérios e autovalores simulados de Horn}
#'   \item{rotacao}{Método de rotação utilizado}
#'   \item{metodo}{Método de extração utilizado}
#'   \item{fa}{Objeto retornado pela função psych::fa}
#' @export
analise_fatorial <- function(dados,
                             n_fatores        = NULL,
                             rotacao          = "varimax",
                             metodo           = "pa",
                             analise_paralela = TRUE,
                             n_sim            = 100,
                             corte_carga      = 0.30,
                             grafico          = TRUE) {

  if (!requireNamespace("psych", quietly = TRUE)) {
    stop("O pacote 'psych' e necessario para a Analise Fatorial Exploratoria. Instale com install.packages('psych').")
  }

  # 1. Validação de dados
  if (!is.data.frame(dados) && !is.matrix(dados)) {
    stop("O argumento 'dados' deve ser um data frame ou matriz.")
  }

  dados_df <- as.data.frame(dados)
  cols_num <- vapply(dados_df, is.numeric, logical(1))
  if (sum(cols_num) < 2) {
    stop("São necessárias pelo menos 2 variáveis numéricas para a Análise Fatorial.")
  }

  df_num <- dados_df[, cols_num, drop = FALSE]
  casos_completos <- stats::complete.cases(df_num)
  if (sum(casos_completos) < nrow(df_num)) {
    df_num <- df_num[casos_completos, , drop = FALSE]
    dados_df <- dados_df[casos_completos, , drop = FALSE]
  }

  n <- nrow(df_num)
  p <- ncol(df_num)
  var_nomes <- colnames(df_num)

  # 2. Diagnóstico de Adequação: KMO e Teste de Bartlett
  mat_cor <- stats::cor(df_num)

  # Bartlett
  det_r <- det(mat_cor)
  if (is.na(det_r) || det_r <= 0) {
    bart_chi2 <- NA; bart_df <- p * (p - 1) / 2; bart_p <- NA
  } else {
    bart_df   <- p * (p - 1) / 2
    bart_chi2 <- -( (n - 1) - (2 * p + 5) / 6 ) * log(det_r)
    bart_p    <- stats::pchisq(bart_chi2, df = bart_df, lower.tail = FALSE)
  }

  # KMO
  inv_r <- tryCatch(solve(mat_cor), error = function(e) NULL)
  if (is.null(inv_r)) {
    kmo_global <- NA
    kmo_ind <- rep(NA, p)
  } else {
    anti_img <- -inv_r / sqrt(outer(diag(inv_r), diag(inv_r)))
    diag(anti_img) <- 1

    r2 <- mat_cor^2
    a2 <- anti_img^2
    diag(r2) <- 0
    diag(a2) <- 0

    kmo_global <- sum(r2) / (sum(r2) + sum(a2))
    kmo_ind <- rowSums(r2) / (rowSums(r2) + rowSums(a2))
  }

  .classifica_kmo <- function(val) {
    if (is.na(val)) return("Indeterminado")
    if (val >= 0.90) return("Excelente")
    if (val >= 0.80) return("Meritório")
    if (val >= 0.70) return("Médio")
    if (val >= 0.60) return("Medíocre")
    if (val >= 0.50) return("Fraco")
    return("Inadmissível")
  }

  kmo_global_nivel <- .classifica_kmo(kmo_global)
  kmo_ind_nivel    <- vapply(kmo_ind, .classifica_kmo, character(1))

  # 3. Autovalores da matriz de correlação e Critérios de Retenção
  ev_cor <- eigen(mat_cor, only.values = TRUE)$values
  var_expl_pca <- (ev_cor / sum(ev_cor)) * 100
  var_acum_pca <- cumsum(var_expl_pca)

  kaiser_sugerido <- max(1, sum(ev_cor > 1.0))

  # Análise Paralela de Horn
  ap_p95 <- rep(NA, p)
  ap_sugerido <- NA
  if (analise_paralela) {
    sim_ev <- matrix(0, nrow = n_sim, ncol = p)
    for (s in seq_len(n_sim)) {
      mat_sim <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
      cor_sim <- stats::cor(mat_sim)
      ev_sim  <- eigen(cor_sim, only.values = TRUE)$values
      sim_ev[s, ] <- ev_sim
    }
    ap_p95 <- apply(sim_ev, 2, stats::quantile, probs = 0.95)
    ap_sugerido <- max(1, sum(ev_cor > ap_p95))
  }

  # Variância acumulada >= 70%
  n_70 <- which(var_acum_pca >= 70)
  n_70_sugerido <- if (length(n_70) > 0) n_70[1] else p

  # Número final de fatores
  n_retidos <- if (!is.null(n_fatores)) {
    min(max(1, as.integer(n_fatores)), p)
  } else {
    kaiser_sugerido
  }

  # 4. Ajuste da Análise Fatorial via psych::fa
  rot_ajuste <- if (n_retidos == 1) "none" else rotacao
  fa_mod <- suppressWarnings(tryCatch(
    psych::fa(df_num, nfactors = n_retidos, rotate = rot_ajuste, fm = metodo),
    error = function(e) stop("Erro ao ajustar o modelo de Análise Fatorial: ", e$message)
  ))

  # Extrair Cargas, Comunalidades e Unicidades
  mat_cargas_raw <- unclass(fa_mod$loadings)
  comunalidades  <- fa_mod$communality
  unicidades     <- fa_mod$uniquenesses

  fatores_nomes <- paste("Fator", seq_len(n_retidos))
  colnames(mat_cargas_raw) <- fatores_nomes

  df_cargas <- data.frame(
    Variavel = var_nomes,
    round(mat_cargas_raw, 4),
    Comunalidade_h2 = round(comunalidades, 4),
    Unicidade_u2    = round(unicidades, 4),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  colnames(df_cargas) <- c("Variavel", fatores_nomes, "Comunalidade (h²)", "Unicidade (u²)")

  # Variância explicada pelos fatores rotacionados
  # fa_mod$Vaccounted tem SS loadings, Proportion Var, Cumulative Var
  vacct <- fa_mod$Vaccounted
  if (!is.null(vacct) && is.matrix(vacct) && nrow(vacct) >= 3) {
    ss_loadings  <- vacct[1, seq_len(n_retidos)]
    prop_var_pct <- vacct[2, seq_len(n_retidos)] * 100
    cum_var_pct  <- vacct[3, seq_len(n_retidos)] * 100
  } else {
    # Fallback se Vaccounted não estiver no formato esperado
    ss_loadings  <- colSums(mat_cargas_raw^2)
    prop_var_pct <- (ss_loadings / p) * 100
    cum_var_pct  <- cumsum(prop_var_pct)
  }

  df_variancia <- data.frame(
    Fator               = fatores_nomes,
    Autovalor_Fatorial  = ss_loadings,
    Var_Expl_Pct        = prop_var_pct,
    Var_Acum_Pct        = cum_var_pct,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Escores fatoriais
  df_escores <- as.data.frame(fa_mod$scores)
  if (ncol(df_escores) == n_retidos) {
    colnames(df_escores) <- fatores_nomes
  }
  rownames(df_escores) <- rownames(df_num)

  # Helpers de formatação
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = ",")
  .fmt_p <- function(x) {
    if (is.na(x)) return("-")
    if (x < 0.001) return("< 0,001")
    formatC(x, format = "f", digits = 3, decimal.mark = ",")
  }
  .ast <- function(p) {
    if (is.na(p) || p >= 0.10) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    return("*")
  }

  # ==========================================================================
  # IMPRESSÃO NO CONSOLE
  # ==========================================================================
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("ANÁLISE FATORIAL EXPLORATÓRIA (AFE)")
  } else {
    cat("\n── ANÁLISE FATORIAL EXPLORATÓRIA (AFE) ──\n")
  }

  # SEÇÃO: ADEQUAÇÃO DA AMOSTRA (PRESSUPOSTOS)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("ADEQUAÇÃO DA AMOSTRA (PRESSUPOSTOS)")
  } else {
    cat("\n▶ ADEQUAÇÃO DA AMOSTRA (PRESSUPOSTOS)\n")
  }

  w_crit <- c(crit = 42, val = 26)
  sep_crit <- paste0("  ", strrep("─", sum(w_crit) + 3))
  cat(sep_crit, "\n")
  cat("  ", .pad("Critério", w_crit["crit"], "left"), .pad("Valor", w_crit["val"], "center"), "\n", sep = "")
  cat(sep_crit, "\n")

  kmo_glob_txt <- if (!is.na(kmo_global)) sprintf("%s (%s)", .fmt_f(kmo_global, 3), kmo_global_nivel) else "Não aplicável"
  bart_txt <- if (!is.na(bart_chi2)) sprintf("%s (df = %d, p %s)%s", .fmt_f(bart_chi2, 2), bart_df, .fmt_p(bart_p), .ast(bart_p)) else "Não aplicável"

  cat("  ", .pad("KMO Global (Kaiser-Meyer-Olkin)", w_crit["crit"], "left"), .pad(kmo_glob_txt, w_crit["val"], "center"), "\n", sep = "")
  cat("  ", .pad("Teste de Esfericidade de Bartlett (χ²)", w_crit["crit"], "left"), .pad(bart_txt, w_crit["val"], "center"), "\n", sep = "")
  cat(sep_crit, "\n")
  cat("  * Bartlett significativo (p < 0,05) indica correlações suficientes para fatoração.\n")
  cat("  * KMO Global: >= 0,90 Excelente | 0,80-0,89 Meritório | 0,70-0,79 Médio | 0,60-0,69 Medíocre | < 0,50 Inadmissível\n\n")

  # Tabela KMO Individual por Variável
  cat("  KMO Individual por Variável:\n")
  w_kmo <- c(var = 28, kmo = 12, niv = 18)
  sep_kmo <- paste0("  ", strrep("─", sum(w_kmo) + 3))
  cat(sep_kmo, "\n")
  cat("  ", .pad("Variável", w_kmo["var"], "left"), .pad("KMO", w_kmo["kmo"], "center"), .pad("Nível do KMO", w_kmo["niv"], "center"), "\n", sep = "")
  cat(sep_kmo, "\n")
  for (i in seq_along(var_nomes)) {
    cat("  ",
        .pad(var_nomes[i], w_kmo["var"], "left"),
        .pad(.fmt_f(kmo_ind[i], 3), w_kmo["kmo"], "center"),
        .pad(kmo_ind_nivel[i], w_kmo["niv"], "center"),
        "\n", sep = "")
  }
  cat(sep_kmo, "\n")
  cat("  Nota: >= 0,90 Excelente | 0,80 - 0,89 Meritório | 0,70 - 0,79 Médio\n")
  cat("        0,60 - 0,69 Medíocre  | 0,50 - 0,59 Fraco     | < 0,50 Inadmissível\n\n")

  # SEÇÃO: CRITÉRIOS DE RETENÇÃO DE FATORES
  if (exists(".print_topico", mode = "function")) {
    .print_topico("CRITÉRIOS DE RETENÇÃO DE FATORES")
  } else {
    cat("\n▶ CRITÉRIOS DE RETENÇÃO DE FATORES\n")
  }

  w_ret <- c(crit = 22, sug = 16, just = 38)
  sep_ret <- paste0("  ", strrep("─", sum(w_ret) + 3))
  cat(sep_ret, "\n")
  cat("  ", .pad("Critério", w_ret["crit"], "left"), .pad("Sugerido", w_ret["sug"], "center"), .pad("Justificativa", w_ret["just"], "left"), "\n", sep = "")
  cat(sep_ret, "\n")

  txt_kaiser_just <- sprintf("Autovalores maiores que 1,0 (%s)",
                             paste0("F", seq_len(kaiser_sugerido), collapse = " e "))
  cat("  ", .pad("Critério de Kaiser", w_ret["crit"], "left"),
      .pad(sprintf("%d fator(es)", kaiser_sugerido), w_ret["sug"], "center"),
      .pad(txt_kaiser_just, w_ret["just"], "left"), "\n", sep = "")

  if (analise_paralela) {
    cat("  ", .pad("Análise Paralela", w_ret["crit"], "left"),
        .pad(sprintf("%d fator(es)", ap_sugerido), w_ret["sug"], "center"),
        .pad("Autovalor empírico > P95 da simulação", w_ret["just"], "left"), "\n", sep = "")
  }

  txt_70_just <- sprintf("Acumulam %s%% da variância", .fmt_f(var_acum_pca[n_70_sugerido], 2))
  cat("  ", .pad("Variância >= 70%", w_ret["crit"], "left"),
      .pad(sprintf("%d fator(es)", n_70_sugerido), w_ret["sug"], "center"),
      .pad(txt_70_just, w_ret["just"], "left"), "\n", sep = "")
  cat(sep_ret, "\n\n")

  # Tabela 2: Comparação Autovalores Obs vs Análise Paralela (Horn)
  if (analise_paralela) {
    cat("  Comparação: Autovalores Observados vs Análise Paralela (Horn):\n")
    w_ap <- c(comp = 14, obs = 18, sim = 22, acum = 20)
    sep_ap <- paste0("  ", strrep("─", sum(w_ap) + 3))
    cat(sep_ap, "\n")
    cat("  ",
        .pad("Fator", w_ap["comp"], "left"),
        .pad("Autovalor Obs.", w_ap["obs"], "center"),
        .pad("Autovalor Sim. (AP)", w_ap["sim"], "center"),
        .pad("Var. Acumulada (%)", w_ap["acum"], "center"),
        "\n", sep = "")
    cat(sep_ap, "\n")
    n_mostrar_ap <- min(p, max(kaiser_sugerido + 2, 5))
    for (i in seq_len(n_mostrar_ap)) {
      cat("  ",
          .pad(paste("Fator", i), w_ap["comp"], "left"),
          .pad(.fmt_f(ev_cor[i], 2), w_ap["obs"], "center"),
          .pad(.fmt_f(ap_p95[i], 2), w_ap["sim"], "center"),
          .pad(paste0(.fmt_f(var_acum_pca[i], 2), "%"), w_ap["acum"], "center"),
          "\n", sep = "")
    }
    cat(sep_ap, "\n")
    cat("  * Reter fatores onde Autovalor Obs. > Autovalor Simulado (P95).\n\n")
  }

  # SEÇÃO: VARIÂNCIA EXPLICADA PELOS FATORES ROTACIONADOS
  tit_rot <- if (rot_ajuste == "none") "SEM ROTAÇÃO" else toupper(rot_ajuste)
  if (exists(".print_topico", mode = "function")) {
    .print_topico(sprintf("VARIÂNCIA EXPLICADA PELOS FATORES ROTACIONADOS (%s)", tit_rot))
  } else {
    cat(sprintf("\n▶ VARIÂNCIA EXPLICADA PELOS FATORES ROTACIONADOS (%s)\n", tit_rot))
  }

  w_vfat <- c(fat = 14, auto = 22, expl = 22, acum = 22)
  sep_vfat <- paste0("  ", strrep("─", sum(w_vfat) + 3))
  cat(sep_vfat, "\n")
  cat("  ",
      .pad("Fator", w_vfat["fat"], "left"),
      .pad("Autovalor Fatorial", w_vfat["auto"], "center"),
      .pad("Var. Explicada (%)", w_vfat["expl"], "center"),
      .pad("Var. Acumulada (%)", w_vfat["acum"], "center"),
      "\n", sep = "")
  cat(sep_vfat, "\n")
  for (i in seq_len(n_retidos)) {
    cat("  ",
        .pad(fatores_nomes[i], w_vfat["fat"], "left"),
        .pad(.fmt_f(ss_loadings[i], 2), w_vfat["auto"], "center"),
        .pad(paste0(.fmt_f(prop_var_pct[i], 2), "%"), w_vfat["expl"], "center"),
        .pad(paste0(.fmt_f(cum_var_pct[i], 2), "%"), w_vfat["acum"], "center"),
        "\n", sep = "")
  }
  cat(sep_vfat, "\n\n")

  # SEÇÃO: CARGAS FATORIAIS E COMUNALIDADES
  if (exists(".print_topico", mode = "function")) {
    .print_topico("CARGAS FATORIAIS E COMUNALIDADES (LOADINGS)")
  } else {
    cat("\n▶ CARGAS FATORIAIS E COMUNALIDADES (LOADINGS)\n")
  }

  w_var_col <- 20
  w_fat_col <- max(10, min(13, floor((46 - w_var_col) / n_retidos)))
  w_com_col <- 18
  w_uni_col <- 16

  largura_total <- w_var_col + n_retidos * (w_fat_col + 1) + w_com_col + w_uni_col + 3
  sep_cargas <- paste0("  ", strrep("─", largura_total))

  cat(sep_cargas, "\n")
  cat("  ", .pad("Variável", w_var_col, "left"))
  for (j in seq_len(n_retidos)) {
    cat(.pad(fatores_nomes[j], w_fat_col, "center"))
  }
  cat(.pad("Comunalidade (h²)", w_com_col, "center"))
  cat(.pad("Unicidade (u²)", w_uni_col, "center"))
  cat("\n")
  cat(sep_cargas, "\n")

  for (i in seq_along(var_nomes)) {
    cat("  ", .pad(var_nomes[i], w_var_col, "left"))
    for (j in seq_len(n_retidos)) {
      val <- mat_cargas_raw[i, j]
      sinal_txt <- if (val >= 0) sprintf("+%s", .fmt_f(val, 3)) else sprintf("-%s", .fmt_f(abs(val), 3))
      cat(.pad(sinal_txt, w_fat_col, "center"))
    }
    cat(.pad(.fmt_f(comunalidades[i], 3), w_com_col, "center"))
    cat(.pad(.fmt_f(unicidades[i], 3), w_uni_col, "center"))
    cat("\n")
  }
  cat(sep_cargas, "\n")
  cat(sprintf("  * Cargas rotacionadas pelo método %s (cargas absolutas > %s consideradas salientes).\n",
              tools::toTitleCase(rot_ajuste), .fmt_f(corte_carga, 2)))
  cat("  * Comunalidade (h²): proporção da variância da variável explicada pelos fatores comuns.\n")
  cat("  * Unicidade (u²): proporção da variância não explicada e erro de medida (1 - h²).\n\n")

  # SEÇÃO: INTERPRETAÇÃO
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO")
  } else {
    cat("\n▶ INTERPRETAÇÃO\n")
  }

  metodo_txt <- switch(metodo,
    pa     = "Eixos Principais (Principal Axis Factoring)",
    minres = "Mínimos Resíduos (Minimum Residual)",
    ml     = "Máxima Verossimilhança (Maximum Likelihood)",
    metodo
  )

  bart_interp <- if (!is.na(bart_p) && bart_p < 0.05) {
    sprintf("O teste de esfericidade de Bartlett foi estatisticamente significativo (p %s),\n  rejeitando a hipótese de que a matriz de correlação é uma identidade.", .fmt_p(bart_p))
  } else {
    "O teste de esfericidade de Bartlett não foi significativo, indicando baixas correlações entre as variáveis."
  }
  kmo_interp <- if (!is.na(kmo_global)) {
    sprintf("O KMO global de %s (%s) indicou adequação satisfatória da amostra para fatoração.",
            .fmt_f(kmo_global, 3), kmo_global_nivel)
  } else ""

  cat(sprintf("  A Análise Fatorial Exploratória foi ajustada pelo método dos %s,\n  utilizando rotação %s.\n  %s\n  %s\n\n",
              metodo_txt, tools::toTitleCase(rot_ajuste), bart_interp, kmo_interp))

  cat(sprintf("  Pelo Critério de Kaiser (autovalor > 1,0), foram retidos %d fator(es),\n  que juntos respondem por %s%% da variabilidade comum entre as %d variáveis analisadas.\n\n",
              n_retidos, .fmt_f(cum_var_pct[n_retidos], 2), p))

  # Detalhes Fator 1
  cargas_f1 <- mat_cargas_raw[, 1]
  vars_pos1 <- names(sort(cargas_f1[cargas_f1 >= corte_carga], decreasing = TRUE))
  vars_neg1 <- names(sort(cargas_f1[cargas_f1 <= -corte_carga]))

  f1_txt <- sprintf("  O Fator 1 é responsável por %s%% da variância explicada.", .fmt_f(prop_var_pct[1], 2))
  if (length(vars_pos1) > 0 || length(vars_neg1) > 0) {
    detalhes1 <- c()
    if (length(vars_pos1) > 0) detalhes1 <- c(detalhes1, sprintf("cargas positivas expressivas em %s", paste(vars_pos1, collapse = ", ")))
    if (length(vars_neg1) > 0) detalhes1 <- c(detalhes1, sprintf("cargas negativas em %s", paste(vars_neg1, collapse = ", ")))
    f1_txt <- paste0(f1_txt, sprintf("\n  Ele é caracterizado principalmente por %s.", paste(detalhes1, collapse = " em contraposição a ")))
  }
  cat(f1_txt, "\n\n")

  # Detalhes Fator 2
  if (n_retidos >= 2) {
    cargas_f2 <- mat_cargas_raw[, 2]
    vars_pos2 <- names(sort(cargas_f2[cargas_f2 >= corte_carga], decreasing = TRUE))
    vars_neg2 <- names(sort(cargas_f2[cargas_f2 <= -corte_carga]))

    f2_txt <- sprintf("  O Fator 2 responde por %s%% da variância explicada.", .fmt_f(prop_var_pct[2], 2))
    if (length(vars_pos2) > 0 || length(vars_neg2) > 0) {
      detalhes2 <- c()
      if (length(vars_pos2) > 0) detalhes2 <- c(detalhes2, sprintf("cargas positivas em %s", paste(vars_pos2, collapse = ", ")))
      if (length(vars_neg2) > 0) detalhes2 <- c(detalhes2, sprintf("cargas negativas em %s", paste(vars_neg2, collapse = ", ")))
      f2_txt <- paste0(f2_txt, sprintf("\n  Destacam-se %s.", paste(detalhes2, collapse = " em oposição a ")))
    }
    cat(f2_txt, "\n\n")
  }

  # Avaliação das Comunalidades
  media_h2 <- mean(comunalidades, na.rm = TRUE)
  baixas_h2 <- names(comunalidades[comunalidades < 0.50])
  if (media_h2 >= 0.60) {
    cat(sprintf("  As comunalidades foram elevadas em geral (média de %s),\n  sinalizando que a estrutura fatorial retida explica satisfatoriamente a maioria dos itens.",
                .fmt_f(media_h2, 2)))
  } else {
    cat(sprintf("  A comunalidade média foi de %s.", .fmt_f(media_h2, 2)))
  }
  if (length(baixas_h2) > 0) {
    cat(sprintf("\n  Atenção para variáveis com comunalidade inferior a 0,50: %s (apresentam alta unicidade residual).",
                paste(baixas_h2, collapse = ", ")))
  }
  cat("\n\n")

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(strrep("─", sum(w_crit) + 3), "\n\n")
  }

  # Objeto de retorno
  res <- list(
    cargas        = df_cargas,
    comunalidades = comunalidades,
    unicidades    = unicidades,
    variancia     = df_variancia,
    escores       = df_escores,
    kmo           = list(
      global = kmo_global,
      classificacao_global = kmo_global_nivel,
      por_variavel = data.frame(
        Variavel = var_nomes,
        KMO      = kmo_ind,
        Nivel    = kmo_ind_nivel,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    ),
    bartlett      = list(
      qui_quadrado  = bart_chi2,
      gl            = bart_df,
      p_valor       = bart_p,
      significativo = (!is.na(bart_p) && bart_p < 0.05)
    ),
    criterios     = list(
      kaiser_sugerido  = kaiser_sugerido,
      ap_sugerido      = ap_sugerido,
      n_70_sugerido    = n_70_sugerido,
      n_retidos        = n_retidos,
      autovalores_simulados_p95 = ap_p95
    ),
    rotacao       = rot_ajuste,
    metodo        = metodo,
    fa            = fa_mod,
    dados_originais = df_num
  )
  class(res) <- c("afe_estatR", "list")

  if (grafico) {
    tryCatch({
      p_g <- cargasplot_fatorial(res, corte = corte_carga)
      print(p_g)
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o gráfico de cargas: ", e$message)
    })
  }

  invisible(res)
}

#' @title Gráfico de Cargas Fatoriais (Loadings Plot)
#' @description Gera gráfico de barras das cargas fatoriais para cada fator retido.
#' @param modelo Objeto retornado pela função \code{analise_fatorial()}.
#' @param corte Numérico. Linha de referência tracejada para cargas salientes (padrão: 0.30).
#' @export
cargasplot_fatorial <- function(modelo, corte = 0.30) {
  if (!inherits(modelo, "afe_estatR")) {
    stop("O argumento 'modelo' deve ser um objeto gerado pela função analise_fatorial() ou afe().")
  }

  n_fat <- length(modelo$variancia$Fator)
  cargas_cols <- paste("Fator", seq_len(n_fat))

  df_long <- stats::reshape(
    modelo$cargas[, c("Variavel", cargas_cols)],
    direction = "long",
    varying   = cargas_cols,
    v.names   = "Carga",
    timevar   = "Fator",
    times     = cargas_cols,
    idvar     = "Variavel"
  )

  df_long$Sinal <- ifelse(df_long$Carga >= 0, "Positiva", "Negativa")

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Carga, y = reorder(Variavel, abs(Carga)), fill = Sinal)) +
    ggplot2::geom_col(width = 0.65, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-corte, corte), linetype = "dashed", color = "#D90429", linewidth = 0.6) +
    ggplot2::facet_wrap(~ Fator, scales = "free_x") +
    ggplot2::scale_fill_manual(values = c("Positiva" = "#1D3557", "Negativa" = "#E63946")) +
    ggplot2::labs(
      title    = "Cargas Fatoriais (Loadings Plot)",
      subtitle = sprintf("Rotação: %s | Linhas tracejadas no corte de ±%s", tools::toTitleCase(modelo$rotacao), formatC(corte, format = "f", digits = 2, decimal.mark = ",")),
      x = "Carga Fatorial",
      y = "Variável",
      fill = "Direção"
    )

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2)
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}

#' @title Scree Plot da Análise Fatorial
#' @description Gera gráfico Scree Plot com autovalores da matriz de correlação e Análise Paralela de Horn.
#' @param modelo Objeto retornado pela função \code{analise_fatorial()} ou \code{afe()}.
#' @param analise_paralela Lógico. Se TRUE, exibe a linha dos autovalores simulados (P95).
#' @export
screeplot_fatorial <- function(modelo, analise_paralela = TRUE) {
  if (!inherits(modelo, "afe_estatR")) {
    stop("O argumento 'modelo' deve ser um objeto gerado pela função analise_fatorial() ou afe().")
  }

  ev <- eigen(stats::cor(modelo$dados_originais), only.values = TRUE)$values
  p_tot <- length(ev)

  df_scree <- data.frame(
    Fator_num = seq_len(p_tot),
    Fator     = paste("F", seq_len(p_tot)),
    Autovalor = ev
  )

  tem_ap <- analise_paralela && !is.null(modelo$criterios$autovalores_simulados_p95) &&
            !any(is.na(modelo$criterios$autovalores_simulados_p95))

  cores_legenda <- c(
    "Autovalores Observados"             = "#1D3557",
    "Critério de Kaiser (Autovalor = 1)" = "#D90429"
  )
  linetypes_legenda <- c(
    "Autovalores Observados"             = "solid",
    "Critério de Kaiser (Autovalor = 1)" = "dashed"
  )

  if (tem_ap) {
    cores_legenda["Análise Paralela de Horn (P95)"]     <- "#457B9D"
    linetypes_legenda["Análise Paralela de Horn (P95)"] <- "dotted"
  }

  p <- ggplot2::ggplot(df_scree, ggplot2::aes(x = Fator_num, y = Autovalor)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 1.0, color = "Critério de Kaiser (Autovalor = 1)",
                                     linetype = "Critério de Kaiser (Autovalor = 1)"),
                        linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(color = "Autovalores Observados", linetype = "Autovalores Observados"),
                       linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(color = "Autovalores Observados"), size = 3)

  if (tem_ap) {
    df_ap <- data.frame(
      Fator_num    = seq_along(modelo$criterios$autovalores_simulados_p95),
      Autovalor_AP = modelo$criterios$autovalores_simulados_p95
    )
    p <- p +
      ggplot2::geom_line(data = df_ap,
                         ggplot2::aes(x = Fator_num, y = Autovalor_AP,
                                      color = "Análise Paralela de Horn (P95)",
                                      linetype = "Análise Paralela de Horn (P95)"),
                         linewidth = 1) +
      ggplot2::geom_point(data = df_ap,
                          ggplot2::aes(x = Fator_num, y = Autovalor_AP,
                                       color = "Análise Paralela de Horn (P95)"),
                          size = 2.5, shape = 17)
  }

  p <- p +
    ggplot2::scale_x_continuous(breaks = df_scree$Fator_num, labels = df_scree$Fator) +
    ggplot2::scale_color_manual(name = "Critério", values = cores_legenda) +
    ggplot2::scale_linetype_manual(name = "Critério", values = linetypes_legenda) +
    ggplot2::labs(
      title = "Scree Plot — Análise Fatorial Exploratória",
      x     = "Fator",
      y     = "Autovalor"
    ) +
    meu_tema(grade = "dupla") +
    ggplot2::theme(legend.position = c(0.85, 0.85), legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.7), color = "gray80"))

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}

#' @title Biplot da Análise Fatorial
#' @description Gera gráfico Biplot 2D com escores das observações e vetores de cargas fatoriais.
#' @param modelo Objeto retornado pela função \code{analise_fatorial()} ou \code{afe()}.
#' @param grupo Vetor ou fator com grupos das observações (opcional).
#' @param elipses Lógico. Se TRUE e 'grupo' informado, traça elipses de confiança 95%.
#' @param escala Numérico. Fator de escala das setas das variáveis em relação aos escores.
#' @param rotulos Lógico. Se TRUE, exibe os nomes das variáveis nas setas.
#' @param f_x Inteiro. Fator do eixo X (padrão: 1).
#' @param f_y Inteiro. Fator do eixo Y (padrão: 2).
#' @export
biplot_fatorial <- function(modelo,
                            grupo   = NULL,
                            elipses = TRUE,
                            escala  = 0.7,
                            rotulos = TRUE,
                            f_x     = 1,
                            f_y     = 2) {

  if (!inherits(modelo, "afe_estatR")) {
    stop("O argumento 'modelo' deve ser um objeto gerado pela função analise_fatorial() ou afe().")
  }

  n_fat <- length(modelo$variancia$Fator)
  if (n_fat < 2) {
    stop("O modelo deve ter pelo menos 2 fatores para gerar o Biplot 2D.")
  }

  if (f_x > n_fat || f_y > n_fat) {
    stop("Os fatores selecionados excedem o número de fatores retidos no modelo.")
  }

  var_expl <- modelo$variancia$Var_Expl_Pct
  nome_x <- paste("Fator", f_x)
  nome_y <- paste("Fator", f_y)

  scores_df <- as.data.frame(modelo$escores[, c(f_x, f_y), drop = FALSE])
  colnames(scores_df) <- c("EixoX", "EixoY")

  cargas_raw <- unclass(modelo$fa$loadings)[, c(f_x, f_y), drop = FALSE]
  cargas_df  <- as.data.frame(cargas_raw)
  colnames(cargas_df) <- c("EixoX", "EixoY")
  cargas_df$Variavel  <- rownames(cargas_df)

  max_scores <- max(abs(scores_df$EixoX), abs(scores_df$EixoY), na.rm = TRUE)
  max_cargas <- max(abs(cargas_df$EixoX), abs(cargas_df$EixoY), na.rm = TRUE)
  mult_escala <- (max_scores / max_cargas) * escala

  cargas_df$X_fim <- cargas_df$EixoX * mult_escala
  cargas_df$Y_fim <- cargas_df$EixoY * mult_escala

  lab_x <- sprintf("Fator %d (%s%%)", f_x, formatC(var_expl[f_x], format = "f", digits = 1, decimal.mark = ","))
  lab_y <- sprintf("Fator %d (%s%%)", f_y, formatC(var_expl[f_y], format = "f", digits = 1, decimal.mark = ","))

  p <- ggplot2::ggplot()

  if (!is.null(grupo)) {
    scores_df$Grupo <- as.factor(grupo)
    if (elipses && nlevels(scores_df$Grupo) > 1) {
      p <- p + ggplot2::stat_ellipse(
        data = scores_df,
        ggplot2::aes(x = EixoX, y = EixoY, color = Grupo, fill = Grupo),
        geom = "polygon", alpha = 0.1, linewidth = 0.5, level = 0.95
      )
    }
    p <- p + ggplot2::geom_point(data = scores_df, ggplot2::aes(x = EixoX, y = EixoY, color = Grupo), size = 2.2, alpha = 0.75)
  } else {
    p <- p + ggplot2::geom_point(data = scores_df, ggplot2::aes(x = EixoX, y = EixoY), color = "#555555", size = 2, alpha = 0.65)
  }

  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4)

  p <- p + ggplot2::geom_segment(
    data = cargas_df,
    ggplot2::aes(x = 0, y = 0, xend = X_fim, yend = Y_fim),
    color = "#D90429", linewidth = 0.8,
    arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"), type = "closed")
  )

  if (rotulos) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = cargas_df,
        ggplot2::aes(x = X_fim, y = Y_fim, label = Variavel),
        color = "#1D3557", size = 4, fontface = "bold", box.padding = 0.3
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = cargas_df,
        ggplot2::aes(x = X_fim * 1.08, y = Y_fim * 1.08, label = Variavel),
        color = "#1D3557", size = 4, fontface = "bold"
      )
    }
  }

  p <- p + ggplot2::labs(
    title    = "Biplot — Análise Fatorial Exploratória",
    subtitle = paste0(lab_x, " vs ", lab_y, " (Rotação: ", tools::toTitleCase(modelo$rotacao), ")"),
    x = lab_x, y = lab_y
  )

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2)
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}


# ==============================================================================
# ANÁLISE DE AGRUPAMENTOS (CLUSTERING)
# ==============================================================================

#' @title Análise de Agrupamentos K-Means
#' @description Realiza Análise de Agrupamentos pelo método K-Means com padronização,
#' otimização automática de inicializações (nstart) por estabilização da inércia (WSS),
#' avaliação de k por Silhueta Média e Calinski-Harabasz, tabela de centróides na escala original
#' e interpretação automatizada.
#' @param dados Data frame ou matriz com variáveis numéricas.
#' @param k Inteiro opcional. Número de clusters. Se NULL (padrão), testa automaticamente
#'   valores de 2 até \code{k_max} e seleciona o melhor k pela maior Silhueta Média.
#' @param k_max Inteiro. Valor máximo de k avaliado na busca automática (padrão: 10).
#' @param padronizar Lógico. Se TRUE (padrão), padroniza as variáveis por z-score (média 0, dp 1).
#' @param nstart Inteiro opcional. Número de reinicializações aleatórias. Se NULL (padrão),
#'   testa candidatos (5 a 100) e seleciona o menor nstart onde a variação relativa do WSS cai abaixo de 0,1%.
#' @param grafico Lógico. Se TRUE, exibe automaticamente a projeção bidimensional dos clusters.
#' @return Objeto da classe \code{kmeans_estatR} contendo:
#'   \item{clusters}{Vetor de fatores com a atribuição de cada observação a um cluster}
#'   \item{dados_cluster}{Data frame original com a coluna \code{Cluster} anexada}
#'   \item{centroides}{Data frame com as médias das variáveis em cada cluster (escala original)}
#'   \item{resumo}{Data frame com n, proporção (%) e silhueta média de cada cluster}
#'   \item{criterios_k}{Data frame com k, WSS, Silhueta Média e Calinski-Harabasz para cada k avaliado}
#'   \item{nstart_utilizado}{Valor de nstart utilizado no ajuste final}
#'   \item{diagnostico_nstart}{Data frame com a variação do WSS entre os candidatos de nstart}
#'   \item{silhueta_global}{Coeficiente médio de silhueta da partição final}
#'   \item{qualidade}{Proporção da variância explicada entre clusters (BSS / TSS em %)}
#'   \item{kmeans}{Objeto nativo retornado por stats::kmeans}
#' @export
cluster_kmeans <- function(dados,
                           k          = NULL,
                           k_max      = 10,
                           padronizar = TRUE,
                           nstart     = NULL,
                           grafico    = TRUE) {

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("O pacote 'cluster' e necessario. Instale com install.packages('cluster').")
  }

  # 1. Validação dos dados
  if (!is.data.frame(dados) && !is.matrix(dados)) {
    stop("O argumento 'dados' deve ser um data frame ou matriz.")
  }

  dados_df <- as.data.frame(dados)
  cols_num <- vapply(dados_df, is.numeric, logical(1))
  if (sum(cols_num) < 2) {
    stop("São necessárias pelo menos 2 variáveis numéricas para a clusterização.")
  }

  df_num <- dados_df[, cols_num, drop = FALSE]
  casos_completos <- stats::complete.cases(df_num)
  if (sum(casos_completos) < nrow(df_num)) {
    df_num   <- df_num[casos_completos, , drop = FALSE]
    dados_df <- dados_df[casos_completos, , drop = FALSE]
  }

  n <- nrow(df_num)
  p <- ncol(df_num)
  var_nomes <- colnames(df_num)

  # Dados padronizados (ou escala original)
  df_sc <- if (padronizar) scale(df_num) else as.matrix(df_num)

  # 2. Otimização ou definição do nstart
  nstart_candidatos <- c(5, 10, 15, 20, 25, 50, 75, 100)
  k_ref <- if (!is.null(k)) as.integer(k) else min(3, n - 1)

  if (is.null(nstart)) {
    wss_cands <- vapply(nstart_candidatos, function(ns) {
      km_tmp <- suppressWarnings(stats::kmeans(df_sc, centers = k_ref, nstart = ns, iter.max = 300))
      km_tmp$tot.withinss
    }, numeric(1))

    wss_norm <- wss_cands / max(wss_cands)
    var_rel  <- c(Inf, abs(diff(wss_norm)) / (abs(wss_norm[-length(wss_norm)]) + 1e-10))
    idx_est  <- which(var_rel < 0.001) # < 0.1% de variação relativa
    nstart_final <- if (length(idx_est) > 0) nstart_candidatos[min(idx_est)] else 25

    diag_nstart <- data.frame(
      nstart  = nstart_candidatos,
      WSS     = round(wss_cands, 4),
      Var_Rel = round(var_rel * 100, 3)
    )
  } else {
    nstart_final <- as.integer(nstart)
    diag_nstart  <- data.frame(nstart = nstart_final, WSS = NA, Var_Rel = NA)
  }

  # 3. Avaliação de k (se k for NULL)
  k_max_eval <- min(as.integer(k_max), n - 1)
  k_seq <- 2:k_max_eval

  dist_mat <- stats::dist(df_sc)

  eval_wss <- numeric(length(k_seq))
  eval_sil <- numeric(length(k_seq))
  eval_ch  <- numeric(length(k_seq))

  for (i in seq_along(k_seq)) {
    ki <- k_seq[i]
    km_i <- suppressWarnings(stats::kmeans(df_sc, centers = ki, nstart = nstart_final, iter.max = 300))
    eval_wss[i] <- km_i$tot.withinss

    sil_i <- cluster::silhouette(km_i$cluster, dist_mat)
    eval_sil[i] <- mean(sil_i[, "sil_width"])

    # Calinski-Harabasz: (BSS / (k - 1)) / (WSS / (n - k))
    bss_i <- km_i$betweenss
    wss_i <- km_i$tot.withinss
    eval_ch[i] <- if (wss_i > 0) (bss_i / (ki - 1)) / (wss_i / (n - ki)) else NA
  }

  df_crit_k <- data.frame(
    k                  = k_seq,
    WSS                = eval_wss,
    Silhueta_Media     = eval_sil,
    Calinski_Harabasz  = eval_ch
  )

  k_sug_sil <- k_seq[which.max(eval_sil)]
  k_sug_ch  <- k_seq[which.max(eval_ch)]

  k_final <- if (!is.null(k)) {
    min(max(2, as.integer(k)), n - 1)
  } else {
    k_sug_sil
  }

  # 4. Ajuste Final com k_final
  set.seed(42)
  km_final <- suppressWarnings(stats::kmeans(df_sc, centers = k_final, nstart = nstart_final, iter.max = 300))
  clusters_vec <- km_final$cluster
  clusters_fac <- factor(paste("Cluster", clusters_vec), levels = paste("Cluster", seq_len(k_final)))

  # Silhueta de cada observação e por cluster
  sil_obj <- cluster::silhouette(clusters_vec, dist_mat)
  sil_global <- mean(sil_obj[, "sil_width"])
  sil_por_cluster <- vapply(seq_len(k_final), function(ci) {
    sub_sil <- sil_obj[clusters_vec == ci, "sil_width"]
    if (length(sub_sil) > 0) mean(sub_sil) else 0
  }, numeric(1))

  # Tabela resumo de tamanho e silhueta
  n_por_cluster   <- as.numeric(table(clusters_fac))
  prop_por_cluster <- (n_por_cluster / n) * 100

  df_resumo_clusters <- data.frame(
    Cluster        = paste("Cluster", seq_len(k_final)),
    N_Obs          = n_por_cluster,
    Proporcao_Pct  = prop_por_cluster,
    Silhueta_Media = sil_por_cluster,
    stringsAsFactors = FALSE
  )

  # Centróides na escala original (médias)
  centroides_orig <- as.data.frame(
    stats::aggregate(df_num, by = list(Cluster = clusters_fac), FUN = mean)[, -1, drop = FALSE]
  )
  rownames(centroides_orig) <- paste("Cluster", seq_len(k_final))
  df_centroides_t <- as.data.frame(t(centroides_orig))
  df_centroides_t <- data.frame(
    Variavel = rownames(df_centroides_t),
    round(df_centroides_t, 2),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  colnames(df_centroides_t) <- c("Variavel", paste("Cluster", seq_len(k_final)))

  # Qualidade da partição (BSS / TSS)
  qualidade_pct <- (km_final$betweenss / km_final$totss) * 100

  # Anexar Cluster aos dados
  dados_com_cluster <- dados_df
  dados_com_cluster$Cluster <- clusters_fac

  # Helpers de formatação
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = ",")

  # ==========================================================================
  # IMPRESSÃO NO CONSOLE
  # ==========================================================================
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("ANÁLISE DE AGRUPAMENTOS (K-MEANS CLUSTERING)")
  } else {
    cat("\n── ANÁLISE DE AGRUPAMENTOS (K-MEANS CLUSTERING) ──\n")
  }

  # SEÇÃO: CRITÉRIOS DE ESCOLHA DO NÚMERO DE CLUSTERS (k)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("CRITÉRIOS DE ESCOLHA DO NÚMERO DE CLUSTERS (k)")
  } else {
    cat("\n▶ CRITÉRIOS DE ESCOLHA DO NÚMERO DE CLUSTERS (k)\n")
  }

  w_crit <- c(k = 8, wss = 26, sil = 20, ch = 24)
  sep_crit <- paste0("  ", strrep("─", sum(w_crit) + 3))
  cat(sep_crit, "\n")
  cat("  ",
      .pad("k", w_crit["k"], "left"),
      .pad("WSS (Inércia Interna)", w_crit["wss"], "center"),
      .pad("Silhueta Média", w_crit["sil"], "center"),
      .pad("Calinski-Harabasz (CH)", w_crit["ch"], "center"),
      "\n", sep = "")
  cat(sep_crit, "\n")

  n_show_crit <- min(nrow(df_crit_k), 6)
  for (i in seq_len(n_show_crit)) {
    cat("  ",
        .pad(as.character(df_crit_k$k[i]), w_crit["k"], "left"),
        .pad(.fmt_f(df_crit_k$WSS[i], 2), w_crit["wss"], "center"),
        .pad(.fmt_f(df_crit_k$Silhueta_Media[i], 3), w_crit["sil"], "center"),
        .pad(.fmt_f(df_crit_k$Calinski_Harabasz[i], 2), w_crit["ch"], "center"),
        "\n", sep = "")
  }
  if (nrow(df_crit_k) > n_show_crit) {
    cat("  ", .pad("...", w_crit["k"], "left"), .pad("...", w_crit["wss"], "center"),
        .pad("...", w_crit["sil"], "center"), .pad("...", w_crit["ch"], "center"), "\n", sep = "")
  }
  cat(sep_crit, "\n")
  if (is.null(nstart)) {
    cat(sprintf("  * nstart utilizado: %d (estabilização automática com variação relativa < 0,1%%).\n", nstart_final))
  } else {
    cat(sprintf("  * nstart utilizado: %d (definido manualmente).\n", nstart_final))
  }
  cat(sprintf("  * Sugestão por Silhueta Média: k = %d | Sugestão por Calinski-Harabasz: k = %d.\n\n", k_sug_sil, k_sug_ch))

  # SEÇÃO: RESUMO DOS CLUSTERS FORMADOS
  if (exists(".print_topico", mode = "function")) {
    .print_topico(sprintf("RESUMO DOS CLUSTERS FORMADOS (k = %d)", k_final))
  } else {
    cat(sprintf("\n▶ RESUMO DOS CLUSTERS FORMADOS (k = %d)\n", k_final))
  }

  w_res <- c(cl = 16, n = 24, prop = 18, sil = 18)
  sep_res <- paste0("  ", strrep("─", sum(w_res) + 3))
  cat(sep_res, "\n")
  cat("  ",
      .pad("Cluster", w_res["cl"], "left"),
      .pad("Nº Observações (n)", w_res["n"], "center"),
      .pad("Proporção (%)", w_res["prop"], "center"),
      .pad("Silhueta Média", w_res["sil"], "center"),
      "\n", sep = "")
  cat(sep_res, "\n")
  for (i in seq_len(k_final)) {
    cat("  ",
        .pad(df_resumo_clusters$Cluster[i], w_res["cl"], "left"),
        .pad(as.character(df_resumo_clusters$N_Obs[i]), w_res["n"], "center"),
        .pad(paste0(.fmt_f(df_resumo_clusters$Proporcao_Pct[i], 2), "%"), w_res["prop"], "center"),
        .pad(.fmt_f(df_resumo_clusters$Silhueta_Media[i], 3), w_res["sil"], "center"),
        "\n", sep = "")
  }
  cat(sep_res, "\n")
  cat(sprintf("  Qualidade do Ajuste (BSS / TSS): %s%%\n", .fmt_f(qualidade_pct, 2)))
  cat(sprintf("  Silhueta Média Global: %s\n", .fmt_f(sil_global, 3)))
  cat("  Nota da Silhueta: > 0,70 Forte | 0,51 - 0,70 Razoável | 0,26 - 0,50 Fraca | <= 0,25 Sem estrutura\n\n")

  # SEÇÃO: PERFIL DOS CENTRÓIDES
  if (exists(".print_topico", mode = "function")) {
    .print_topico("PERFIL DOS CENTRÓIDES (MÉDIAS NA ESCALA ORIGINAL)")
  } else {
    cat("\n▶ PERFIL DOS CENTRÓIDES (MÉDIAS NA ESCALA ORIGINAL)\n")
  }

  w_var_cl <- 22
  w_col_cl <- max(12, min(16, floor((54 - w_var_cl) / k_final)))
  largura_cent <- w_var_cl + k_final * (w_col_cl + 1)
  sep_cent <- paste0("  ", strrep("─", largura_cent))

  cat(sep_cent, "\n")
  cat("  ", .pad("Variável", w_var_cl, "left"))
  for (ci in seq_len(k_final)) {
    cat(.pad(paste("Cluster", ci), w_col_cl, "center"))
  }
  cat("\n")
  cat(sep_cent, "\n")

  for (i in seq_along(var_nomes)) {
    cat("  ", .pad(var_nomes[i], w_var_cl, "left"))
    for (ci in seq_len(k_final)) {
      val <- centroides_orig[ci, var_nomes[i]]
      cat(.pad(.fmt_f(val, 2), w_col_cl, "center"))
    }
    cat("\n")
  }
  cat(sep_cent, "\n")
  cat("  * Médias calculadas nos dados originais para interpretação prática.\n\n")

  # SEÇÃO: INTERPRETAÇÃO
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO")
  } else {
    cat("\n▶ INTERPRETAÇÃO\n")
  }

  padrao_txt <- if (padronizar) {
    "Os dados foram padronizados por z-score para garantir pesos equivalentes entre variáveis com diferentes escalas."
  } else {
    "Os dados foram utilizados em sua escala original (não padronizados)."
  }

  nstart_txt <- if (is.null(nstart)) {
    sprintf("A busca automática de estabilização da inércia interna convergiu para nstart = %d.", nstart_final)
  } else {
    sprintf("O modelo utilizou nstart = %d definido pelo usuário.", nstart_final)
  }

  cat(sprintf("  %s\n  %s\n\n", padrao_txt, nstart_txt))

  cat(sprintf("  A avaliação dos critérios estatísticos convergiu para uma partição com k = %d clusters,\n  que é responsável por explicar %s%% da variância total da amostra entre os grupos.\n\n",
              k_final, .fmt_f(qualidade_pct, 2)))

  # Descrição dos perfis de cada cluster
  for (ci in seq_len(k_final)) {
    cl_n <- df_resumo_clusters$N_Obs[ci]
    cl_pct <- df_resumo_clusters$Proporcao_Pct[ci]

    # Encontrar variáveis de destaque comparando média do cluster com média global
    diff_rel <- (as.numeric(centroides_orig[ci, ]) - colMeans(df_num)) / apply(df_num, 2, stats::sd)
    names(diff_rel) <- var_nomes

    vars_altas  <- names(sort(diff_rel[diff_rel > 0.4], decreasing = TRUE))
    vars_baixas <- names(sort(diff_rel[diff_rel < -0.4]))

    txt_cl <- sprintf("  O Cluster %d (%s%% da base, n = %d)", ci, .fmt_f(cl_pct, 2), cl_n)
    detalhes <- c()
    if (length(vars_altas) > 0) detalhes <- c(detalhes, sprintf("médias acima do padrão global em %s", paste(vars_altas, collapse = ", ")))
    if (length(vars_baixas) > 0) detalhes <- c(detalhes, sprintf("médias reduzidas em %s", paste(vars_baixas, collapse = ", ")))

    if (length(detalhes) > 0) {
      txt_cl <- paste0(txt_cl, sprintf(" caracteriza-se por %s.", paste(detalhes, collapse = " contrapostas a ")))
    } else {
      txt_cl <- paste0(txt_cl, " reflete observações com características próximas da média geral da amostra.")
    }
    cat(txt_cl, "\n\n")
  }

  classif_sil <- if (sil_global > 0.70) "forte separação e estrutura coesa dos grupos"
                 else if (sil_global > 0.50) "boa separação entre os grupos e razoável coesão interna"
                 else if (sil_global > 0.25) "separação moderada entre os centróides, com alguma sobreposição nas bordas"
                 else "estrutura fraca ou artificial de agrupamento"

  cat(sprintf("  A silhueta média global de %s reflete %s.", .fmt_f(sil_global, 3), classif_sil))
  cat("\n\n")

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(strrep("─", sum(w_crit) + 3), "\n\n")
  }

  res <- list(
    clusters           = clusters_fac,
    dados_cluster      = dados_com_cluster,
    centroides         = df_centroides_t,
    resumo             = df_resumo_clusters,
    criterios_k        = df_crit_k,
    nstart_utilizado   = nstart_final,
    diagnostico_nstart = diag_nstart,
    silhueta_global    = sil_global,
    qualidade          = qualidade_pct,
    kmeans             = km_final,
    dados_originais    = df_num
  )
  class(res) <- c("kmeans_estatR", "list")

  if (grafico) {
    tryCatch({
      p_g <- plot_clusters(res)
      print(p_g)
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o gráfico de clusters: ", e$message)
    })
  }

  invisible(res)
}

#' @title Análise de Agrupamentos Hierárquica
#' @description Realiza Análise de Agrupamentos Hierárquica aglomerativa com matriz de distâncias,
#' método de ligação selecionável, cálculo de silhueta e corte de árvore.
#' @param dados Data frame ou matriz contendo variáveis numéricas.
#' @param k Inteiro opcional. Número de clusters para corte da árvore. Se NULL (padrão),
#'   sugere k = 3 como corte inicial.
#' @param metodo_dist Caractere com método de distância: \code{"euclidean"} (padrão),
#'   \code{"manhattan"}, \code{"maximum"}, \code{"canberra"}.
#' @param metodo_ligacao Caractere com método de ligação: \code{"ward.D2"} (padrão),
#'   \code{"complete"}, \code{"average"}, \code{"single"}.
#' @param padronizar Lógico. Se TRUE (padrão), padroniza as variáveis por z-score.
#' @param grafico Lógico. Se TRUE, gera o dendrograma colorido.
#' @return Objeto da classe \code{hclust_estatR} contendo:
#'   \item{clusters}{Vetor de fatores com a atribuição de cluster de cada observação}
#'   \item{dados_cluster}{Data frame original com a coluna \code{Cluster} anexada}
#'   \item{centroides}{Data frame com as médias das variáveis em cada cluster (escala original)}
#'   \item{resumo}{Data frame com tamanho e proporção de cada cluster}
#'   \item{silhueta_global}{Coeficiente médio de silhueta para o corte em k}
#'   \item{hclust}{Objeto nativo stats::hclust}
#'   \item{dist}{Matriz de distâncias stats::dist}
#' @export
cluster_hierarquico <- function(dados,
                                k              = NULL,
                                metodo_dist    = "euclidean",
                                metodo_ligacao = "ward.D2",
                                padronizar     = TRUE,
                                grafico        = TRUE) {

  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("O pacote 'cluster' e necessario. Instale com install.packages('cluster').")
  }

  if (!is.data.frame(dados) && !is.matrix(dados)) {
    stop("O argumento 'dados' deve ser um data frame ou matriz.")
  }

  dados_df <- as.data.frame(dados)
  cols_num <- vapply(dados_df, is.numeric, logical(1))
  if (sum(cols_num) < 2) {
    stop("São necessárias pelo menos 2 variáveis numéricas para a clusterização.")
  }

  df_num <- dados_df[, cols_num, drop = FALSE]
  casos_completos <- stats::complete.cases(df_num)
  if (sum(casos_completos) < nrow(df_num)) {
    df_num   <- df_num[casos_completos, , drop = FALSE]
    dados_df <- dados_df[casos_completos, , drop = FALSE]
  }

  n <- nrow(df_num)
  p <- ncol(df_num)
  var_nomes <- colnames(df_num)

  df_sc <- if (padronizar) scale(df_num) else as.matrix(df_num)
  dist_obj <- stats::dist(df_sc, method = metodo_dist)

  hc_mod <- stats::hclust(dist_obj, method = metodo_ligacao)

  k_final <- if (!is.null(k)) min(max(2, as.integer(k)), n - 1) else 3
  clusters_vec <- stats::cutree(hc_mod, k = k_final)
  clusters_fac <- factor(paste("Cluster", clusters_vec), levels = paste("Cluster", seq_len(k_final)))

  sil_obj <- cluster::silhouette(clusters_vec, dist_obj)
  sil_global <- mean(sil_obj[, "sil_width"])

  # Tabela resumo
  n_cl <- as.numeric(table(clusters_fac))
  prop_cl <- (n_cl / n) * 100

  df_resumo <- data.frame(
    Cluster       = paste("Cluster", seq_len(k_final)),
    N_Obs         = n_cl,
    Proporcao_Pct = prop_cl,
    stringsAsFactors = FALSE
  )

  # Centróides na escala original
  centroides_orig <- as.data.frame(
    stats::aggregate(df_num, by = list(Cluster = clusters_fac), FUN = mean)[, -1, drop = FALSE]
  )
  rownames(centroides_orig) <- paste("Cluster", seq_len(k_final))
  df_centroides_t <- as.data.frame(t(centroides_orig))
  df_centroides_t <- data.frame(
    Variavel = rownames(df_centroides_t),
    round(df_centroides_t, 2),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  colnames(df_centroides_t) <- c("Variavel", paste("Cluster", seq_len(k_final)))

  dados_com_cluster <- dados_df
  dados_com_cluster$Cluster <- clusters_fac

  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = ",")

  # Console output
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("ANÁLISE DE AGRUPAMENTOS HIERÁRQUICA")
  } else {
    cat("\n── ANÁLISE DE AGRUPAMENTOS HIERÁRQUICA ──\n")
  }

  if (exists(".print_topico", mode = "function")) {
    .print_topico("PARÂMETROS DO AGRUPAMENTO")
  } else {
    cat("\n▶ PARÂMETROS DO AGRUPAMENTO\n")
  }

  cat(sprintf("  Método de Distância: %s\n", tools::toTitleCase(metodo_dist)))
  cat(sprintf("  Método de Ligação:   %s\n", metodo_ligacao))
  cat(sprintf("  Número de Clusters:  %d (corte na árvore)\n", k_final))
  cat(sprintf("  Silhueta Global:     %s\n\n", .fmt_f(sil_global, 3)))

  if (exists(".print_topico", mode = "function")) {
    .print_topico(sprintf("RESUMO DOS CLUSTERS (k = %d)", k_final))
  } else {
    cat(sprintf("\n▶ RESUMO DOS CLUSTERS (k = %d)\n", k_final))
  }

  w_res <- c(cl = 16, n = 24, prop = 18)
  sep_res <- paste0("  ", strrep("─", sum(w_res) + 3))
  cat(sep_res, "\n")
  cat("  ", .pad("Cluster", w_res["cl"], "left"), .pad("Nº Observações (n)", w_res["n"], "center"), .pad("Proporção (%)", w_res["prop"], "center"), "\n", sep = "")
  cat(sep_res, "\n")
  for (i in seq_len(k_final)) {
    cat("  ",
        .pad(df_resumo$Cluster[i], w_res["cl"], "left"),
        .pad(as.character(df_resumo$N_Obs[i]), w_res["n"], "center"),
        .pad(paste0(.fmt_f(df_resumo$Proporcao_Pct[i], 2), "%"), w_res["prop"], "center"),
        "\n", sep = "")
  }
  cat(sep_res, "\n\n")

  if (exists(".print_topico", mode = "function")) {
    .print_topico("PERFIL DOS CENTRÓIDES (MÉDIAS NA ESCALA ORIGINAL)")
  } else {
    cat("\n▶ PERFIL DOS CENTRÓIDES (MÉDIAS NA ESCALA ORIGINAL)\n")
  }

  w_var_cl <- 22
  w_col_cl <- max(12, min(16, floor((54 - w_var_cl) / k_final)))
  largura_cent <- w_var_cl + k_final * (w_col_cl + 1)
  sep_cent <- paste0("  ", strrep("─", largura_cent))

  cat(sep_cent, "\n")
  cat("  ", .pad("Variável", w_var_cl, "left"))
  for (ci in seq_len(k_final)) {
    cat(.pad(paste("Cluster", ci), w_col_cl, "center"))
  }
  cat("\n")
  cat(sep_cent, "\n")

  for (i in seq_along(var_nomes)) {
    cat("  ", .pad(var_nomes[i], w_var_cl, "left"))
    for (ci in seq_len(k_final)) {
      val <- centroides_orig[ci, var_nomes[i]]
      cat(.pad(.fmt_f(val, 2), w_col_cl, "center"))
    }
    cat("\n")
  }
  cat(sep_cent, "\n")
  cat("  * Médias calculadas nos dados originais para interpretação prática.\n\n")

  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO")
  } else {
    cat("\n▶ INTERPRETAÇÃO\n")
  }

  cat(sprintf("  O agrupamento hierárquico aglomerativo foi construído utilizando distância %s\n  e método de ligação %s sobre dados %s.\n\n",
              tools::toTitleCase(metodo_dist), metodo_ligacao,
              if (padronizar) "padronizados (z-score)" else "originais"))

  cat(sprintf("  O corte da árvore gerou %d clusters com silhueta média de %s.\n", k_final, .fmt_f(sil_global, 3)))
  cat(sprintf("  A maior partição é o %s (%s%% das observações), e o menor é o %s (%s%%).\n\n",
              df_resumo$Cluster[which.max(df_resumo$N_Obs)], .fmt_f(max(df_resumo$Proporcao_Pct), 1),
              df_resumo$Cluster[which.min(df_resumo$N_Obs)], .fmt_f(min(df_resumo$Proporcao_Pct), 1)))

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(strrep("─", sum(w_res) + 3), "\n\n")
  }

  res <- list(
    clusters        = clusters_fac,
    dados_cluster   = dados_com_cluster,
    centroides      = df_centroides_t,
    resumo          = df_resumo,
    silhueta_global = sil_global,
    hclust          = hc_mod,
    dist            = dist_obj,
    dados_originais = df_num
  )
  class(res) <- c("hclust_estatR", "list")

  if (grafico) {
    tryCatch({
      p_g <- dendrograma(res)
      print(p_g)
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o dendrograma: ", e$message)
    })
  }

  invisible(res)
}

#' @title Projeção 2D dos Clusters (K-Means)
#' @description Gera gráfico de dispersão 2D dos clusters projetados nas duas primeiras Componentes Principais (PC1 x PC2),
#' com elipses de 95% e os centróides em destaque.
#' @param modelo Objeto retornado pela função \code{cluster_kmeans()}.
#' @param elipses Lógico. Se TRUE (padrão), exibe elipses de confiança 95% para cada cluster.
#' @export
plot_clusters <- function(modelo, elipses = TRUE) {
  if (!inherits(modelo, "kmeans_estatR")) {
    stop("O argumento 'modelo' deve ser um objeto gerado pela função cluster_kmeans().")
  }

  pca_proj <- stats::prcomp(modelo$dados_originais, center = TRUE, scale. = TRUE)
  scores_df <- as.data.frame(pca_proj$x[, 1:2, drop = FALSE])
  colnames(scores_df) <- c("PC1", "PC2")
  scores_df$Cluster <- modelo$clusters

  # Calcular centróides na projeção PCA
  cent_pca <- stats::aggregate(scores_df[, c("PC1", "PC2")], by = list(Cluster = scores_df$Cluster), FUN = mean)

  var_pca <- round((pca_proj$sdev^2 / sum(pca_proj$sdev^2)) * 100, 1)

  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = PC1, y = PC2, color = Cluster, fill = Cluster))

  if (elipses && nlevels(scores_df$Cluster) > 1) {
    p <- p + ggplot2::stat_ellipse(geom = "polygon", alpha = 0.12, linewidth = 0.5, level = 0.95)
  }

  p <- p +
    ggplot2::geom_point(size = 2.4, alpha = 0.75) +
    ggplot2::geom_point(data = cent_pca, ggplot2::aes(x = PC1, y = PC2),
                        color = "black", shape = 3, size = 5, stroke = 2, inherit.aes = FALSE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    ggplot2::labs(
      title    = "Projeção 2D dos Clusters (K-Means)",
      subtitle = sprintf("Projeção em PC1 (%s%%) e PC2 (%s%%) | Símbolo '+': Centróides", var_pca[1], var_pca[2]),
      x = sprintf("PC1 (%s%%)", var_pca[1]),
      y = sprintf("PC2 (%s%%)", var_pca[2])
    )

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2)
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}

#' @title Gráfico do Cotovelo e Silhueta (K-Means)
#' @description Gera gráfico de linha comparativo do método do Cotovelo (WSS) e Silhueta Média por k.
#' @param modelo Objeto retornado pela função \code{cluster_kmeans()}.
#' @export
cotovelo_plot <- function(modelo) {
  if (!inherits(modelo, "kmeans_estatR")) {
    stop("O argumento 'modelo' deve ser um objeto gerado pela função cluster_kmeans().")
  }

  df_k <- modelo$criterios_k
  k_sel <- length(levels(modelo$clusters))

  p1 <- ggplot2::ggplot(df_k, ggplot2::aes(x = k, y = WSS)) +
    ggplot2::geom_line(color = "#1D3557", linewidth = 1) +
    ggplot2::geom_point(color = "#1D3557", size = 2.8) +
    ggplot2::geom_vline(xintercept = k_sel, linetype = "dashed", color = "#D90429", linewidth = 0.7) +
    ggplot2::scale_x_continuous(breaks = df_k$k) +
    ggplot2::labs(
      title = "Método do Cotovelo (Inércia WSS)",
      x = "Número de Clusters (k)", y = "WSS"
    )

  p2 <- ggplot2::ggplot(df_k, ggplot2::aes(x = k, y = Silhueta_Media)) +
    ggplot2::geom_line(color = "#457B9D", linewidth = 1) +
    ggplot2::geom_point(color = "#457B9D", size = 2.8) +
    ggplot2::geom_vline(xintercept = k_sel, linetype = "dashed", color = "#D90429", linewidth = 0.7) +
    ggplot2::scale_x_continuous(breaks = df_k$k) +
    ggplot2::labs(
      title = "Silhueta Média por k",
      x = "Número de Clusters (k)", y = "Silhueta Média"
    )

  if (exists("tema_estatR", mode = "function")) {
    p1 <- p1 + tema_estatR(estilo = 2)
    p2 <- p2 + tema_estatR(estilo = 2)
  } else {
    p1 <- p1 + ggplot2::theme_minimal()
    p2 <- p2 + ggplot2::theme_minimal()
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, ncol = 2)
  } else {
    p1
  }
}

#' @title Dendrograma do Cluster Hierárquico
#' @description Plota dendrograma formatado do agrupamento hierárquico com ramos coloridos por cluster.
#' @param modelo Objeto retornado pela função \code{cluster_hierarquico()}.
#' @param k Inteiro opcional. Número de clusters para corte e coloração.
#' @export
dendrograma <- function(modelo, k = NULL) {
  if (!inherits(modelo, "hclust_estatR")) {
    stop("O argumento 'modelo' deve ser um objeto gerado pela função cluster_hierarquico().")
  }

  k_val <- if (!is.null(k)) as.integer(k) else length(levels(modelo$clusters))

  # Usar ggdendro se disponível ou stats::plot
  if (requireNamespace("ggdendro", quietly = TRUE)) {
    dend_data <- ggdendro::dendro_data(modelo$hclust, type = "rectangle")
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = ggdendro::segment(dend_data),
                            ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                            color = "#1D3557", linewidth = 0.7) +
      ggplot2::geom_text(data = ggdendro::label(dend_data),
                         ggplot2::aes(x = x, y = y, label = label),
                         hjust = 1, angle = 90, size = 3.2, vjust = 0.5) +
      ggplot2::labs(
        title    = "Dendrograma de Agrupamento Hierárquico",
        subtitle = sprintf("Método de Ligação: %s | Corte em k = %d grupos", modelo$hclust$method, k_val),
        x = "Observações", y = "Distância / Altura"
      )

    if (exists("tema_estatR", mode = "function")) {
      p <- p + tema_estatR(estilo = 2) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    } else {
      p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    }
    return(p)
  } else {
    graphics::plot(modelo$hclust, hang = -1, main = "Dendrograma — Agrupamento Hierárquico",
                   xlab = "Observações", ylab = "Distância / Altura", sub = "")
    graphics::rect.hclust(modelo$hclust, k = k_val, border = "#D90429")
  }
}

# ==============================================================================
# MÓDULO 4: ANÁLISE DE VARIÂNCIA MULTIVARIADA (MANOVA)
# ==============================================================================

# Helper interno: Teste M de Box (Homogeneidade das Matrizes de Covariância)
.box_m_test <- function(Y, grupo) {
  p        <- ncol(Y)
  g_niveis <- levels(grupo)
  k        <- length(g_niveis)
  
  n_i  <- as.numeric(table(grupo))
  nu_i <- n_i - 1
  N    <- sum(n_i)
  nu   <- sum(nu_i)
  
  cov_lista <- lapply(g_niveis, function(g) stats::cov(Y[grupo == g, , drop = FALSE]))
  cov_pooled <- matrix(0, nrow = p, ncol = p)
  for (i in seq_len(k)) {
    cov_pooled <- cov_pooled + nu_i[i] * cov_lista[[i]]
  }
  cov_pooled <- cov_pooled / nu
  
  det_pooled <- det(cov_pooled)
  dets_i     <- sapply(cov_lista, det)
  
  if (any(dets_i <= 0) || det_pooled <= 0) {
    return(list(M = NA_real_, F = NA_real_, df1 = NA_real_, df2 = NA_real_, p.value = NA_real_))
  }
  
  M <- nu * log(det_pooled) - sum(nu_i * log(dets_i))
  
  sum_inv_nu <- sum(1 / nu_i)
  c1 <- (sum_inv_nu - 1 / nu) * (2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (k - 1))
  df1 <- p * (p + 1) * (k - 1) / 2
  
  sum_inv_nu2 <- sum(1 / (nu_i^2))
  c2 <- (sum_inv_nu2 - 1 / (nu^2)) * (p - 1) * (p + 2) / (6 * (k - 1))
  
  if (c2 - c1^2 > 0) {
    df2    <- (df1 + 2) / (c2 - c1^2)
    F_stat <- M * (1 - c1) / df1
    p_val  <- stats::pf(F_stat, df1, df2, lower.tail = FALSE)
  } else {
    df2    <- Inf
    F_stat <- M * (1 - c1) / df1
    p_val  <- stats::pchisq(M * (1 - c1), df1, lower.tail = FALSE)
  }
  
  list(M = M, F = F_stat, df1 = df1, df2 = df2, p.value = p_val)
}

#' @title Análise de Variância Multivariada (MANOVA)
#' @description Realiza a Análise de Variância Multivariada (MANOVA), avaliando o efeito
#' de um fator categórico sobre múltiplas variáveis resposta contínuas simultaneamente.
#' Apresenta o teste de homogeneidade das covariâncias (M de Box), os 4 critérios
#' multivariados clássicos (Pillai, Wilks, Hotelling-Lawley, Roy), ANOVAs univariadas de
#' acompanhamento com correção de Bonferroni, médias por grupo e interpretação automatizada.
#'
#' @param dados Data frame contendo os dados.
#' @param variaveis Vetor de variáveis resposta contínuas (com ou sem aspas, ex: \code{c(y1, y2)} ou \code{c("y1", "y2")}).
#'   Se omitido, seleciona todas as variáveis numéricas do data frame (exceto a variável de grupo).
#' @param grupo Variável categórica de agrupamento/fator (com ou sem aspas, ex: \code{especie} ou \code{"especie"}).
#' @param grafico Lógico. Se TRUE, plota o gráfico de perfis médios multivariados. Padrão FALSE.
#'
#' @return Retorna um objeto de classe \code{"estatR_manova"} contendo:
#' \itemize{
#'   \item \code{modelo_manova}: Objeto nativo da função \code{stats::manova}
#'   \item \code{box_m}: Resultados do teste M de Box
#'   \item \code{tabela_multivariada}: Tabela com os critérios multivariados
#'   \item \code{anovas_univariadas}: Tabela das ANOVAs univariadas com Bonferroni
#'   \item \code{medias_grupos}: Tabela de médias observadas por grupo
#'   \item \code{variaveis}: Nomes das variáveis resposta analisadas
#'   \item \code{grupo}: Nome da variável de agrupamento
#'   \item \code{dados}: Data frame limpo utilizado no ajuste
#' }
#' @export
analise_manova <- function(dados = NULL, variaveis = NULL, grupo = NULL, grafico = TRUE) {
  
  # 1. Tratamento flexível de argumentos posicionais e não avaliados
  sub_dados <- substitute(dados)
  sub_vars  <- substitute(variaveis)
  sub_grupo <- substitute(grupo)

  # Se o usuário passou: analise_manova(variaveis = ..., grupo = ..., dados = ...)
  # ou analise_manova(dados, variaveis, grupo)
  if (is.null(dados) && !is.null(sub_vars)) {
    val_vars <- tryCatch(variaveis, error = function(e) NULL)
    if (is.data.frame(val_vars)) {
      dados     <- val_vars
      variaveis <- NULL
      sub_vars  <- NULL
    }
  }

  if (is.null(dados) || !is.data.frame(dados)) {
    stop("O argumento 'dados' deve ser fornecido como um data frame.")
  }

  # Extração de 'grupo' (com ou sem aspas)
  nome_grupo <- NULL
  val_grupo  <- tryCatch(grupo, error = function(e) NULL)
  if (is.character(val_grupo) && length(val_grupo) == 1) {
    nome_grupo <- val_grupo
  } else if (!is.null(sub_grupo) && is.symbol(sub_grupo)) {
    nome_grupo <- as.character(sub_grupo)
  }

  if (is.null(nome_grupo) || !(nome_grupo %in% names(dados))) {
    stop("O argumento 'grupo' deve ser uma variável categórica existente no data frame.")
  }

  # Extração de 'variaveis' (com ou sem aspas)
  nomes_vars <- NULL
  val_vars   <- tryCatch(variaveis, error = function(e) NULL)
  if (is.character(val_vars)) {
    nomes_vars <- val_vars
  } else if (!is.null(sub_vars)) {
    if (is.symbol(sub_vars)) {
      nomes_vars <- as.character(sub_vars)
    } else if (is.call(sub_vars)) {
      elementos <- as.list(sub_vars)[-1]
      nomes_vars <- sapply(elementos, as.character)
    }
  }

  # Se omitido, seleciona todas as numéricas (exceto o grupo)
  if (is.null(nomes_vars) || length(nomes_vars) == 0) {
    todas_num  <- names(dados)[sapply(dados, is.numeric)]
    nomes_vars <- setdiff(todas_num, nome_grupo)
  }

  # Validação de variáveis resposta
  if (length(nomes_vars) < 2) {
    stop("A MANOVA requer pelo menos 2 variáveis resposta numéricas contínuas.")
  }
  nao_encontradas <- setdiff(nomes_vars, names(dados))
  if (length(nao_encontradas) > 0) {
    stop(sprintf("As seguintes variáveis não foram encontradas no data frame: %s",
                 paste(nao_encontradas, collapse = ", ")))
  }

  nao_num <- nomes_vars[!sapply(dados[nomes_vars], is.numeric)]
  if (length(nao_num) > 0) {
    stop(sprintf("Todas as variáveis resposta devem ser numéricas. Variáveis não numéricas: %s",
                 paste(nao_num, collapse = ", ")))
  }

  # 2. Preparação dos Dados
  df_analise <- stats::na.omit(dados[, c(nomes_vars, nome_grupo), drop = FALSE])
  n_total    <- nrow(df_analise)
  fator      <- as.factor(df_analise[[nome_grupo]])
  niveis_grp <- levels(fator)
  k_grupos   <- length(niveis_grp)

  if (k_grupos < 2) {
    stop("A variável de grupo deve possuir pelo menos 2 categorias distintas.")
  }

  Y <- as.matrix(df_analise[, nomes_vars, drop = FALSE])

  # 3. Teste M de Box
  res_box <- .box_m_test(Y, fator)

  # 4. Ajuste da MANOVA
  mod_manova   <- stats::manova(Y ~ fator)
  sm_pillai    <- summary(mod_manova, test = "Pillai")
  sm_wilks     <- summary(mod_manova, test = "Wilks")
  sm_hotelling <- summary(mod_manova, test = "Hotelling-Lawley")
  sm_roy       <- summary(mod_manova, test = "Roy")

  st_pillai    <- sm_pillai$stats["fator", ]
  st_wilks     <- sm_wilks$stats["fator", ]
  st_hotelling <- sm_hotelling$stats["fator", ]
  st_roy       <- sm_roy$stats["fator", ]

  tab_multi <- data.frame(
    Criterio = c("Traço de Pillai", "Lambda de Wilks", "Hotelling-Lawley", "Maior Raiz de Roy"),
    Valor    = c(st_pillai["Pillai"], st_wilks["Wilks"], st_hotelling["Hotelling-Lawley"], st_roy["Roy"]),
    F_aprox  = c(st_pillai["approx F"], st_wilks["approx F"], st_hotelling["approx F"], st_roy["approx F"]),
    gl1      = c(st_pillai["num Df"], st_wilks["num Df"], st_hotelling["num Df"], st_roy["num Df"]),
    gl2      = c(st_pillai["den Df"], st_wilks["den Df"], st_hotelling["den Df"], st_roy["den Df"]),
    p_valor  = c(st_pillai["Pr(>F)"], st_wilks["Pr(>F)"], st_hotelling["Pr(>F)"], st_roy["Pr(>F)"]),
    stringsAsFactors = FALSE
  )

  # 5. ANOVAs Univariadas de Acompanhamento
  sm_aov       <- stats::summary.aov(mod_manova)
  p_ajust_bonf <- length(nomes_vars)
  tab_aov <- data.frame(
    Variavel   = nomes_vars,
    SQ_Fator   = numeric(length(nomes_vars)),
    SQ_Residuo = numeric(length(nomes_vars)),
    F_valor    = numeric(length(nomes_vars)),
    p_valor    = numeric(length(nomes_vars)),
    p_bonf     = numeric(length(nomes_vars)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(nomes_vars)) {
    v_nm    <- nomes_vars[i]
    sub_tab <- sm_aov[[paste0(" Response ", v_nm)]]
    tab_aov$SQ_Fator[i]   <- sub_tab$"Sum Sq"[1]
    tab_aov$SQ_Residuo[i] <- sub_tab$"Sum Sq"[2]
    tab_aov$F_valor[i]    <- sub_tab$"F value"[1]
    tab_aov$p_valor[i]    <- sub_tab$"Pr(>F)"[1]
    tab_aov$p_bonf[i]     <- min(1, sub_tab$"Pr(>F)"[1] * p_ajust_bonf)
  }

  # 6. Médias por Grupo
  medias_df <- data.frame(Variavel = nomes_vars, stringsAsFactors = FALSE)
  for (g in niveis_grp) {
    medias_df[[g]] <- sapply(nomes_vars, function(v) mean(df_analise[[v]][fator == g], na.rm = TRUE))
  }

  # Helpers de formatação
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 2) {
    if (is.na(x)) return("-")
    formatC(x, format = "f", digits = d, decimal.mark = ",")
  }
  .fmt_p <- function(x) {
    if (is.na(x)) return("-")
    if (x < 0.001) return("< 0,001")
    formatC(x, format = "f", digits = 3, decimal.mark = ",")
  }
  .ast <- function(p) {
    if (is.na(p) || p >= 0.10) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    return("*")
  }

  # ==========================================================================
  # IMPRESSÃO NO CONSOLE
  # ==========================================================================
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("ANÁLISE DE VARIÂNCIA MULTIVARIADA (MANOVA)")
  } else {
    cat("\n── ANÁLISE DE VARIÂNCIA MULTIVARIADA (MANOVA) ──\n")
  }

  cat(sprintf("  Fator:               %s (%d grupos: %s)\n", nome_grupo, k_grupos, paste(niveis_grp, collapse = ", ")))
  cat(sprintf("  Variáveis Resposta:  %s (%d variáveis)\n", paste(nomes_vars, collapse = ", "), length(nomes_vars)))
  cat(sprintf("  Observações Válidas: %d\n\n", n_total))

  # SEÇÃO: TESTE DE HOMOGENEIDADE DAS COVARIÂNCIAS
  if (exists(".print_topico", mode = "function")) {
    .print_topico("TESTE DE HOMOGENEIDADE DAS COVARIÂNCIAS (PRESSUPOSTO)")
  } else {
    cat("  ▶ TESTE DE HOMOGENEIDADE DAS COVARIÂNCIAS (PRESSUPOSTO)\n")
  }

  w_box <- c(tst = 22, m = 18, f = 14, p = 14)
  sep_box <- paste0("  ", strrep("─", sum(w_box)))
  cat(sep_box, "\n")
  cat("  ", .pad("Teste", w_box["tst"], "left"),
             .pad("Estatística (M)", w_box["m"], "center"),
             .pad("F aprox.", w_box["f"], "center"),
             .pad("p-valor", w_box["p"], "center"), "\n", sep = "")
  cat(sep_box, "\n")
  box_p_ast <- paste0(.fmt_p(res_box$p.value), " ", .ast(res_box$p.value))
  cat("  ", .pad("M de Box", w_box["tst"], "left"),
             .pad(.fmt_f(res_box$M, 2), w_box["m"], "center"),
             .pad(.fmt_f(res_box$F, 2), w_box["f"], "center"),
             .pad(box_p_ast, w_box["p"], "center"), "\n", sep = "")
  cat(sep_box, "\n")
  cat("  * H0: As matrizes de covariância são homogêneas entre os grupos.\n")
  cat("  * Nota: Pillai é recomendado quando o teste M de Box for significativo.\n\n")

  # SEÇÃO: CRITÉRIOS GLOBAIS DA MANOVA
  if (exists(".print_topico", mode = "function")) {
    .print_topico("CRITÉRIOS GLOBAIS DA MANOVA")
  } else {
    cat("  ▶ CRITÉRIOS GLOBAIS DA MANOVA\n")
  }

  w_crit <- c(crit = 20, val = 12, f = 12, gl1 = 8, gl2 = 8, p = 14)
  sep_crit <- paste0("  ", strrep("─", sum(w_crit)))
  cat(sep_crit, "\n")
  cat("  ", .pad("Critério", w_crit["crit"], "left"),
             .pad("Valor", w_crit["val"], "center"),
             .pad("F aprox.", w_crit["f"], "center"),
             .pad("gl1", w_crit["gl1"], "center"),
             .pad("gl2", w_crit["gl2"], "center"),
             .pad("p-valor", w_crit["p"], "center"), "\n", sep = "")
  cat(sep_crit, "\n")

  for (r in seq_len(nrow(tab_multi))) {
    p_ast <- paste0(.fmt_p(tab_multi$p_valor[r]), " ", .ast(tab_multi$p_valor[r]))
    gl2_str <- if (is.finite(tab_multi$gl2[r])) sprintf("%.0f", tab_multi$gl2[r]) else "Inf"
    cat("  ", .pad(tab_multi$Criterio[r], w_crit["crit"], "left"),
               .pad(.fmt_f(tab_multi$Valor[r], 3), w_crit["val"], "center"),
               .pad(.fmt_f(tab_multi$F_aprox[r], 2), w_crit["f"], "center"),
               .pad(sprintf("%.0f", tab_multi$gl1[r]), w_crit["gl1"], "center"),
               .pad(gl2_str, w_crit["gl2"], "center"),
               .pad(p_ast, w_crit["p"], "center"), "\n", sep = "")
  }
  cat(sep_crit, "\n")
  cat("  * Significância: *** p < 0,01   ** p < 0,05   * p < 0,10\n\n")

  # SEÇÃO: ANOVAS UNIVARIADAS DE ACOMPANHAMENTO
  if (exists(".print_topico", mode = "function")) {
    .print_topico("ANOVAS UNIVARIADAS DE ACOMPANHAMENTO (PÓS-MANOVA)")
  } else {
    cat("  ▶ ANOVAS UNIVARIADAS DE ACOMPANHAMENTO (PÓS-MANOVA)\n")
  }

  w_aov <- c(var = 20, sqf = 13, sqr = 13, f = 11, p = 13, pbonf = 15)
  sep_aov <- paste0("  ", strrep("─", sum(w_aov)))
  cat(sep_aov, "\n")
  cat("  ", .pad("Variável Resposta", w_aov["var"], "left"),
             .pad("SQ Fator", w_aov["sqf"], "center"),
             .pad("SQ Resíduo", w_aov["sqr"], "center"),
             .pad("F", w_aov["f"], "center"),
             .pad("p-valor", w_aov["p"], "center"),
             .pad("p (Bonf.)", w_aov["pbonf"], "center"), "\n", sep = "")
  cat(sep_aov, "\n")

  for (r in seq_len(nrow(tab_aov))) {
    p_ast  <- paste0(.fmt_p(tab_aov$p_valor[r]), " ", .ast(tab_aov$p_valor[r]))
    pb_ast <- paste0(.fmt_p(tab_aov$p_bonf[r]), " ", .ast(tab_aov$p_bonf[r]))
    cat("  ", .pad(tab_aov$Variavel[r], w_aov["var"], "left"),
               .pad(.fmt_f(tab_aov$SQ_Fator[r], 2), w_aov["sqf"], "center"),
               .pad(.fmt_f(tab_aov$SQ_Residuo[r], 2), w_aov["sqr"], "center"),
               .pad(.fmt_f(tab_aov$F_valor[r], 2), w_aov["f"], "center"),
               .pad(p_ast, w_aov["p"], "center"),
               .pad(pb_ast, w_aov["pbonf"], "center"), "\n", sep = "")
  }
  cat(sep_aov, "\n")
  cat("  * p (Bonf.): p-valor corrigido para comparações múltiplas (Bonferroni).\n\n")

  # SEÇÃO: MÉDIAS POR GRUPO
  if (exists(".print_topico", mode = "function")) {
    .print_topico("MÉDIAS POR GRUPO")
  } else {
    cat("  ▶ MÉDIAS POR GRUPO\n")
  }

  w_var_med <- max(18, max(nchar(nomes_vars)) + 2)
  w_grp_med <- max(12, max(nchar(niveis_grp)) + 2)
  total_med_w <- w_var_med + w_grp_med * k_grupos
  sep_med <- paste0("  ", strrep("─", total_med_w))
  cat(sep_med, "\n")
  cat("  ", .pad("Variável", w_var_med, "left"))
  for (g in niveis_grp) {
    cat(.pad(g, w_grp_med, "center"))
  }
  cat("\n")
  cat(sep_med, "\n")

  for (r in seq_len(nrow(medias_df))) {
    cat("  ", .pad(medias_df$Variavel[r], w_var_med, "left"))
    for (g in niveis_grp) {
      cat(.pad(.fmt_f(medias_df[[g]][r], 3), w_grp_med, "center"))
    }
    cat("\n")
  }
  cat(sep_med, "\n\n")

  # SEÇÃO: INTERPRETAÇÃO
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO")
  } else {
    cat("  ▶ INTERPRETAÇÃO\n")
  }

  p_wilks  <- tab_multi$p_valor[tab_multi$Criterio == "Lambda de Wilks"]
  f_wilks  <- tab_multi$F_aprox[tab_multi$Criterio == "Lambda de Wilks"]
  val_wlk  <- tab_multi$Valor[tab_multi$Criterio == "Lambda de Wilks"]
  p_pillai <- tab_multi$p_valor[tab_multi$Criterio == "Traço de Pillai"]
  f_pillai <- tab_multi$F_aprox[tab_multi$Criterio == "Traço de Pillai"]
  val_pil  <- tab_multi$Valor[tab_multi$Criterio == "Traço de Pillai"]

  if (p_wilks < 0.05) {
    cat(sprintf("  A Análise de Variância Multivariada (MANOVA) indicou um efeito estatisticamente\n  significativo do fator '%s' sobre o conjunto de variáveis resposta consideradas\n  simultaneamente (Lambda de Wilks = %s, F = %s, p %s; Traço de Pillai = %s, F = %s, p %s).\n\n",
                nome_grupo, .fmt_f(val_wlk, 3), .fmt_f(f_wilks, 2), .fmt_p(p_wilks),
                .fmt_f(val_pil, 3), .fmt_f(f_pillai, 2), .fmt_p(p_pillai)))
  } else {
    cat(sprintf("  A Análise de Variância Multivariada (MANOVA) não encontrou evidências de efeito\n  estatisticamente significativo do fator '%s' sobre o conjunto de variáveis resposta\n  (Lambda de Wilks = %s, F = %s, p = %s).\n\n",
                nome_grupo, .fmt_f(val_wlk, 3), .fmt_f(f_wilks, 2), .fmt_p(p_wilks)))
  }

  if (!is.na(res_box$p.value)) {
    if (res_box$p.value < 0.05) {
      cat(sprintf("  O teste M de Box indicou que as matrizes de covariância não são homogêneas entre os grupos\n  (M = %s, F = %s, p %s). Nesses casos, o Traço de Pillai deve ser priorizado como o critério\n  mais robusto a desvios de homocedasticidade multivariada.\n\n",
                  .fmt_f(res_box$M, 2), .fmt_f(res_box$F, 2), .fmt_p(res_box$p.value)))
    } else {
      cat(sprintf("  O teste M de Box indicou homogeneidade das matrizes de covariância entre os grupos\n  (M = %s, F = %s, p = %s), atendendo plenamente ao pressuposto clássico da MANOVA.\n\n",
                  .fmt_f(res_box$M, 2), .fmt_f(res_box$F, 2), .fmt_p(res_box$p.value)))
    }
  }

  sig_aov <- tab_aov$Variavel[tab_aov$p_bonf < 0.05]
  if (length(sig_aov) == length(nomes_vars)) {
    cat("  Nas análises univariadas de acompanhamento com correção de Bonferroni, todas as variáveis resposta\n  apresentaram diferenças estatisticamente significativas entre os grupos (p < 0,05).\n")
  } else if (length(sig_aov) > 0) {
    cat(sprintf("  Nas análises univariadas de acompanhamento com correção de Bonferroni, as variáveis que\n  apresentaram diferenças estatisticamente significativas entre os grupos foram:\n  %s (p < 0,05).\n",
                paste(paste0("'", sig_aov, "'"), collapse = ", ")))
  } else {
    cat("  Nas análises univariadas de acompanhamento com correção de Bonferroni, nenhuma variável resposta\n  individual manteve significância estatística após o ajuste de comparações múltiplas.\n")
  }

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(paste0(strrep("─", 78), "\n\n"))
  }

  resultado <- list(
    modelo_manova       = mod_manova,
    box_m               = res_box,
    tabela_multivariada = tab_multi,
    anovas_univariadas  = tab_aov,
    medias_grupos       = medias_df,
    variaveis           = nomes_vars,
    grupo               = nome_grupo,
    dados               = df_analise
  )
  class(resultado) <- c("estatR_manova", "list")

  if (isTRUE(grafico)) {
    print(plot_manova(resultado))
  }

  invisible(resultado)
}

#' @title Gráfico de Perfis Multivariados da MANOVA
#' @description Plota os perfis médios padronizados (escores Z) de cada grupo
#' ao longo das variáveis resposta analisadas na MANOVA.
#' @param modelo Objeto retornado pela função \code{analise_manova}.
#' @return Objeto ggplot2.
#' @export
plot_manova <- function(modelo) {
  if (!inherits(modelo, "estatR_manova")) {
    stop("O objeto fornecido deve ser retornado pela função 'analise_manova'.")
  }

  df_dados   <- modelo$dados
  nomes_vars <- modelo$variaveis
  nome_grupo <- modelo$grupo

  # Padronizar variáveis em escore Z para comparar na mesma escala
  df_pad <- df_dados
  for (v in nomes_vars) {
    s_val <- stats::sd(df_pad[[v]], na.rm = TRUE)
    m_val <- mean(df_pad[[v]], na.rm = TRUE)
    df_pad[[v]] <- if (s_val > 0) (df_pad[[v]] - m_val) / s_val else 0
  }

  # Calcular médias de Z por grupo
  niveis_grp <- levels(as.factor(df_dados[[nome_grupo]]))
  perfis_lista <- lapply(niveis_grp, function(g) {
    sub_df <- df_pad[df_pad[[nome_grupo]] == g, , drop = FALSE]
    medias <- sapply(nomes_vars, function(v) mean(sub_df[[v]], na.rm = TRUE))
    data.frame(
      Variavel = factor(nomes_vars, levels = nomes_vars),
      Grupo    = as.character(g),
      Media_Z  = medias,
      stringsAsFactors = FALSE
    )
  })
  df_plot <- do.call(rbind, perfis_lista)
  df_plot$Grupo <- factor(df_plot$Grupo, levels = niveis_grp)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Variavel, y = Media_Z, color = Grupo, group = Grupo)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title    = "Perfis Médios Multivariados (Escores Z)",
      subtitle = sprintf("Comparação dos grupos do fator: %s", nome_grupo),
      x        = "Variáveis Resposta",
      y        = "Média Padronizada (Z)",
      color    = nome_grupo
    )

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}

# ==============================================================================
# MÓDULO 5: ANÁLISE DISCRIMINANTE (LDA / QDA)
# ==============================================================================

# Helper interno: cálculo de métricas de classificação
.calc_metricas_classif <- function(tab) {
  N   <- sum(tab)
  acc <- sum(diag(tab)) / N
  err <- 1 - acc
  
  # Kappa de Cohen
  p_o   <- acc
  p_e   <- sum(rowSums(tab) * colSums(tab)) / (N^2)
  kappa <- if (p_e < 1) (p_o - p_e) / (1 - p_e) else 1
  
  # Macro Sensibilidade, Precisão e F1
  k        <- nrow(tab)
  sens_vec <- numeric(k)
  prec_vec <- numeric(k)
  f1_vec   <- numeric(k)
  
  for (i in seq_len(k)) {
    tp <- tab[i, i]
    fn <- sum(tab[i, -i])
    fp <- sum(tab[-i, i])
    
    sens_vec[i] <- if (tp + fn > 0) tp / (tp + fn) else NA
    prec_vec[i] <- if (tp + fp > 0) tp / (tp + fp) else NA
    
    if (!is.na(sens_vec[i]) && !is.na(prec_vec[i]) && (sens_vec[i] + prec_vec[i] > 0)) {
      f1_vec[i] <- 2 * (sens_vec[i] * prec_vec[i]) / (sens_vec[i] + prec_vec[i])
    } else {
      f1_vec[i] <- 0
    }
  }
  
  list(
    acuracia      = acc,
    taxa_erro     = err,
    kappa         = kappa,
    sensibilidade = mean(sens_vec, na.rm = TRUE),
    precisao      = mean(prec_vec, na.rm = TRUE),
    f1_score      = mean(f1_vec, na.rm = TRUE)
  )
}

#' @title Análise Discriminante (LDA / QDA)
#' @description Realiza a Análise Discriminante Linear de Fisher (LDA) ou Quadrática (QDA),
#' avaliando a capacidade de preditores numéricos contínuos em separar e classificar categorias.
#' Apresenta probabilidades a priori, autovalores e correlações canônicas, coeficientes das
#' funções discriminantes, matriz de confusão, métricas completas de classificação
#' (Acurácia, Taxa de Erro, Kappa, Macro Sensibilidade, Precisão, F1-Score com validação cruzada LOO),
#' teste de homogeneidade de covariâncias (M de Box) e interpretação automatizada.
#'
#' @param dados Data frame contendo os dados.
#' @param variaveis Vetor de variáveis preditoras contínuas (com ou sem aspas, ex: \code{c(x1, x2)} ou \code{c("x1", "x2")}).
#'   Se omitido, seleciona todas as variáveis numéricas do data frame (exceto a variável de grupo).
#' @param grupo Variável categórica de agrupamento/resposta (com ou sem aspas, ex: \code{classe} ou \code{"classe"}).
#' @param metodo Método discriminante: \code{"lda"} (padrão, Linear de Fisher) ou \code{"qda"} (Quadrática).
#' @param validacao_cruzada Lógico. Se TRUE, realiza validação cruzada leave-one-out (LOO). Padrão TRUE.
#' @param grafico Lógico. Se TRUE, plota o gráfico discriminante (biplot 2D ou histograma). Padrão FALSE.
#'
#' @return Retorna um objeto de classe \code{"estatR_lda"} contendo:
#' \itemize{
#'   \item \code{modelo_mass}: Objeto ajustado da função \code{MASS::lda} ou \code{MASS::qda}
#'   \item \code{metodo}: Método utilizado ("lda" ou "qda")
#'   \item \code{autovalores}: Tabela de autovalores, variância explicada e correlação canônica
#'   \item \code{coeficientes}: Coeficientes lineares das funções discriminantes (scaling)
#'   \item \code{escores}: Escores discriminantes projetados para cada observação
#'   \item \code{confusao_resub}: Matriz de confusão por resubstituição
#'   \item \code{confusao_cv}: Matriz de confusão por validação cruzada LOO
#'   \item \code{metricas}: Tabela comparativa de métricas de classificação
#'   \item \code{box_m}: Resultados do teste M de Box
#'   \item \code{medias_grupo}: Médias das variáveis preditoras por grupo
#'   \item \code{priori}: Probabilidades a priori dos grupos
#'   \item \code{variaveis}: Nomes das variáveis preditoras
#'   \item \code{grupo}: Nome da variável de agrupamento
#'   \item \code{dados}: Data frame limpo utilizado
#' }
#' @export
analise_discriminante <- function(dados = NULL, variaveis = NULL, grupo = NULL,
                                  metodo = c("lda", "qda"),
                                  validacao_cruzada = TRUE,
                                  grafico = TRUE) {
  metodo <- match.arg(metodo)

  # 1. Tratamento flexível de argumentos posicionais e não avaliados
  sub_dados <- substitute(dados)
  sub_vars  <- substitute(variaveis)
  sub_grupo <- substitute(grupo)

  # Se o usuário passou: analise_discriminante(variaveis = ..., grupo = ..., dados = ...)
  # ou analise_discriminante(dados, variaveis, grupo)
  if (is.null(dados) && !is.null(sub_vars)) {
    val_vars <- tryCatch(variaveis, error = function(e) NULL)
    if (is.data.frame(val_vars)) {
      dados     <- val_vars
      variaveis <- NULL
      sub_vars  <- NULL
    }
  }

  if (is.null(dados) || !is.data.frame(dados)) {
    stop("O argumento 'dados' deve ser fornecido como um data frame.")
  }

  # Extração de 'grupo' (com ou sem aspas)
  nome_grupo <- NULL
  val_grupo  <- tryCatch(grupo, error = function(e) NULL)
  if (is.character(val_grupo) && length(val_grupo) == 1) {
    nome_grupo <- val_grupo
  } else if (!is.null(sub_grupo) && is.symbol(sub_grupo)) {
    nome_grupo <- as.character(sub_grupo)
  }

  if (is.null(nome_grupo) || !(nome_grupo %in% names(dados))) {
    stop("O argumento 'grupo' deve ser uma variável categórica existente no data frame.")
  }

  # Extração de 'variaveis' (com ou sem aspas)
  nomes_vars <- NULL
  val_vars   <- tryCatch(variaveis, error = function(e) NULL)
  if (is.character(val_vars)) {
    nomes_vars <- val_vars
  } else if (!is.null(sub_vars)) {
    if (is.symbol(sub_vars)) {
      nomes_vars <- as.character(sub_vars)
    } else if (is.call(sub_vars)) {
      elementos  <- as.list(sub_vars)[-1]
      nomes_vars <- sapply(elementos, as.character)
    }
  }

  # Se omitido, seleciona todas as numéricas (exceto o grupo)
  if (is.null(nomes_vars) || length(nomes_vars) == 0) {
    todas_num  <- names(dados)[sapply(dados, is.numeric)]
    nomes_vars <- setdiff(todas_num, nome_grupo)
  }

  # Validação de variáveis preditoras
  if (length(nomes_vars) < 2) {
    stop("A Análise Discriminante requer pelo menos 2 variáveis preditoras numéricas contínuas.")
  }
  nao_encontradas <- setdiff(nomes_vars, names(dados))
  if (length(nao_encontradas) > 0) {
    stop(sprintf("As seguintes variáveis não foram encontradas no data frame: %s",
                 paste(nao_encontradas, collapse = ", ")))
  }

  nao_num <- nomes_vars[!sapply(dados[nomes_vars], is.numeric)]
  if (length(nao_num) > 0) {
    stop(sprintf("Todas as variáveis preditoras devem ser numéricas. Variáveis não numéricas: %s",
                 paste(nao_num, collapse = ", ")))
  }

  # 2. Preparação dos Dados
  df_analise <- stats::na.omit(dados[, c(nomes_vars, nome_grupo), drop = FALSE])
  n_total    <- nrow(df_analise)
  fator      <- as.factor(df_analise[[nome_grupo]])
  niveis_grp <- levels(fator)
  k_grupos   <- length(niveis_grp)

  if (k_grupos < 2) {
    stop("A variável de grupo deve possuir pelo menos 2 categorias distintas.")
  }

  X_mat <- as.matrix(df_analise[, nomes_vars, drop = FALSE])

  # 3. Teste M de Box
  res_box <- .box_m_test(X_mat, fator)

  # 4. Ajuste do Modelo Discriminante (MASS)
  if (metodo == "lda") {
    mod_mass <- MASS::lda(X_mat, fator)
    pred_res <- stats::predict(mod_mass, X_mat)
    cls_res  <- pred_res$class
    escores  <- as.data.frame(pred_res$x)
  } else {
    mod_mass <- MASS::qda(X_mat, fator)
    pred_res <- stats::predict(mod_mass, X_mat)
    cls_res  <- pred_res$class
    escores  <- NULL
  }

  if (isTRUE(validacao_cruzada)) {
    mod_cv <- if (metodo == "lda") MASS::lda(X_mat, fator, CV = TRUE) else MASS::qda(X_mat, fator, CV = TRUE)
    cls_cv <- mod_cv$class
  } else {
    cls_cv <- NULL
  }

  # 5. Autovalores, Variância Explicada e Correlação Canônica (específico do LDA)
  tab_autoval <- NULL
  if (metodo == "lda") {
    n_funcs <- min(ncol(X_mat), k_grupos - 1)
    
    # Extração de autovalores via decomposição MANOVA E^-1 H
    fit_manova <- stats::manova(X_mat ~ fator)
    E_mat      <- stats::summary.manova(fit_manova)$SS$Residuals
    H_mat      <- stats::summary.manova(fit_manova)$SS$fator
    eig_raw    <- Re(eigen(solve(E_mat) %*% H_mat)$values)[seq_len(n_funcs)]
    
    # Previne valores negativos/espúrios numéricos
    eig_raw  <- pmax(eig_raw, 0)
    var_exp  <- (eig_raw / sum(eig_raw)) * 100
    var_acum <- cumsum(var_exp)
    r_canon  <- sqrt(eig_raw / (1 + eig_raw))
    
    tab_autoval <- data.frame(
      Funcao        = paste0("LD", seq_len(n_funcs)),
      Autovalor     = eig_raw,
      Var_Expl_Pct  = var_exp,
      Var_Acum_Pct  = var_acum,
      Corr_Canonica = r_canon,
      stringsAsFactors = FALSE
    )
  }

  # 6. Matrizes de Confusão e Métricas
  conf_resub <- table(Observado = fator, Previsto = cls_res)
  met_resub  <- .calc_metricas_classif(conf_resub)

  if (isTRUE(validacao_cruzada)) {
    conf_cv <- table(Observado = fator, Previsto = cls_cv)
    met_cv  <- .calc_metricas_classif(conf_cv)
  } else {
    conf_cv <- NULL
    met_cv  <- NULL
  }

  # Tabela comparativa de métricas
  df_metricas <- data.frame(
    Metrica = c("Acurácia Global", "Taxa de Erro", "Kappa de Cohen",
                "Sensibilidade Média (Recall)", "Precisão Média", "F1-Score Médio"),
    Resubstituicao = c(met_resub$acuracia, met_resub$taxa_erro, met_resub$kappa,
                       met_resub$sensibilidade, met_resub$precisao, met_resub$f1_score),
    Validacao_Cruzada = if (!is.null(met_cv)) c(met_cv$acuracia, met_cv$taxa_erro, met_cv$kappa,
                                                met_cv$sensibilidade, met_cv$precisao, met_cv$f1_score) else rep(NA_real_, 6),
    stringsAsFactors = FALSE
  )

  # 7. Tabela de Probabilidades a Priori
  df_priori <- data.frame(
    Grupo    = names(mod_mass$prior),
    N        = as.numeric(table(fator)),
    A_Priori = as.numeric(mod_mass$prior),
    stringsAsFactors = FALSE
  )

  # Helpers de formatação
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 2) {
    if (is.na(x)) return("-")
    formatC(x, format = "f", digits = d, decimal.mark = ",")
  }
  .fmt_p <- function(x) {
    if (is.na(x)) return("-")
    if (x < 0.001) return("< 0,001")
    formatC(x, format = "f", digits = 3, decimal.mark = ",")
  }
  .fmt_pct <- function(x, d = 1) {
    if (is.na(x)) return("-")
    paste0(formatC(x * 100, format = "f", digits = d, decimal.mark = ","), "%")
  }
  .ast <- function(p) {
    if (is.na(p) || p >= 0.10) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    return("*")
  }

  # ==========================================================================
  # IMPRESSÃO NO CONSOLE
  # ==========================================================================
  titulo_met <- if (metodo == "lda") "ANÁLISE DISCRIMINANTE LINEAR (LDA)" else "ANÁLISE DISCRIMINANTE QUADRÁTICA (QDA)"
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo(titulo_met)
  } else {
    cat(sprintf("\n── %s ──\n", titulo_met))
  }

  cat(sprintf("  Fator / Grupo:        %s (%d grupos: %s)\n", nome_grupo, k_grupos, paste(niveis_grp, collapse = ", ")))
  cat(sprintf("  Variáveis Preditoras: %s (%d variáveis)\n", paste(nomes_vars, collapse = ", "), length(nomes_vars)))
  cat(sprintf("  Observações Válidas:  %d\n", n_total))
  cat(sprintf("  Método:               %s\n\n", if (metodo == "lda") "LDA (Linear de Fisher)" else "QDA (Quadrática)"))

  # SEÇÃO: PROBABILIDADES A PRIORI E GRUPOS
  if (exists(".print_topico", mode = "function")) {
    .print_topico("PROBABILIDADES A PRIORI E GRUPOS")
  } else {
    cat("  ▶ PROBABILIDADES A PRIORI E GRUPOS\n")
  }

  w_priori <- c(grp = max(18, max(nchar(niveis_grp)) + 2), n = 12, pr = 18)
  sep_priori <- paste0("  ", strrep("─", sum(w_priori)))
  cat(sep_priori, "\n")
  cat("  ", .pad("Grupo", w_priori["grp"], "left"),
             .pad("N", w_priori["n"], "center"),
             .pad("Prob. a Priori", w_priori["pr"], "center"), "\n", sep = "")
  cat(sep_priori, "\n")
  for (r in seq_len(nrow(df_priori))) {
    cat("  ", .pad(df_priori$Grupo[r], w_priori["grp"], "left"),
               .pad(as.character(df_priori$N[r]), w_priori["n"], "center"),
               .pad(.fmt_f(df_priori$A_Priori[r], 4), w_priori["pr"], "center"), "\n", sep = "")
  }
  cat(sep_priori, "\n\n")

  # SEÇÃO: FUNÇÕES DISCRIMINANTES E AUTOVALORES (somente LDA)
  if (metodo == "lda" && !is.null(tab_autoval)) {
    if (exists(".print_topico", mode = "function")) {
      .print_topico("FUNÇÕES DISCRIMINANTES E AUTOVALORES")
    } else {
      cat("  ▶ FUNÇÕES DISCRIMINANTES E AUTOVALORES\n")
    }

    w_eig <- c(fn = 12, eig = 14, ve = 16, va = 16, rc = 16)
    sep_eig <- paste0("  ", strrep("─", sum(w_eig)))
    cat(sep_eig, "\n")
    cat("  ", .pad("Função", w_eig["fn"], "left"),
               .pad("Autovalor", w_eig["eig"], "center"),
               .pad("Var. Expl. (%)", w_eig["ve"], "center"),
               .pad("Var. Acum. (%)", w_eig["va"], "center"),
               .pad("Corr. Canônica", w_eig["rc"], "center"), "\n", sep = "")
    cat(sep_eig, "\n")
    for (r in seq_len(nrow(tab_autoval))) {
      cat("  ", .pad(tab_autoval$Funcao[r], w_eig["fn"], "left"),
                 .pad(.fmt_f(tab_autoval$Autovalor[r], 3), w_eig["eig"], "center"),
                 .pad(paste0(.fmt_f(tab_autoval$Var_Expl_Pct[r], 1), "%"), w_eig["ve"], "center"),
                 .pad(paste0(.fmt_f(tab_autoval$Var_Acum_Pct[r], 1), "%"), w_eig["va"], "center"),
                 .pad(.fmt_f(tab_autoval$Corr_Canonica[r], 3), w_eig["rc"], "center"), "\n", sep = "")
    }
    cat(sep_eig, "\n\n")

    # SEÇÃO: COEFICIENTES DAS FUNÇÕES DISCRIMINANTES (SCALING)
    if (exists(".print_topico", mode = "function")) {
      .print_topico("COEFICIENTES DAS FUNÇÕES DISCRIMINANTES (SCALING)")
    } else {
      cat("  ▶ COEFICIENTES DAS FUNÇÕES DISCRIMINANTES (SCALING)\n")
    }

    coefs_lda <- as.data.frame(mod_mass$scaling)
    w_var_coef <- max(20, max(nchar(nomes_vars)) + 2)
    w_ld_col   <- 14
    total_coef_w <- w_var_coef + w_ld_col * ncol(coefs_lda)
    sep_coef <- paste0("  ", strrep("─", total_coef_w))
    cat(sep_coef, "\n")
    cat("  ", .pad("Variável Preditora", w_var_coef, "left"))
    for (nm_ld in colnames(coefs_lda)) {
      cat(.pad(nm_ld, w_ld_col, "center"))
    }
    cat("\n")
    cat(sep_coef, "\n")
    for (r in seq_len(nrow(coefs_lda))) {
      cat("  ", .pad(rownames(coefs_lda)[r], w_var_coef, "left"))
      for (c_idx in seq_len(ncol(coefs_lda))) {
        cat(.pad(.fmt_f(coefs_lda[r, c_idx], 3), w_ld_col, "center"))
      }
      cat("\n")
    }
    cat(sep_coef, "\n\n")
  }

  # SEÇÃO: MATRIZ DE CONFUSÃO (RESUBSTITUIÇÃO)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("MATRIZ DE CONFUSÃO (RESUBSTITUIÇÃO)")
  } else {
    cat("  ▶ MATRIZ DE CONFUSÃO (RESUBSTITUIÇÃO)\n")
  }

  w_rot_conf <- max(22, max(nchar(niveis_grp)) + 2)
  w_col_conf <- max(12, max(nchar(niveis_grp)) + 2)
  total_conf_w <- w_rot_conf + w_col_conf * k_grupos
  sep_conf <- paste0("  ", strrep("─", total_conf_w))
  cat(sep_conf, "\n")
  cat("  ", .pad("Observado \\ Previsto", w_rot_conf, "left"))
  for (g in niveis_grp) {
    cat(.pad(g, w_col_conf, "center"))
  }
  cat("\n")
  cat(sep_conf, "\n")
  for (r in seq_len(nrow(conf_resub))) {
    cat("  ", .pad(rownames(conf_resub)[r], w_rot_conf, "left"))
    for (c_idx in seq_len(ncol(conf_resub))) {
      cat(.pad(as.character(conf_resub[r, c_idx]), w_col_conf, "center"))
    }
    cat("\n")
  }
  cat(sep_conf, "\n\n")

  # SEÇÃO: MÉTRICAS DE CLASSIFICAÇÃO E DESEMPENHO
  if (exists(".print_topico", mode = "function")) {
    .print_topico("MÉTRICAS DE CLASSIFICAÇÃO E DESEMPENHO")
  } else {
    cat("  ▶ MÉTRICAS DE CLASSIFICAÇÃO E DESEMPENHO\n")
  }

  w_met <- c(m = 32, resub = 18, cv = 24)
  sep_met <- paste0("  ", strrep("─", sum(w_met)))
  cat(sep_met, "\n")
  cat("  ", .pad("Métrica", w_met["m"], "left"),
             .pad("Resubstituição", w_met["resub"], "center"),
             .pad("Validação Cruzada (CV)", w_met["cv"], "center"), "\n", sep = "")
  cat(sep_met, "\n")
  for (r in seq_len(nrow(df_metricas))) {
    v_resub <- if (r %in% c(1, 2)) .fmt_pct(df_metricas$Resubstituicao[r], 1) else .fmt_f(df_metricas$Resubstituicao[r], 3)
    v_cv    <- if (is.na(df_metricas$Validacao_Cruzada[r])) "-" else if (r %in% c(1, 2)) .fmt_pct(df_metricas$Validacao_Cruzada[r], 1) else .fmt_f(df_metricas$Validacao_Cruzada[r], 3)
    cat("  ", .pad(df_metricas$Metrica[r], w_met["m"], "left"),
               .pad(v_resub, w_met["resub"], "center"),
               .pad(v_cv, w_met["cv"], "center"), "\n", sep = "")
  }
  cat(sep_met, "\n\n")

  # SEÇÃO: TESTE DE HOMOGENEIDADE DAS COVARIÂNCIAS (PRESSUPOSTO)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("TESTE DE HOMOGENEIDADE DAS COVARIÂNCIAS (PRESSUPOSTO)")
  } else {
    cat("  ▶ TESTE DE HOMOGENEIDADE DAS COVARIÂNCIAS (PRESSUPOSTO)\n")
  }

  w_box <- c(tst = 22, m = 18, f = 14, p = 14)
  sep_box <- paste0("  ", strrep("─", sum(w_box)))
  cat(sep_box, "\n")
  cat("  ", .pad("Teste", w_box["tst"], "left"),
             .pad("Estatística (M)", w_box["m"], "center"),
             .pad("F aprox.", w_box["f"], "center"),
             .pad("p-valor", w_box["p"], "center"), "\n", sep = "")
  cat(sep_box, "\n")
  box_p_ast <- paste0(.fmt_p(res_box$p.value), " ", .ast(res_box$p.value))
  cat("  ", .pad("M de Box", w_box["tst"], "left"),
             .pad(.fmt_f(res_box$M, 2), w_box["m"], "center"),
             .pad(.fmt_f(res_box$F, 2), w_box["f"], "center"),
             .pad(box_p_ast, w_box["p"], "center"), "\n", sep = "")
  cat(sep_box, "\n")
  cat("  * H0: As matrizes de covariância são homogêneas (pressuposto do LDA).\n")
  cat("  * Nota: Se M de Box for significativo, QDA pode ser uma alternativa recomendada.\n\n")

  # SEÇÃO: INTERPRETAÇÃO
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO")
  } else {
    cat("  ▶ INTERPRETAÇÃO\n")
  }

  # 1º Parágrafo: Acurácia e desempenho
  cv_txt <- if (!is.null(met_cv)) sprintf(" e %s sob validação cruzada leave-one-out", .fmt_pct(met_cv$acuracia, 1)) else ""
  overfit_txt <- if (!is.null(met_cv) && (met_resub$acuracia - met_cv$acuracia > 0.05)) {
    "A queda de desempenho observada na validação cruzada sugere atenção para um possível sobreajuste (overfitting)."
  } else {
    "indicando excelente capacidade preditiva e consistência na generalização do modelo."
  }
  cat(sprintf("  A Análise Discriminante (%s) obteve uma Acurácia Global de %s por resubstituição%s,\n  com coeficiente Kappa de %s, %s\n\n",
              toupper(metodo), .fmt_pct(met_resub$acuracia, 1), cv_txt, .fmt_f(met_resub$kappa, 3), overfit_txt))

  # 2º Parágrafo: Funções discriminantes (se LDA)
  if (metodo == "lda" && !is.null(tab_autoval)) {
    ld1_pct <- tab_autoval$Var_Expl_Pct[1]
    ld1_rc  <- tab_autoval$Corr_Canonica[1]
    cat(sprintf("  A primeira função discriminante (LD1) concentrou %s da variabilidade entre os grupos com correlação\n  canônica de %s, sendo a dimensão predominante para discriminar as categorias de '%s'.\n\n",
                paste0(.fmt_f(ld1_pct, 1), "%"), .fmt_f(ld1_rc, 3), nome_grupo))

    # Variável de maior peso em LD1
    coefs_ld1 <- mod_mass$scaling[, 1]
    var_maior <- names(which.max(abs(coefs_ld1)))
    peso_maior <- coefs_ld1[var_maior]
    cat(sprintf("  Dentre os preditores analisados, a variável '%s' apresentou o maior peso absoluto na primeira função\n  discriminante (coeficiente = %s), evidenciando-se como o indicador mais determinante na separação dos grupos.\n\n",
                var_maior, .fmt_f(peso_maior, 3)))
  }

  # 3º Parágrafo: Pressuposto do M de Box
  if (!is.na(res_box$p.value)) {
    if (res_box$p.value < 0.05) {
      cat(sprintf("  O teste M de Box indicou rejeição da homogeneidade das matrizes de covariância (M = %s, p %s).\n  Embora o LDA seja razoavelmente robusto com grupos balanceados, a Análise Discriminante Quadrática (QDA)\n  ou procedimentos de regularização devem ser considerados caso haja interesse em fronteiras não-lineares.\n",
                  .fmt_f(res_box$M, 2), .fmt_p(res_box$p.value)))
    } else {
      cat(sprintf("  O teste M de Box indicou homogeneidade das matrizes de covariância entre os grupos (M = %s, p = %s),\n  atendendo com sucesso ao pressuposto fundamental da Análise Discriminante Linear.\n",
                  .fmt_f(res_box$M, 2), .fmt_p(res_box$p.value)))
    }
  }

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat(paste0(strrep("─", 78), "\n\n"))
  }

  # Montagem do objeto de retorno
  if (metodo == "lda") {
    escores$Grupo <- fator
  }

  resultado <- list(
    modelo_mass    = mod_mass,
    metodo         = metodo,
    autovalores    = tab_autoval,
    coeficientes   = if (metodo == "lda") mod_mass$scaling else NULL,
    escores        = escores,
    confusao_resub = conf_resub,
    confusao_cv    = conf_cv,
    metricas       = df_metricas,
    box_m          = res_box,
    medias_grupo   = mod_mass$means,
    priori         = mod_mass$prior,
    variaveis      = nomes_vars,
    grupo          = nome_grupo,
    dados          = df_analise
  )
  class(resultado) <- c("estatR_lda", "list")

  if (isTRUE(grafico)) {
    print(biplot_lda(resultado))
  }

  invisible(resultado)
}

#' @title Gráfico da Análise Discriminante (Biplot / Dispersão)
#' @description Plota a projeção das observações no espaço discriminante.
#' Para 3 ou mais grupos (2+ funções discriminantes), gera um Biplot 2D (LD1 x LD2) com elipses
#' de confiança de 95% e setas das variáveis originais. Para 2 grupos (1 função),
#' plota a distribuição/histograma dos escores de LD1 separados por cor.
#' @param modelo Objeto retornado pela função \code{analise_discriminante}.
#' @param elipses Lógico. Se TRUE, exibe elipses de confiança de 95% para cada grupo. Padrão TRUE.
#' @return Objeto ggplot2.
#' @export
biplot_lda <- function(modelo, elipses = TRUE) {
  if (!inherits(modelo, "estatR_lda")) {
    stop("O objeto fornecido deve ser retornado pela função 'analise_discriminante'.")
  }

  if (modelo$metodo != "lda") {
    stop("O biplot discriminante está disponível apenas para modelos lineares (metodo = 'lda').")
  }

  escores    <- modelo$escores
  nome_grupo <- modelo$grupo
  n_funcs    <- if (!is.null(modelo$autovalores)) nrow(modelo$autovalores) else (ncol(escores) - 1)

  # Caso 1: 2 ou mais funções discriminantes -> Biplot 2D (LD1 x LD2)
  if (n_funcs >= 2) {
    loads_mat <- modelo$coeficientes[, 1:2, drop = FALSE]
    loads_df  <- as.data.frame(loads_mat)
    names(loads_df) <- c("LD1", "LD2")
    loads_df$Variavel <- rownames(loads_df)

    # Escala para projetar setas esteticamente com os pontos
    escala_setas <- (max(abs(escores[, c("LD1", "LD2")])) /
                     max(abs(loads_df[, c("LD1", "LD2")]))) * 0.70

    pct_ld1 <- modelo$autovalores$Var_Expl_Pct[1]
    pct_ld2 <- modelo$autovalores$Var_Expl_Pct[2]
    rot_x   <- sprintf("LD1 (%.1f%%)", pct_ld1)
    rot_y   <- sprintf("LD2 (%.1f%%)", pct_ld2)

    p <- ggplot2::ggplot()

    if (isTRUE(elipses)) {
      p <- p + ggplot2::stat_ellipse(
        data = escores,
        ggplot2::aes(x = LD1, y = LD2, color = Grupo, fill = Grupo),
        geom = "polygon", alpha = 0.08, linewidth = 0.6, level = 0.95
      )
    }

    p <- p +
      ggplot2::geom_point(
        data = escores,
        ggplot2::aes(x = LD1, y = LD2, color = Grupo),
        size = 2.5, alpha = 0.80
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
      ggplot2::geom_segment(
        data = loads_df,
        ggplot2::aes(x = 0, y = 0, xend = LD1 * escala_setas, yend = LD2 * escala_setas),
        color = "gray20", linewidth = 0.8,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.22, "cm"), type = "closed")
      )

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = loads_df,
        ggplot2::aes(x = LD1 * escala_setas * 1.1, y = LD2 * escala_setas * 1.1, label = Variavel),
        color = "gray10", size = 3.5, fontface = "bold", box.padding = 0.35
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = loads_df,
        ggplot2::aes(x = LD1 * escala_setas * 1.15, y = LD2 * escala_setas * 1.15, label = Variavel),
        color = "gray10", size = 3.5, fontface = "bold"
      )
    }

    p <- p + ggplot2::labs(
      title    = "Biplot — Análise Discriminante Linear (LDA)",
      subtitle = sprintf("%s vs %s | Elipses de confiança 95%% por grupo", rot_x, rot_y),
      x        = rot_x,
      y        = rot_y,
      color    = nome_grupo,
      fill     = nome_grupo
    )

  } else {
    # Caso 2: Apenas 1 função discriminante (2 grupos) -> Histograma de distribuição de LD1
    p <- ggplot2::ggplot(escores, ggplot2::aes(x = LD1, fill = Grupo, color = Grupo)) +
      ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 20) +
      ggplot2::labs(
        title    = "Distribuição do Escore Discriminante (LD1)",
        subtitle = sprintf("Separação entre os grupos do fator: %s", nome_grupo),
        x        = "Escore Discriminante (LD1)",
        y        = "Frequência",
        fill     = nome_grupo,
        color    = nome_grupo
      )
  }

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}

