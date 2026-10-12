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
                grafico          = FALSE) {

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
        .pad(paste0(.fmt_f(var_expl_pct[i], 2), "%"), w_var["expl"], "center"),
        .pad(paste0(.fmt_f(var_acum_pct[i], 2), "%"), w_var["acum"], "center"),
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
          .pad(paste0(.fmt_f(var_acum_pct[i], 2), "%"), w_ap["acum"], "center"),
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

  p <- ggplot2::ggplot(df_scree, ggplot2::aes(x = Componente_num, y = Autovalor)) +
    ggplot2::geom_line(color = "#1D3557", linewidth = 1) +
    ggplot2::geom_point(color = "#1D3557", size = 3) +
    ggplot2::geom_hline(yintercept = 1.0, linetype = "dashed", color = "#D90429", linewidth = 0.8) +
    ggplot2::scale_x_continuous(breaks = df_scree$Componente_num, labels = df_scree$Componente) +
    ggplot2::labs(
      title    = "Scree Plot — ACP",
      subtitle = "Autovalores Observados vs Critério de Kaiser (linha vermelha: autovalor = 1)",
      x = "Componente Principal",
      y = "Autovalor"
    )

  if (analise_paralela && !is.null(modelo$criterios$autovalores_simulados_p95) &&
      !any(is.na(modelo$criterios$autovalores_simulados_p95))) {
    df_ap <- data.frame(
      Componente_num = seq_along(modelo$criterios$autovalores_simulados_p95),
      Autovalor_AP   = modelo$criterios$autovalores_simulados_p95
    )
    p <- p +
      ggplot2::geom_line(data = df_ap, ggplot2::aes(x = Componente_num, y = Autovalor_AP),
                         color = "#457B9D", linetype = "dotted", linewidth = 1) +
      ggplot2::geom_point(data = df_ap, ggplot2::aes(x = Componente_num, y = Autovalor_AP),
                          color = "#457B9D", size = 2.5, shape = 17) +
      ggplot2::labs(subtitle = "Linha azul contínua: Observado | Linha vermelha: Kaiser (1,0) | Linha pontilhada: Horn (P95 Simulado)")
  }

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2)
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}
