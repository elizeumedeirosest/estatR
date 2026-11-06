# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: MODELAGEM
# Funções de preparação de dados para modelagem e machine learning
# ─────────────────────────────────────────────────────────────────────────────

# ── Helper interno ────────────────────────────────────────────────────────────
.mod_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
  paste0(s, strrep(" ", sp))
}

.mod_sep <- function(w) paste0("  ", strrep("\u2500", w))


# ─────────────────────────────────────────────────────────────────────────────
#' @title Divisão Treino / Teste
#' @description Divide um conjunto de dados em partições de treino e teste para
#'   validação de modelos. Suporta três modos:
#'   \itemize{
#'     \item \strong{Aleatório} (padrão): sorteio simples sem reposição.
#'     \item \strong{Estratificado}: o sorteio é feito dentro de cada nível
#'       do estrato, mantendo as proporções dos grupos em ambas as partições.
#'       Ideal quando há classes desbalanceadas.
#'     \item \strong{Temporal}: corte cronológico sem embaralhamento. Os primeiros
#'       \code{proporcao * 100\%} das linhas vão para o treino e o restante
#'       para o teste. Essencial para séries temporais (evita data leakage).
#'   }
#' @param dados Data frame com os dados a dividir.
#' @param proporcao Proporção das observações destinadas ao treino (padrão \code{0.8}).
#' @param estrato Nome da coluna de estratificação (sem aspas), opcional.
#'   Quando fornecido, o sorteio é feito dentro de cada grupo do estrato.
#' @param temporal Lógico. Se TRUE, realiza corte cronológico (sem embaralhamento).
#'   Não compatível com \code{estrato}.
#' @param semente Valor inteiro para reprodutibilidade do sorteio (padrão \code{NULL}).
#' @return Lista com dois data frames: \code{$treino} e \code{$teste}.
#' @export
dividir_dados <- function(dados, proporcao = 0.8, estrato = NULL,
                           temporal = FALSE, semente = NULL) {

  if (!is.data.frame(dados)) stop("'dados' deve ser um data frame.")
  if (proporcao <= 0 || proporcao >= 1) stop("'proporcao' deve estar entre 0 e 1 (exclusive).")

  estrato_sub  <- substitute(estrato)
  estrato_nome <- if (!is.null(estrato_sub) && deparse(estrato_sub) != "NULL")
    sub(".*\\$", "", deparse(estrato_sub)) else NULL

  if (temporal && !is.null(estrato_nome)) {
    stop("'temporal' e 'estrato' n\u00e3o podem ser usados simultaneamente.")
  }

  n <- nrow(dados)

  # ── Modo temporal ─────────────────────────────────────────────────────────
  if (temporal) {
    idx_treino <- 1:floor(n * proporcao)
    idx_teste  <- (floor(n * proporcao) + 1):n

    treino <- dados[idx_treino, , drop = FALSE]
    teste  <- dados[idx_teste,  , drop = FALSE]

    w_sep <- 67
  .print_titulo("DIVIS\u00c3O TREINO / TESTE (TEMPORAL) %s")
    cat(sprintf("  Modo: Corte cronol\u00f3gico (sem embaralhamento)\n"))
    cat(sprintf("  Propor\u00e7\u00e3o: %.0f%% treino / %.0f%% teste\n\n",
                proporcao * 100, (1 - proporcao) * 100))

    wc <- c(part = 12, n = 12, pct = 14)
    hdr <- paste0(.mod_pad("Parti\u00e7\u00e3o", wc["part"], "left"),
                  .mod_pad("N", wc["n"], "center"),
                  .mod_pad("% Real", wc["pct"], "center"))
    cat(.mod_sep(nchar(hdr)), "\n")
    cat("  ", hdr, "\n", sep = "")
    cat(.mod_sep(nchar(hdr)), "\n")
    cat("  ", .mod_pad("Treino", wc["part"], "left"),
        .mod_pad(nrow(treino), wc["n"], "center"),
        .mod_pad(sprintf("%.1f%%", nrow(treino)/n*100), wc["pct"], "center"), "\n", sep = "")
    cat("  ", .mod_pad("Teste", wc["part"], "left"),
        .mod_pad(nrow(teste), wc["n"], "center"),
        .mod_pad(sprintf("%.1f%%", nrow(teste)/n*100), wc["pct"], "center"), "\n", sep = "")
    cat("  ", .mod_pad("Total", wc["part"], "left"),
        .mod_pad(n, wc["n"], "center"),
        .mod_pad("100,0%", wc["pct"], "center"), "\n", sep = "")
    cat(.mod_sep(nchar(hdr)), "\n")
    cat("  Nota: as linhas foram mantidas na ordem original (corte cronol\u00f3gico).\n")
  .print_rodape()

    return(invisible(list(treino = treino, teste = teste)))
  }

  if (!is.null(semente)) set.seed(semente)

  # ── Modo estratificado ────────────────────────────────────────────────────
  if (!is.null(estrato_nome)) {
    if (!(estrato_nome %in% names(dados)))
      stop(sprintf("Coluna '%s' n\u00e3o encontrada.", estrato_nome))

    grupos <- as.factor(dados[[estrato_nome]])
    niveis <- levels(grupos)

    idx_treino <- c()
    resumo_list <- list()

    for (niv in niveis) {
      idx_grupo <- which(grupos == niv)
      n_grupo   <- length(idx_grupo)
      n_treino  <- max(1, round(n_grupo * proporcao))
      sel       <- sample(idx_grupo, n_treino)
      idx_treino <- c(idx_treino, sel)

      resumo_list[[niv]] <- data.frame(
        estrato  = niv,
        n_total  = n_grupo,
        n_treino = n_treino,
        n_teste  = n_grupo - n_treino
      )
    }

    idx_teste <- setdiff(1:n, idx_treino)
    treino    <- dados[sort(idx_treino), , drop = FALSE]
    teste     <- dados[sort(idx_teste),  , drop = FALSE]

    resumo <- do.call(rbind, resumo_list)

    w_sep <- 67
  .print_titulo("DIVIS\u00c3O TREINO / TESTE (ESTRATIFICADA) %s")
    cat(sprintf("  Modo: Estratificado por '%s'\n", estrato_nome))
    cat(sprintf("  Propor\u00e7\u00e3o: %.0f%% treino / %.0f%% teste\n\n",
                proporcao * 100, (1 - proporcao) * 100))

    wc <- c(est = 16, ntot = 12, ntr = 12, nte = 12)
    hdr <- paste0(.mod_pad("Estrato",  wc["est"],  "left"),
                  .mod_pad("N Total",  wc["ntot"], "center"),
                  .mod_pad("N Treino", wc["ntr"],  "center"),
                  .mod_pad("N Teste",  wc["nte"],  "center"))
    cat(.mod_sep(nchar(hdr)), "\n")
    cat("  ", hdr, "\n", sep = "")
    cat(.mod_sep(nchar(hdr)), "\n")
    for (i in 1:nrow(resumo)) {
      cat("  ",
          .mod_pad(resumo$estrato[i],  wc["est"],  "left"),
          .mod_pad(resumo$n_total[i],  wc["ntot"], "center"),
          .mod_pad(resumo$n_treino[i], wc["ntr"],  "center"),
          .mod_pad(resumo$n_teste[i],  wc["nte"],  "center"),
          "\n", sep = "")
    }
    cat(.mod_sep(nchar(hdr)), "\n")
    cat("  ",
        .mod_pad("Total",          wc["est"],  "left"),
        .mod_pad(n,                wc["ntot"], "center"),
        .mod_pad(nrow(treino),     wc["ntr"],  "center"),
        .mod_pad(nrow(teste),      wc["nte"],  "center"),
        "\n", sep = "")
    cat(.mod_sep(nchar(hdr)), "\n")
    cat("  Nota: o sorteio foi realizado dentro de cada estrato.\n")
  .print_rodape()

    return(invisible(list(treino = treino, teste = teste)))
  }

  # ── Modo aleatório (padrão) ────────────────────────────────────────────────
  idx_treino <- sample(1:n, size = floor(n * proporcao))
  idx_teste  <- setdiff(1:n, idx_treino)

  treino <- dados[sort(idx_treino), , drop = FALSE]
  teste  <- dados[sort(idx_teste),  , drop = FALSE]

  w_sep <- 67
  .print_titulo("DIVIS\u00c3O TREINO / TESTE (ALEAT\u00d3RIA) %s")
  cat(sprintf("  Modo: Sorteio aleat\u00f3rio\n"))
  cat(sprintf("  Propor\u00e7\u00e3o: %.0f%% treino / %.0f%% teste\n\n",
              proporcao * 100, (1 - proporcao) * 100))

  wc <- c(part = 12, n = 12, pct = 14)
  hdr <- paste0(.mod_pad("Parti\u00e7\u00e3o", wc["part"], "left"),
                .mod_pad("N", wc["n"], "center"),
                .mod_pad("% Real", wc["pct"], "center"))
  cat(.mod_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.mod_sep(nchar(hdr)), "\n")
  cat("  ", .mod_pad("Treino", wc["part"], "left"),
      .mod_pad(nrow(treino), wc["n"], "center"),
      .mod_pad(sprintf("%.1f%%", nrow(treino)/n*100), wc["pct"], "center"), "\n", sep = "")
  cat("  ", .mod_pad("Teste", wc["part"], "left"),
      .mod_pad(nrow(teste), wc["n"], "center"),
      .mod_pad(sprintf("%.1f%%", nrow(teste)/n*100), wc["pct"], "center"), "\n", sep = "")
  cat("  ", .mod_pad("Total", wc["part"], "left"),
      .mod_pad(n, wc["n"], "center"),
      .mod_pad("100,0%", wc["pct"], "center"), "\n", sep = "")
  cat(.mod_sep(nchar(hdr)), "\n")
  if (!is.null(semente)) cat(sprintf("  Semente utilizada: %s\n", semente))
  .print_rodape()

  invisible(list(treino = treino, teste = teste))
}

# ─────────────────────────────────────────────────────────────────────────────
# COMPARAÇÃO E RANQUEAMENTO UNIFICADO DE MODELOS
# ─────────────────────────────────────────────────────────────────────────────

#' @title Comparacao e Ranqueamento Unificado de Modelos
#' @description Compara dois ou mais modelos estatisticos ajustados,
#'   organizando um ranking pelo criterio escolhido (AICc, AIC, BIC, RMSE ou R2),
#'   incluindo Delta-AIC e Pesos de Akaike.
#'
#' @param ... Dois ou mais modelos ajustados ou uma lista de modelos.
#' @param criterio Criterio de ordenacao: 'aicc' (padrao), 'aic', 'bic', 'rmse' ou 'r2'.
#' @param nomes Vetor opcional com nomes descritivos para os modelos.
#' @return Objeto invisivel da classe estatR_comparacao_modelos.
#' @export
comparar_modelos <- function(..., criterio = c("aicc", "aic", "bic", "rmse", "r2"), nomes = NULL) {
  criterio <- match.arg(criterio)
  args <- list(...)

  # Suporte a passagem de lista unica de modelos
  if (length(args) == 1L && is.list(args[[1]]) && !inherits(args[[1]], c("lm", "glm", "coxph", "Arima", "ets"))) {
    modelos <- args[[1]]
    nomes_arg <- names(modelos)
  } else {
    modelos <- args
    nomes_arg <- as.character(match.call(expand.dots = FALSE)$...)
  }

  if (length(modelos) < 2L) {
    stop("Informe pelo menos 2 modelos para realizar a comparacao.", call. = FALSE)
  }

  if (!is.null(nomes) && length(nomes) == length(modelos)) {
    nomes_modelos <- nomes
  } else if (!is.null(nomes_arg) && all(nzchar(nomes_arg))) {
    nomes_modelos <- nomes_arg
  } else {
    nomes_modelos <- paste0("Modelo_", seq_along(modelos))
  }

  if (inherits(modelos[[1]], "estatR_variograma")) {
    if (exists(".comparar_modelos_variograma", mode = "function")) {
      return(.comparar_modelos_variograma(modelos, nomes = nomes_modelos))
    }
  }
  # Extracao de metricas de cada modelo
  linhas_info <- list()
  for (i in seq_along(modelos)) {
    m <- modelos[[i]]
    if (inherits(m, "estatR_regressao_espacial")) {
      ll   <- m$loglik
      aic  <- m$aic
      bic  <- m$bic
      res  <- m$residuos
      k_p  <- if (!is.null(m$parametros$estimativa)) length(m$parametros$estimativa) else 2L
      n_obs <- if (!is.null(res)) length(res) else 100L
      aicc <- if (!is.na(aic) && !is.na(k_p) && (n_obs - k_p - 1L) > 0) {
        aic + (2 * k_p * (k_p + 1L)) / (n_obs - k_p - 1L)
      } else aic
      rmse <- if (!is.null(res)) sqrt(mean(res^2, na.rm = TRUE)) else NA_real_
      r2   <- if (!is.null(m$r2_pseudo)) m$r2_pseudo else NA_real_
      linhas_info[[i]] <- list(
        Modelo = nomes_modelos[i],
        k      = as.integer(k_p),
        LogLik = ll,
        AIC    = aic,
        AICc   = aicc,
        BIC    = bic,
        RMSE   = rmse,
        R2     = r2
      )
      next
    }
    # Caso venha de estatR com $modelo embutido
    if (is.list(m) && !is.null(m$modelo)) m_core <- m$modelo else m_core <- m

    ll   <- tryCatch(as.numeric(stats::logLik(m_core)), error = function(e) NA_real_)
    aic  <- tryCatch(stats::AIC(m_core), error = function(e) if (!is.null(m_core$aic)) m_core$aic else NA_real_)
    bic  <- tryCatch(stats::BIC(m_core), error = function(e) if (!is.null(m_core$bic)) m_core$bic else NA_real_)
    k_p  <- tryCatch(attr(stats::logLik(m_core), "df"), error = function(e) length(stats::coef(m_core)))
    if (is.null(k_p) || length(k_p) == 0L || is.na(k_p)) k_p <- 2L

    res  <- tryCatch(as.numeric(stats::residuals(m_core)), error = function(e) NULL)
    n_obs <- if (!is.null(res)) length(res[!is.na(res)]) else tryCatch(stats::nobs(m_core), error = function(e) 100L)

    aicc <- if (!is.na(aic) && !is.na(k_p) && (n_obs - k_p - 1L) > 0) {
      aic + (2 * k_p * (k_p + 1L)) / (n_obs - k_p - 1L)
    } else aic

    rmse <- if (!is.null(res)) sqrt(mean(res^2, na.rm = TRUE)) else NA_real_

    r2 <- tryCatch({
      s <- summary(m_core)
      if (!is.null(s$r.squared)) s$r.squared
      else if (!is.null(s$pseudo.r.squared)) s$pseudo.r.squared
      else NA_real_
    }, error = function(e) NA_real_)

    linhas_info[[i]] <- list(
      Modelo = nomes_modelos[i],
      k      = as.integer(k_p),
      LogLik = ll,
      AIC    = aic,
      AICc   = aicc,
      BIC    = bic,
      RMSE   = rmse,
      R2     = r2
    )
  }

  df_comp <- do.call(rbind.data.frame, linhas_info)

  # Ordenacao
  decrescente <- (criterio == "r2")
  col_ord <- switch(criterio,
    "aicc" = "AICc", "aic" = "AIC", "bic" = "BIC",
    "rmse" = "RMSE", "r2" = "R2"
  )

  ord_idx <- order(df_comp[[col_ord]], decreasing = decrescente, na.last = TRUE)
  df_comp <- df_comp[ord_idx, ]

  # Delta AICc e Pesos de Akaike
  if (all(!is.na(df_comp$AICc))) {
    min_aicc <- min(df_comp$AICc, na.rm = TRUE)
    df_comp$Delta_AICc <- df_comp$AICc - min_aicc
    w_raw <- exp(-0.5 * df_comp$Delta_AICc)
    df_comp$Peso_Akaike <- w_raw / sum(w_raw)
  } else {
    df_comp$Delta_AICc  <- NA_real_
    df_comp$Peso_Akaike <- NA_real_
  }

  df_comp$Ranking <- paste0(seq_len(nrow(df_comp)), "º")

  # Formatacao para impressao
  pad <- function(s, w, align = "left") {
    s <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right") return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }
  sep_l <- function(w) paste0("  ", strrep("─", w))

  .print_titulo("COMPARAÇÃO E RANQUEAMENTO DE MODELOS")
  cat(sprintf("  Quantidade de Modelos: %d\n", nrow(df_comp)))
  cat(sprintf("  Critério de Seleção:  %s (%s)\n\n", toupper(criterio), if (decrescente) "Maior é melhor" else "Menor é melhor"))

  .print_topico(sprintf("RANKING GERAL (Ordenado por %s)", toupper(criterio)))

  w_cols <- c(rank = 7, mod = 20, k = 5, aicc = 10, bic = 10, delta = 10, peso = 10, rmse = 10)
  w_total <- sum(w_cols) + length(w_cols) - 1

  hdr_txt <- paste0(
    pad("Rank",   w_cols["rank"],  "center"),
    pad("Modelo", w_cols["mod"],   "left"),
    pad("k",      w_cols["k"],     "center"),
    pad("AICc",   w_cols["aicc"],  "center"),
    pad("BIC",    w_cols["bic"],   "center"),
    pad("Δ AICc",w_cols["delta"], "center"),
    pad("Peso Wi",w_cols["peso"],  "center"),
    pad("RMSE",   w_cols["rmse"],  "center")
  )

  cat(sep_l(w_total), "\n")
  cat("  ", hdr_txt, "\n", sep = "")
  cat(sep_l(w_total), "\n")

  for (r in seq_len(nrow(df_comp))) {
    aicc_s  <- if (is.na(df_comp$AICc[r])) "-" else .fmt_num(df_comp$AICc[r], 1)
    bic_s   <- if (is.na(df_comp$BIC[r])) "-" else .fmt_num(df_comp$BIC[r], 1)
    delta_s <- if (is.na(df_comp$Delta_AICc[r])) "-" else .fmt_num(df_comp$Delta_AICc[r], 2)
    peso_s  <- if (is.na(df_comp$Peso_Akaike[r])) "-" else sprintf("%.1f%%", df_comp$Peso_Akaike[r] * 100)
    rmse_s  <- if (is.na(df_comp$RMSE[r])) "-" else .fmt_num(df_comp$RMSE[r], 3)

    cat("  ",
        pad(df_comp$Ranking[r], w_cols["rank"],  "center"),
        pad(df_comp$Modelo[r],  w_cols["mod"],   "left"),
        pad(df_comp$k[r],       w_cols["k"],     "center"),
        pad(aicc_s,             w_cols["aicc"],  "center"),
        pad(bic_s,              w_cols["bic"],   "center"),
        pad(delta_s,            w_cols["delta"], "center"),
        pad(peso_s,             w_cols["peso"],  "center"),
        pad(rmse_s,             w_cols["rmse"],  "center"),
        "\n", sep = "")
  }
  cat(sep_l(w_total), "\n\n")

  .print_topico("INTERPRETAÇÃO")
  melhor_mod <- df_comp$Modelo[1]
  cat(sprintf("  O modelo '%s' apresentou o melhor desempenho segundo o critério %s.\n", melhor_mod, toupper(criterio)))
  if (!is.na(df_comp$Peso_Akaike[1])) {
    cat(sprintf("  A probabilidade de '%s' ser o melhor modelo entre os avaliados é de %.1f%% (Peso de Akaike).\n",
                melhor_mod, df_comp$Peso_Akaike[1] * 100))
  }
  cat("  Modelos com Δ AICc < 2 possuem suporte substancial; Δ AICc > 10 não possuem suporte.\n\n")
  .print_rodape()

  res <- list(
    ranking     = df_comp,
    melhor      = melhor_mod,
    criterio    = criterio,
    modelos     = modelos
  )
  class(res) <- c("estatR_comparacao_modelos", "list")
  invisible(res)
}
