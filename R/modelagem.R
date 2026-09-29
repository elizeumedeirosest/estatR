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
    cat(sprintf("\n\u2500\u2500 DIVIS\u00c3O TREINO / TESTE (TEMPORAL) %s\n", strrep("\u2500", w_sep - 35)))
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
    cat(strrep("\u2500", w_sep), "\n\n")

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
    cat(sprintf("\n\u2500\u2500 DIVIS\u00c3O TREINO / TESTE (ESTRATIFICADA) %s\n", strrep("\u2500", w_sep - 39)))
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
    cat(strrep("\u2500", w_sep), "\n\n")

    return(invisible(list(treino = treino, teste = teste)))
  }

  # ── Modo aleatório (padrão) ────────────────────────────────────────────────
  idx_treino <- sample(1:n, size = floor(n * proporcao))
  idx_teste  <- setdiff(1:n, idx_treino)

  treino <- dados[sort(idx_treino), , drop = FALSE]
  teste  <- dados[sort(idx_teste),  , drop = FALSE]

  w_sep <- 67
  cat(sprintf("\n\u2500\u2500 DIVIS\u00c3O TREINO / TESTE (ALEAT\u00d3RIA) %s\n", strrep("\u2500", w_sep - 36)))
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
  cat(strrep("\u2500", w_sep), "\n\n")

  invisible(list(treino = treino, teste = teste))
}
