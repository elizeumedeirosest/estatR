# MÃ³dulo de imputaÃ§Ã£o de dados para o pacote estatR
# Implementa 17 mÃ©todos de imputaÃ§Ã£o com saÃ­da formatada no console

# â”€â”€ FUNÃ‡Ã•ES AUXILIARES INTERNAS â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Calcula a moda de um vetor (retorna o primeiro em caso de empate)
#' @param x Vetor (numÃ©rico, character ou factor)
#' @noRd
.calcular_moda_imp <- function(x) {
  x_sem_na <- x[!is.na(x)]
  if (length(x_sem_na) == 0L) return(NA)
  tab  <- table(x_sem_na)
  nome <- names(tab)[which.max(tab)]
  # Preserva o tipo original
  if (is.numeric(x)) as.numeric(nome) else nome
}

#' Last Observation Carried Forward (LOCF) em um vetor
#' @param x Vetor com possÃ­veis NAs
#' @noRd
.imputar_locf <- function(x) {
  ultimo_valido <- NA
  for (i in seq_along(x)) {
    if (!is.na(x[i])) {
      ultimo_valido <- x[i]
    } else if (!is.na(ultimo_valido)) {
      x[i] <- ultimo_valido
    }
  }
  x
}

#' Next Observation Carried Backward (NOCB) em um vetor
#' @param x Vetor com possÃ­veis NAs
#' @noRd
.imputar_nocb <- function(x) {
  rev(.imputar_locf(rev(x)))
}

#' ImputaÃ§Ã£o KNN manual com distÃ¢ncia euclidiana normalizada
#' @param dados Data frame com colunas numÃ©ricas a imputar
#' @param k NÃºmero de vizinhos mais prÃ³ximos
#' @noRd
.imputar_knn <- function(dados, k = 5) {
  dados_num <- dados[, vapply(dados, is.numeric, logical(1L)), drop = FALSE]
  if (ncol(dados_num) == 0L) return(dados)

  n  <- nrow(dados_num)
  nc <- ncol(dados_num)

  # Normaliza cada coluna pelo desvio padrÃ£o (ignora NAs)
  sds <- vapply(dados_num, function(col) {
    s <- stats::sd(col, na.rm = TRUE)
    if (is.na(s) || s == 0) 1 else s
  }, numeric(1L))

  dados_norm <- sweep(dados_num, 2L, sds, "/")

  for (j in seq_len(nc)) {
    col_orig <- names(dados_num)[j]
    idx_na   <- which(is.na(dados_num[[j]]))
    if (length(idx_na) == 0L) next

    for (i in idx_na) {
      # Linhas com valor observado na coluna j e pelo menos uma coluna em comum
      linhas_com_valor <- which(!is.na(dados_num[[j]]))
      if (length(linhas_com_valor) == 0L) next

      # DistÃ¢ncia euclidiana nas demais colunas (excluindo j)
      outras_cols <- setdiff(seq_len(nc), j)
      if (length(outras_cols) == 0L) {
        # Somente uma coluna numÃ©rica: usa mÃ©dia
        dados[[col_orig]][i] <- mean(dados_num[[j]], na.rm = TRUE)
        next
      }

      linha_i <- as.numeric(dados_norm[i, outras_cols, drop = FALSE])
      dists <- vapply(linhas_com_valor, function(l) {
        linha_l <- as.numeric(dados_norm[l, outras_cols, drop = FALSE])
        # Pares com ambos observados
        mask <- !is.na(linha_i) & !is.na(linha_l)
        if (!any(mask)) return(Inf)
        sqrt(sum((linha_i[mask] - linha_l[mask])^2))
      }, numeric(1L))

      # Seleciona os k vizinhos mais prÃ³ximos (finitos)
      k_ef  <- min(k, sum(is.finite(dists)))
      if (k_ef == 0L) next
      idx_viz <- linhas_com_valor[order(dists)[seq_len(k_ef)]]
      dados[[col_orig]][i] <- mean(dados_num[[j]][idx_viz], na.rm = TRUE)
    }
  }
  dados
}

# â”€â”€ FUNÃ‡ÃƒO PRINCIPAL â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' @title ImputaÃ§Ã£o de Dados Ausentes
#' @description Realiza imputaÃ§Ã£o de valores ausentes (NA) em um data frame com
#'   17 mÃ©todos disponÃ­veis: estatÃ­sticas descritivas simples e por grupo,
#'   mÃ©todos de sÃ©ries temporais (LOCF/NOCB), interpolaÃ§Ã£o, regressÃ£o linear,
#'   k-vizinhos mais prÃ³ximos (KNN) e imputaÃ§Ã£o mÃºltipla (MICE). O resultado Ã©
#'   exibido com diagnÃ³stico completo antes e depois da imputaÃ§Ã£o.
#'
#' @param dados Data frame contendo os dados a serem imputados.
#' @param metodo Character. MÃ©todo de imputaÃ§Ã£o a utilizar. OpÃ§Ãµes:
#'   \code{"media"}, \code{"mediana"}, \code{"moda"}, \code{"min"}, \code{"max"},
#'   \code{"media_grupo"}, \code{"mediana_grupo"}, \code{"moda_grupo"},
#'   \code{"min_grupo"}, \code{"max_grupo"},
#'   \code{"locf"}, \code{"nocb"},
#'   \code{"interpolacao_linear"}, \code{"interpolacao_spline"},
#'   \code{"regressao"}, \code{"knn"}, \code{"mice"}.
#' @param variaveis Character vector opcional. Nomes das colunas a imputar.
#'   Se \code{NULL} (padrÃ£o), todas as colunas com pelo menos um NA sÃ£o usadas.
#' @param grupo Character. Nome da coluna de agrupamento, necessÃ¡rio para os
#'   mÃ©todos \code{*_grupo}.
#' @param k Inteiro. NÃºmero de vizinhos para o mÃ©todo \code{"knn"} (padrÃ£o: 5).
#' @param n_imp Inteiro. NÃºmero de imputaÃ§Ãµes para o mÃ©todo \code{"mice"} (padrÃ£o: 5).
#' @param semente Inteiro opcional. Semente aleatÃ³ria para reprodutibilidade.
#'
#' @return Objeto invisÃ­vel da classe \code{c("estatR_imputacao", "list")} com:
#'   \item{dados_imputados}{Data frame com os valores ausentes imputados.}
#'   \item{relatorio}{Data frame com variÃ¡vel, n_ausentes_antes, n_ausentes_depois
#'     e n_imputados.}
#'   \item{metodo}{MÃ©todo de imputaÃ§Ã£o utilizado.}
#'   \item{variaveis_imputadas}{Nomes das variÃ¡veis efetivamente imputadas.}
#'   \item{n_total_imputados}{Total de valores imputados em todas as variÃ¡veis.}
#'   \item{grupo}{Nome da coluna de agrupamento (ou NULL).}
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' dados_ex <- data.frame(
#'   x = c(1, 2, NA, 4, 5),
#'   y = c(NA, 2, 3, 4, NA),
#'   g = c("A", "A", "B", "B", "A")
#' )
#'
#' # ImputaÃ§Ã£o pela mÃ©dia
#' res <- imputar(dados_ex, metodo = "media")
#' res$dados_imputados
#'
#' # ImputaÃ§Ã£o pela mediana por grupo
#' res2 <- imputar(dados_ex, metodo = "mediana_grupo", grupo = "g")
#' res2$dados_imputados
imputar <- function(dados,
                    metodo    = c("media", "mediana", "moda", "min", "max",
                                  "media_grupo", "mediana_grupo", "moda_grupo",
                                  "min_grupo", "max_grupo",
                                  "locf", "nocb",
                                  "interpolacao_linear", "interpolacao_spline",
                                  "regressao", "knn", "mice"),
                    variaveis = NULL,
                    grupo     = NULL,
                    k         = 5L,
                    n_imp     = 5L,
                    semente   = NULL) {

  # â”€â”€ 1. ValidaÃ§Ãµes â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  if (!is.data.frame(dados))
    stop("'dados' deve ser um data.frame.")

  metodo <- match.arg(metodo)

  # MÃ©todos que exigem coluna de grupo
  metodos_grupo <- c("media_grupo", "mediana_grupo", "moda_grupo",
                     "min_grupo",  "max_grupo")
  if (metodo %in% metodos_grupo && is.null(grupo))
    stop("O argumento 'grupo' \u00e9 obrigat\u00f3rio para o m\u00e9todo '", metodo, "'.")

  if (!is.null(grupo) && !grupo %in% names(dados))
    stop("A coluna de grupo '", grupo, "' n\u00e3o foi encontrada em 'dados'.")

  if (!is.null(semente)) set.seed(semente)

  dados_novo <- dados

  # â”€â”€ 2. SeleÃ§Ã£o de variÃ¡veis a imputar â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  if (is.null(variaveis)) {
    variaveis <- names(dados)[vapply(dados, function(col) anyNA(col), logical(1L))]
  } else {
    ausentes <- setdiff(variaveis, names(dados))
    if (length(ausentes) > 0L)
      stop("Vari\u00e1veis n\u00e3o encontradas em 'dados': ",
           paste(ausentes, collapse = ", "))
    variaveis <- variaveis[vapply(variaveis,
                                  function(v) anyNA(dados[[v]]), logical(1L))]
  }

  # Remove a coluna de grupo da lista de variÃ¡veis a imputar
  if (!is.null(grupo)) variaveis <- setdiff(variaveis, grupo)

  if (length(variaveis) == 0L) {
    message("Nenhuma vari\u00e1vel com valores ausentes encontrada.")
    res <- list(
      dados_imputados     = dados,
      relatorio           = data.frame(variavel = character(0),
                                       n_ausentes_antes  = integer(0),
                                       n_ausentes_depois = integer(0),
                                       n_imputados       = integer(0)),
      metodo              = metodo,
      variaveis_imputadas = character(0),
      n_total_imputados   = 0L,
      grupo               = grupo
    )
    class(res) <- c("estatR_imputacao", "list")
    return(invisible(res))
  }

  # â”€â”€ 3. DiagnÃ³stico antes da imputaÃ§Ã£o â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  n_antes <- vapply(variaveis, function(v) sum(is.na(dados[[v]])), integer(1L))
  n_obs   <- nrow(dados)

  # â”€â”€ 4. Tipo de cada variÃ¡vel â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  tipo_var <- vapply(variaveis, function(v) {
    if (is.numeric(dados[[v]])) "Num\u00e9rica" else "Categ\u00f3rica"
  }, character(1L))

  # â”€â”€ 5. Executa a imputaÃ§Ã£o â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

  # Helper: estatÃ­stica simples (toda a coluna)
  .stat_simples <- function(col, stat) {
    switch(stat,
      media   = mean(col, na.rm = TRUE),
      mediana = stats::median(col, na.rm = TRUE),
      moda    = .calcular_moda_imp(col),
      min     = min(col, na.rm = TRUE),
      max     = max(col, na.rm = TRUE)
    )
  }

  # RÃ³tulo amigÃ¡vel do mÃ©todo (para exibiÃ§Ã£o)
  .rotulo_metodo <- function(met, grp = NULL) {
    base <- switch(met,
      media                = "m\u00e9dia",
      mediana              = "mediana",
      moda                 = "moda",
      min                  = "m\u00ednimo",
      max                  = "m\u00e1ximo",
      media_grupo          = paste0("m\u00e9dia por grupo (", grp, ")"),
      mediana_grupo        = paste0("mediana por grupo (", grp, ")"),
      moda_grupo           = paste0("moda por grupo (", grp, ")"),
      min_grupo            = paste0("m\u00ednimo por grupo (", grp, ")"),
      max_grupo            = paste0("m\u00e1ximo por grupo (", grp, ")"),
      locf                 = "LOCF",
      nocb                 = "NOCB",
      interpolacao_linear  = "interpola\u00e7\u00e3o linear",
      interpolacao_spline  = "interpola\u00e7\u00e3o spline",
      regressao            = "regress\u00e3o linear (lm)",
      knn                  = paste0("KNN (k = ", k, ")"),
      mice                 = paste0("MICE (", n_imp, " imputa\u00e7\u00f5es)")
    )
    base
  }

  # --- mÃ©todos simples: media, mediana, moda, min, max -------------------------
  if (metodo %in% c("media", "mediana", "moda", "min", "max")) {
    for (v in variaveis) {
      col     <- dados_novo[[v]]
      idx_na  <- which(is.na(col))
      if (is.numeric(col)) {
        val <- .stat_simples(col, metodo)
      } else {
        # VariÃ¡veis character/factor: sempre moda
        val <- .calcular_moda_imp(col)
      }
      dados_novo[[v]][idx_na] <- val
    }

  # --- mÃ©todos por grupo -------------------------------------------------------
  } else if (metodo %in% metodos_grupo) {
    stat_base <- sub("_grupo$", "", metodo)   # ex: "media_grupo" -> "media"

    for (v in variaveis) {
      col     <- dados_novo[[v]]
      grp_col <- dados_novo[[grupo]]
      idx_na  <- which(is.na(col))
      if (length(idx_na) == 0L) next

      for (i in idx_na) {
        grp_val   <- grp_col[i]
        mask_grp  <- grp_col == grp_val & !is.na(grp_col)
        col_grp   <- col[mask_grp]

        if (is.numeric(col)) {
          val <- if (all(is.na(col_grp))) .stat_simples(col, stat_base)
                 else .stat_simples(col_grp, stat_base)
        } else {
          val <- if (all(is.na(col_grp))) .calcular_moda_imp(col)
                 else .calcular_moda_imp(col_grp)
        }
        dados_novo[[v]][i] <- val
      }
    }

  # --- LOCF -------------------------------------------------------------------
  } else if (metodo == "locf") {
    for (v in variaveis) {
      dados_novo[[v]] <- .imputar_locf(dados_novo[[v]])
    }

  # --- NOCB -------------------------------------------------------------------
  } else if (metodo == "nocb") {
    for (v in variaveis) {
      dados_novo[[v]] <- .imputar_nocb(dados_novo[[v]])
    }

  # --- InterpolaÃ§Ã£o linear -----------------------------------------------------
  } else if (metodo == "interpolacao_linear") {
    for (v in variaveis) {
      col    <- dados_novo[[v]]
      if (!is.numeric(col)) {
        # Para nÃ£o-numÃ©ricas: moda
        dados_novo[[v]][is.na(col)] <- .calcular_moda_imp(col)
        next
      }
      idx_ok <- which(!is.na(col))
      idx_na <- which(is.na(col))
      if (length(idx_ok) < 2L || length(idx_na) == 0L) next
      # InterpolaÃ§Ã£o: approx com regra = 2 (extrapola nas bordas pelo extremo)
      interp <- stats::approx(x = idx_ok, y = col[idx_ok],
                              xout = idx_na, rule = 2)
      dados_novo[[v]][idx_na] <- interp$y
    }

  # --- InterpolaÃ§Ã£o spline -----------------------------------------------------
  } else if (metodo == "interpolacao_spline") {
    for (v in variaveis) {
      col    <- dados_novo[[v]]
      if (!is.numeric(col)) {
        dados_novo[[v]][is.na(col)] <- .calcular_moda_imp(col)
        next
      }
      idx_ok <- which(!is.na(col))
      idx_na <- which(is.na(col))
      if (length(idx_ok) < 3L || length(idx_na) == 0L) {
        # Fallback: linear quando pontos insuficientes para spline
        if (length(idx_ok) >= 2L) {
          interp <- stats::approx(x = idx_ok, y = col[idx_ok],
                                  xout = idx_na, rule = 2)
          dados_novo[[v]][idx_na] <- interp$y
        }
        next
      }
      interp <- stats::spline(x = idx_ok, y = col[idx_ok],
                              xout = idx_na, method = "natural")
      dados_novo[[v]][idx_na] <- interp$y
    }

  # --- RegressÃ£o linear -------------------------------------------------------
  } else if (metodo == "regressao") {
    # Armazena coeficientes dos modelos (para exibiÃ§Ã£o posterior)
    modelos_regressao <<- list()

    vars_num <- variaveis[vapply(variaveis, function(v) is.numeric(dados_novo[[v]]),
                                 logical(1L))]
    vars_cat <- setdiff(variaveis, vars_num)

    # VariÃ¡veis categÃ³ricas: moda
    for (v in vars_cat) {
      col <- dados_novo[[v]]
      dados_novo[[v]][is.na(col)] <- .calcular_moda_imp(col)
    }

    # VariÃ¡veis numÃ©ricas: regressÃ£o
    for (v in vars_num) {
      idx_na <- which(is.na(dados_novo[[v]]))
      if (length(idx_na) == 0L) next

      # Preditores: demais variÃ¡veis numÃ©ricas sem NA nas linhas de treino
      preditores <- setdiff(names(dados_novo)[vapply(dados_novo, is.numeric,
                                                      logical(1L))], v)
      if (length(preditores) == 0L) {
        # Fallback: mÃ©dia
        dados_novo[[v]][idx_na] <- mean(dados_novo[[v]], na.rm = TRUE)
        next
      }

      # Linhas de treino: sem NA em v nem nos preditores
      mask_treino <- !is.na(dados_novo[[v]]) &
        apply(dados_novo[, preditores, drop = FALSE], 1L,
              function(r) !anyNA(r))

      if (sum(mask_treino) < (length(preditores) + 1L)) {
        # Poucos dados: fallback mÃ©dia
        dados_novo[[v]][idx_na] <- mean(dados_novo[[v]], na.rm = TRUE)
        next
      }

      formula_str <- paste(v, "~", paste(preditores, collapse = " + "))
      mod <- tryCatch(
        stats::lm(stats::as.formula(formula_str),
                  data = dados_novo[mask_treino, , drop = FALSE]),
        error = function(e) NULL
      )

      if (is.null(mod)) {
        dados_novo[[v]][idx_na] <- mean(dados_novo[[v]], na.rm = TRUE)
        next
      }

      modelos_regressao[[v]] <<- mod

      # Linhas de previsÃ£o: sem NA nos preditores
      mask_pred <- idx_na[apply(dados_novo[idx_na, preditores, drop = FALSE],
                                1L, function(r) !anyNA(r))]
      if (length(mask_pred) > 0L) {
        pred <- stats::predict(mod, newdata = dados_novo[mask_pred, , drop = FALSE])
        dados_novo[[v]][mask_pred] <- pred
      }

      # Eventuais NA restantes (preditores ausentes): mÃ©dia
      ainda_na <- which(is.na(dados_novo[[v]]))
      if (length(ainda_na) > 0L)
        dados_novo[[v]][ainda_na] <- mean(dados_novo[[v]], na.rm = TRUE)
    }

  # --- KNN manual -------------------------------------------------------------
  } else if (metodo == "knn") {
    dados_novo <- .imputar_knn(dados_novo, k = k)

    # VariÃ¡veis categÃ³ricas: moda
    for (v in variaveis) {
      col <- dados_novo[[v]]
      if (!is.numeric(col) && anyNA(col))
        dados_novo[[v]][is.na(col)] <- .calcular_moda_imp(dados[[v]])
    }

  # --- MICE -------------------------------------------------------------------
  } else if (metodo == "mice") {
    if (!requireNamespace("mice", quietly = TRUE))
      stop("O pacote 'mice' n\u00e3o est\u00e1 instalado. ",
           "Instale com: install.packages('mice')")

    res_mice <- mice::mice(dados_novo, m = n_imp, seed = semente,
                           printFlag = FALSE)
    dados_novo <- mice::complete(res_mice, action = 1L)
  }

  # â”€â”€ 6. DiagnÃ³stico depois da imputaÃ§Ã£o â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  n_depois <- vapply(variaveis, function(v) sum(is.na(dados_novo[[v]])),
                     integer(1L))
  n_imput  <- n_antes - n_depois

  # â”€â”€ 7. SaÃ­da no console â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

  # Helpers locais de formataÃ§Ã£o de tabela
  pad <- function(s, w, align = "left") {
    s  <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0L) return(s)
    if (align == "right")  return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp / 2L)), s,
                                         strrep(" ", ceiling(sp / 2L))))
    paste0(s, strrep(" ", sp))
  }

  sep_linha <- function(w) paste0("  ", strrep("\u2500", w))

  # RÃ³tulo do mÃ©todo para exibiÃ§Ã£o
  rotulo <- .rotulo_metodo(metodo, grupo)
  rotulo_cap <- paste0(toupper(substring(rotulo, 1L, 1L)),
                       substring(rotulo, 2L))

  # CabeÃ§alho
  cat("\n")
  .print_titulo("IMPUTA\u00c7\u00c3O DE DADOS")

  # â”€â”€ ConfiguraÃ§Ã£o â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  .print_topico("CONFIGURA\u00c7\u00c3O")
  cat(sprintf("  M\u00e9todo          : %s\n", rotulo_cap))
  if (!is.null(grupo))
    cat(sprintf("  Vari\u00e1vel grupo  : %s\n", grupo))
  cat(sprintf("  Vari\u00e1veis       : %s\n", paste(variaveis, collapse = ", ")))
  cat(sprintf("  Observa\u00e7\u00f5es     : %d\n", n_obs))
  cat("\n")

  # â”€â”€ DiagnÃ³stico antes â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  .print_topico("DIAGN\u00d3STICO DE VALORES AUSENTES (ANTES)")

  w_d <- c(var = 20L, tipo = 14L, naus = 14L, pct = 12L)
  hdr_d <- paste0(
    pad("Vari\u00e1vel",    w_d["var"],  "left"),
    pad("Tipo",        w_d["tipo"], "left"),
    pad("N ausentes",  w_d["naus"], "right"),
    pad("% ausente",   w_d["pct"],  "right")
  )
  w_tot_d <- sum(w_d)

  cat(sep_linha(w_tot_d), "\n")
  cat("  ", hdr_d, "\n", sep = "")
  cat(sep_linha(w_tot_d), "\n")

  for (i in seq_along(variaveis)) {
    cat("  ",
        pad(variaveis[i],      w_d["var"],  "left"),
        pad(tipo_var[i],       w_d["tipo"], "left"),
        pad(n_antes[i],        w_d["naus"], "right"),
        pad(.fmt_pct(n_antes[i] / n_obs), w_d["pct"], "right"),
        "\n", sep = "")
  }
  cat(sep_linha(w_tot_d), "\n\n")

  # â”€â”€ ImputaÃ§Ã£o realizada â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  .print_topico("IMPUTA\u00c7\u00c3O")

  w_i <- c(var = 20L, imp = 14L, sub = 40L)
  hdr_i <- paste0(
    pad("Vari\u00e1vel",        w_i["var"], "left"),
    pad("Imputados",       w_i["imp"], "right"),
    pad("Substitu\u00eddo por", w_i["sub"], "left")
  )
  w_tot_i <- sum(w_i)

  cat(sep_linha(w_tot_i), "\n")
  cat("  ", hdr_i, "\n", sep = "")
  cat(sep_linha(w_tot_i), "\n")

  for (i in seq_along(variaveis)) {
    # RÃ³tulo do que foi substituÃ­do
    if (!is.numeric(dados[[variaveis[i]]]) &&
        !metodo %in% c("locf", "nocb",
                       "interpolacao_linear", "interpolacao_spline",
                       "regressao", "knn", "mice")) {
      sub_txt <- "moda (vari\u00e1vel categ\u00f3rica)"
    } else {
      sub_txt <- rotulo
    }

    cat("  ",
        pad(variaveis[i],  w_i["var"], "left"),
        pad(n_imput[i],    w_i["imp"], "right"),
        "  ",
        pad(sub_txt,       w_i["sub"] - 2L, "left"),
        "\n", sep = "")
  }
  cat(sep_linha(w_tot_i), "\n\n")

  # â”€â”€ DiagnÃ³stico depois â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  .print_topico("DIAGN\u00d3STICO DE VALORES AUSENTES (DEPOIS)")

  w_p <- c(var = 20L, naus = 14L, pct = 12L)
  hdr_p <- paste0(
    pad("Vari\u00e1vel",   w_p["var"],  "left"),
    pad("N ausentes", w_p["naus"], "right"),
    pad("% ausente",  w_p["pct"],  "right")
  )
  w_tot_p <- sum(w_p)

  cat(sep_linha(w_tot_p), "\n")
  cat("  ", hdr_p, "\n", sep = "")
  cat(sep_linha(w_tot_p), "\n")

  for (i in seq_along(variaveis)) {
    cat("  ",
        pad(variaveis[i],    w_p["var"],  "left"),
        pad(n_depois[i],     w_p["naus"], "right"),
        pad(.fmt_pct(n_depois[i] / n_obs), w_p["pct"], "right"),
        "\n", sep = "")
  }
  cat(sep_linha(w_tot_p), "\n\n")

  # â”€â”€ Notas â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  cat("  Notas:\n")

  tem_cat <- any(!vapply(variaveis, function(v) is.numeric(dados[[v]]),
                         logical(1L)))
  tem_num <- any(vapply(variaveis, function(v) is.numeric(dados[[v]]),
                        logical(1L)))

  if (tem_cat && metodo %in% c("media", "mediana", "min", "max",
                               "media_grupo", "mediana_grupo",
                               "min_grupo", "max_grupo",
                               "regressao", "knn")) {
    cat("  Vari\u00e1veis de texto (character/factor): imputadas pela moda.\n")
  }

  if (metodo == "regressao" && tem_num) {
    cat("  Vari\u00e1veis num\u00e9ricas: imputadas por regress\u00e3o linear (lm).\n")
    cat("  Signific\u00e2ncia: *** p < 0,001   ** p < 0,01   * p < 0,05\n")
  }

  if (metodo %in% c("locf", "nocb")) {
    cat("  NAs nas bordas sem valor anterior/posterior permanecem como NA.\n")
  }

  if (metodo == "mice") {
    cat("  Imputado com o pacote 'mice' (imputaÃ§Ã£o mÃºltipla); retornada a 1Âª imputaÃ§Ã£o completa.\n")
  }

  cat("\n")
  .print_rodape()

  # â”€â”€ 8. Monta e retorna o objeto â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  relatorio <- data.frame(
    variavel          = variaveis,
    n_ausentes_antes  = as.integer(n_antes),
    n_ausentes_depois = as.integer(n_depois),
    n_imputados       = as.integer(n_imput),
    stringsAsFactors  = FALSE
  )

  res <- list(
    dados_imputados     = dados_novo,
    relatorio           = relatorio,
    metodo              = metodo,
    variaveis_imputadas = variaveis,
    n_total_imputados   = sum(n_imput),
    grupo               = grupo
  )
  class(res) <- c("estatR_imputacao", "list")
  invisible(res)
}
