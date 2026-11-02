# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# MÃ“DULO: VALIDAÃ‡ÃƒO DE DADOS
# Regras de integridade, verificaÃ§Ã£o e diagnÃ³stico de consistÃªncia
# â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

# â”€â”€ REGRAS DE VALIDAÃ‡ÃƒO (FUNÃ‡Ã•ES AUXILIARES EXPORTADAS) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' Regra: Valor entre mÃ­nimo e mÃ¡ximo
#' @param min Limite inferior (inclusivo)
#' @param max Limite superior (inclusivo)
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
entre <- function(min, max) {
  structure(
    list(
      nome = sprintf("entre(%s, %s)", min, max),
      fn = function(x) !is.na(x) & x >= min & x <= max
    ),
    class = "regra_validacao"
  )
}

#' Regra: Valor estritamente maior que um limiar
#' @param x Limiar mÃ­nimo
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
maior_que <- function(x) {
  structure(
    list(
      nome = sprintf("maior_que(%s)", x),
      fn = function(val) !is.na(val) & val > x
    ),
    class = "regra_validacao"
  )
}

#' Regra: Valor estritamente menor que um limiar
#' @param x Limiar mÃ¡ximo
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
menor_que <- function(x) {
  structure(
    list(
      nome = sprintf("menor_que(%s)", x),
      fn = function(val) !is.na(val) & val < x
    ),
    class = "regra_validacao"
  )
}

#' Regra: Sem valores nulos (ausentes / NA)
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
sem_nulos <- function() {
  structure(
    list(
      nome = "sem_nulos()",
      fn = function(x) !is.na(x)
    ),
    class = "regra_validacao"
  )
}

#' Regra: Sem valores duplicados na coluna
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
sem_duplicatas <- function() {
  structure(
    list(
      nome = "sem_duplicatas()",
      fn = function(x) !duplicated(x) & !duplicated(x, fromLast = TRUE)
    ),
    class = "regra_validacao"
  )
}

#' Regra: Valores pertencentes a um conjunto de nÃ­veis permitidos
#' @param v Vetor de valores permitidos
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
niveis_permitidos <- function(v) {
  rotulo <- if (length(v) <= 3) paste(v, collapse = ", ") else paste0(paste(v[1:3], collapse = ", "), ", ...")
  structure(
    list(
      nome = sprintf("niveis_permitidos(c(%s))", rotulo),
      fn = function(x) !is.na(x) & (x %in% v)
    ),
    class = "regra_validacao"
  )
}

#' Regra: Valor corresponde a um padrÃ£o Regex
#' @param regex ExpressÃ£o regular
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
padrao <- function(regex) {
  structure(
    list(
      nome = sprintf("padrao(\"%s\")", regex),
      fn = function(x) !is.na(x) & grepl(regex, as.character(x))
    ),
    class = "regra_validacao"
  )
}

#' Regra: VariÃ¡vel deve ser de tipo numÃ©rico
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
tipo_numerico <- function() {
  structure(
    list(
      nome = "tipo_numerico()",
      fn = function(x) is.numeric(x)
    ),
    class = "regra_validacao"
  )
}

#' Regra: VariÃ¡vel deve ser de tipo texto / categÃ³rica
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
tipo_texto <- function() {
  structure(
    list(
      nome = "tipo_texto()",
      fn = function(x) is.character(x) || is.factor(x)
    ),
    class = "regra_validacao"
  )
}

#' Regra: Comprimento do texto entre a e b caracteres
#' @param a Comprimento mÃ­nimo
#' @param b Comprimento mÃ¡ximo
#' @return Objeto de regra de validaÃ§Ã£o
#' @export
comprimento_entre <- function(a, b) {
  structure(
    list(
      nome = sprintf("comprimento_entre(%d, %d)", a, b),
      fn = function(x) !is.na(x) & nchar(as.character(x)) >= a & nchar(as.character(x)) <= b
    ),
    class = "regra_validacao"
  )
}


# â”€â”€ FUNÃ‡ÃƒO PRINCIPAL: VALIDAR_DADOS â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€

#' @title ValidaÃ§Ã£o de ConsistÃªncia e Integridade de Dados
#' @description Valida um conjunto de dados segundo regras lÃ³gicas expressas
#'   em formato de fÃ³rmula (\code{variavel ~ regra}). Identifica observaÃ§Ãµes violadas,
#'   quantifica inconsistÃªncias e gera relatÃ³rio formatado de auditoria.
#'
#' @param dados Um \code{data.frame} contendo as observaÃ§Ãµes a validar.
#' @param ... FÃ³rmulas no formato \code{variavel ~ regra}, onde regra pode ser
#'   \code{entre()}, \code{maior_que()}, \code{menor_que()}, \code{sem_nulos()},
#'   \code{sem_duplicatas()}, \code{niveis_permitidos()}, \code{padrao()}, etc.
#'
#' @return Retorna um objeto invisÃ­vel da classe \code{estatR_validacao} com
#'   a tabela de violaÃ§Ãµes, lista de linhas afetadas e indicador booleano de conformidade.
#' @export
validar_dados <- function(dados, ...) {
  if (!is.data.frame(dados)) {
    stop("O argumento 'dados' deve ser um data.frame.", call. = FALSE)
  }

  regras_args <- list(...)
  if (length(regras_args) == 0L) {
    stop("Nenhuma regra de valida\u00e7\u00e3o fornecida. Use formulas do tipo: var ~ regra().", call. = FALSE)
  }

  n_obs <- nrow(dados)

  # Helpers locais de formataÃ§Ã£o
  pad <- function(s, w, align = "left") {
    s <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right")  return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }

  sep_linha <- function(w) paste0("  ", strrep("\u2500", w))

  var_nomes       <- character()
  regra_nomes     <- character()
  n_violacoes     <- integer()
  linhas_violadas <- list()
  status_regra    <- character()

  for (i in seq_along(regras_args)) {
    item <- regras_args[[i]]
    if (!inherits(item, "formula")) {
      stop(sprintf("O argumento %d deve ser uma formula no padrao: variavel ~ regra().", i), call. = FALSE)
    }

    var_sym <- item[[2]]
    var_str <- as.character(var_sym)

    if (!var_str %in% names(dados)) {
      var_nomes       <- c(var_nomes, var_str)
      regra_nomes     <- c(regra_nomes, "Coluna n\u00e3o encontrada")
      n_violacoes     <- c(n_violacoes, n_obs)
      linhas_violadas[[length(linhas_violadas) + 1L]] <- seq_len(n_obs)
      status_regra    <- c(status_regra, "Erro: Coluna inexistente")
      next
    }

    # Avalia a regra
    regra_obj <- eval(item[[3]], envir = parent.frame())
    if (!inherits(regra_obj, "regra_validacao")) {
      stop(sprintf("O lado direito da formula %d n\u00e3o produziu uma regra v\u00e1lida.", i), call. = FALSE)
    }

    col_val <- dados[[var_str]]
    valido <- tryCatch(regra_obj$fn(col_val), error = function(e) rep(FALSE, n_obs))

    if (length(valido) == 1L && length(valido) != n_obs) {
      valido <- rep(valido, n_obs)
    }

    violadas <- which(!valido)
    n_v      <- length(violadas)

    var_nomes       <- c(var_nomes, var_str)
    regra_nomes     <- c(regra_nomes, regra_obj$nome)
    n_violacoes     <- c(n_violacoes, n_v)
    linhas_violadas[[length(linhas_violadas) + 1L]] <- violadas
    status_regra    <- c(status_regra, if (n_v == 0L) "[OK] Conforme" else "[!] Viola\u00e7\u00e3o")
  }

  total_violacoes <- sum(n_violacoes)
  todas_linhas_afetadas <- sort(unique(unlist(linhas_violadas)))
  n_linhas_afetadas <- length(todas_linhas_afetadas)
  vars_com_problema <- unique(var_nomes[n_violacoes > 0L])

  # â”€â”€ SaÃ­da no Console â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  cat("\n")
  .print_header("VALIDA\u00c7\u00c3O DE DADOS")

  # ConfiguraÃ§Ã£o
  if (exists(".print_topico", mode = "function")) {
    .print_topico("CONFIGURA\u00c7\u00c3O")
  } else {
    cat("  ▶ CONFIGURAÇÃO\n")
  }
  cat(sprintf("  Observa\u00e7\u00f5es analisadas : %d\n", n_obs))
  cat(sprintf("  Regras verificadas     : %d\n\n", length(regras_args)))

  # Tabela de Resultados
  if (exists(".print_topico", mode = "function")) {
    .print_topico("RESULTADOS")
  } else {
    cat("  ▶ RESULTADOS\n")
  }

  linhas_resumo <- character(length(var_nomes))
  for (k in seq_along(var_nomes)) {
    lh <- linhas_violadas[[k]]
    if (length(lh) == 0L) {
      linhas_resumo[k] <- "\u2014"
    } else if (length(lh) <= 8L) {
      linhas_resumo[k] <- paste0("{", paste(lh, collapse = ", "), "}")
    } else {
      linhas_resumo[k] <- paste0("{", paste(lh[1:8], collapse = ", "), ", ... e mais ", length(lh) - 8L, "}")
    }
  }

  w_tb <- c(var = 14, regra = 28, viol = 12, linhas = 26)
  w_tot <- sum(w_tb) + 3

  hdr_t <- paste0(
    pad("Vari\u00e1vel",     w_tb["var"],    "left"),
    pad("Regra",        w_tb["regra"],  "left"),
    pad("N Viola\u00e7\u00f5es", w_tb["viol"],   "right"), " ",
    pad("Linhas afetadas", w_tb["linhas"], "left")
  )

  cat(sep_linha(w_tot), "\n")
  cat("  ", hdr_t, "\n", sep = "")
  cat(sep_linha(w_tot), "\n")

  for (k in seq_along(var_nomes)) {
    cat("  ",
        pad(var_nomes[k],        w_tb["var"],    "left"),
        pad(regra_nomes[k],      w_tb["regra"],  "left"),
        pad(n_violacoes[k],      w_tb["viol"],   "right"), " ",
        pad(linhas_resumo[k],    w_tb["linhas"], "left"),
        "\n", sep = "")
  }
  cat(sep_linha(w_tot), "\n\n")

  # Resumo
  if (exists(".print_topico", mode = "function")) {
    .print_topico("RESUMO")
  } else {
    cat("  ▶ RESUMO\n")
  }
  cat(sprintf("  Total de viola\u00e7\u00f5es : %d\n", total_violacoes))
  cat(sprintf("  Linhas afetadas    : %d (%.1f%% do banco)\n", n_linhas_afetadas, (n_linhas_afetadas / n_obs) * 100))
  cat(sprintf("  Vari\u00e1veis afetadas : %d de %d\n", length(vars_com_problema), length(unique(var_nomes))))

  if (total_violacoes == 0L) {
    cat("  Status             : [OK] Todos os dados est\u00e3o conformes!\n\n")
  } else {
    cat("  Status             : [!] Dados com problemas detectados\n\n")
  }

  .print_footer()

  res_df <- data.frame(
    variavel        = var_nomes,
    regra           = regra_nomes,
    n_violacoes     = as.integer(n_violacoes),
    linhas_afetadas = I(linhas_violadas),
    status          = status_regra,
    stringsAsFactors = FALSE
  )

  res <- list(
    violacoes             = res_df,
    n_total_violacoes     = total_violacoes,
    linhas_afetadas_total = todas_linhas_afetadas,
    n_linhas_afetadas     = n_linhas_afetadas,
    variaveis_afetadas    = vars_com_problema,
    valido                = (total_violacoes == 0L),
    dados                 = dados
  )

  class(res) <- c("estatR_validacao", "list")
  invisible(res)
}
