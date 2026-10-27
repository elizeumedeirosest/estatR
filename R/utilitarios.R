# Funções auxiliares para formatação e impressão no console

#' Imprime um cabeçalho formatado no console
#' @param titulo Texto do cabeçalho
#' @param largura Largura total da linha (padrão 80)
#' @noRd
.print_header <- function(titulo, largura = 85) {
  tit_formatado <- paste0(" ", titulo, " ")
  n_tracos_inicio <- 2
  n_tracos_fim <- largura - nchar(tit_formatado) - n_tracos_inicio
  
  if (n_tracos_fim < 0) n_tracos_fim <- 0
  
  linha_inicio <- paste0(rep("─", n_tracos_inicio), collapse = "")
  linha_fim <- paste0(rep("─", n_tracos_fim), collapse = "")
  
  cat(paste0("\033[1;34m", linha_inicio, tit_formatado, linha_fim, "\033[0m\n"))
}

#' Imprime um rodapé formatado no console
#' @param largura Largura total da linha
#' @noRd
.print_footer <- function(largura = 85) {
  linha <- paste0(rep("─", largura), collapse = "")
  cat(paste0("\033[1;34m", linha, "\033[0m\n"))
}

#' Formata número com vírgula e casas decimais fixas
#' @noRd
.fmt_num <- function(x, decimais = 2) {
  if (is.na(x) || is.nan(x)) return("-")
  format(round(x, decimais), nsmall = decimais, decimal.mark = ",", big.mark = ".")
}

#' Formata porcentagem
#' @noRd
.fmt_pct <- function(x, decimais = 1) {
  if (is.na(x) || is.nan(x)) return("-")
  paste0(.fmt_num(x * 100, decimais), "%")
}

#' Formata p-valor
#' @noRd
.fmt_p <- function(p) {
  if (is.na(p) || is.nan(p)) return("-")
  if (p < 0.001) return("< 0,001")
  format(round(p, 3), nsmall = 3, decimal.mark = ",")
}

.pad_string <- function(x, width, align = c("right", "left", "center")) {
  align <- match.arg(align)
  x <- as.character(x)
  
  sapply(x, function(val) {
    val_clean <- gsub("\033\\[[0-9;]*m", "", val) # remove cores para contar largura real
    n_espacos <- width - nchar(val_clean)
    if (n_espacos <= 0) return(val)
    
    if (align == "right") {
      paste0(paste0(rep(" ", n_espacos), collapse = ""), val)
    } else if (align == "left") {
      paste0(val, paste0(rep(" ", n_espacos), collapse = ""))
    } else {
      left_p <- floor(n_espacos / 2)
      right_p <- ceiling(n_espacos / 2)
      paste0(paste0(rep(" ", left_p), collapse = ""), val, paste0(rep(" ", right_p), collapse = ""))
    }
  }, USE.NAMES = FALSE)
}

.th_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
  paste0(s, strrep(" ", sp))
}

.th_sep <- function(w) paste0("  ", strrep("\u2500", w))


#' Calcula a Moda de um vetor (numérico ou caractere/fator)
#' @noRd
.calcular_moda <- function(x) {
  x_sem_na <- x[!is.na(x)]
  if (length(x_sem_na) == 0) return("Sem dados")
  
  tab <- table(x_sem_na)
  max_freq <- max(tab)
  
  # Se todos aparecem 1 vez, não há moda
  if (max_freq == 1 && length(tab) > 1) return("Amodal")
  
  modas <- names(tab)[tab == max_freq]
  
  if (length(modas) > 3) {
    return("Multimodal (>3 modas)")
  }
  
  paste(modas, collapse = ", ")
}

#' Calcula Assimetria (Skewness)
#' @noRd
.calcular_assimetria <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA)
  m3 <- sum((x - mean(x))^3) / n
  s3 <- (sqrt(sum((x - mean(x))^2) / n))^3
  if (s3 == 0) return(0)
  m3 / s3
}

#' Calcula Curtose (Kurtosis)
#' @noRd
.calcular_curtose <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) return(NA)
  m4 <- sum((x - mean(x))^4) / n
  s4 <- (sum((x - mean(x))^2) / n)^2
  if (s4 == 0) return(0)
  (m4 / s4) - 3 # Curtose excedente
}

# ── FORMATADORES VISUAIS PADRONIZADOS ────────────────────────────────────────
#' @noRd
.print_titulo <- function(titulo, w_sep = 78) {
  tit_formatado <- paste0("── ", titulo, " ")
  n_tracos <- w_sep - nchar(tit_formatado)
  if (n_tracos < 0) n_tracos <- 0
  linha_completa <- paste0(tit_formatado, strrep("─", n_tracos))
  cat(sprintf("\n\033[1;34m%s\033[0m\n\n", linha_completa))
}

#' @noRd
.print_topico <- function(texto) {
  cat(sprintf("  \033[1m▶ %s\033[0m\n", texto))
}

#' @noRd
.print_rodape <- function(w_sep = 78) {
  cat(sprintf("\033[1;34m%s\033[0m\n\n", strrep("─", w_sep)))
}

# ── SELEÇÃO AUTOMÁTICA DE FUNÇÃO DE LIGAÇÃO (GLMs) ──────────────────────────
#' @noRd
.auto_ligacao <- function(formula, dados, familia_fn, candidatas) {
  resultados <- lapply(candidatas, function(lk) {
    m <- tryCatch(
      suppressWarnings(stats::glm(formula, data = dados, family = familia_fn(link = lk))),
      error = function(e) NULL
    )
    convergiu <- !is.null(m) && isTRUE(m$converged)
    list(
      ligacao = lk,
      aic     = if (convergiu) stats::AIC(m) else NA_real_,
      desvio  = if (convergiu) m$deviance    else NA_real_
    )
  })
  aics <- sapply(resultados, function(r) if (is.na(r$aic)) Inf else r$aic)
  list(ligacao = candidatas[which.min(aics)], resultados = resultados)
}

#' @noRd
.imprimir_selecao_ligacao <- function(resultados, ligacao_escolhida, pad_fn = NULL, sep_fn = NULL) {
  if (exists(".print_topico", mode = "function")) {
    .print_topico("SELEÇÃO DA FUNÇÃO DE LIGAÇÃO")
  } else {
    cat("  SELEÇÃO DA FUNÇÃO DE LIGAÇÃO\n")
  }
  
  if (is.null(pad_fn)) {
    pad_fn <- function(s, w, align = "left") {
      s <- as.character(s)
      sp <- w - nchar(s)
      if (sp <= 0) return(s)
      if (align == "right")  return(paste0(strrep(" ", sp), s))
      if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
      paste0(s, strrep(" ", sp))
    }
  }
  
  w <- c(lig = 18, aic = 16, dev = 16)
  sep_linha <- if (!is.null(sep_fn)) sep_fn(sum(w)) else paste0("  ", strrep("─", sum(w)))
  
  cat(sep_linha, "\n")
  cat("  ", pad_fn("Ligação", w["lig"], "left"),
             pad_fn("AIC",     w["aic"], "center"),
             pad_fn("Desvio",  w["dev"], "center"), "\n", sep = "")
  cat(sep_linha, "\n")
  for (r in resultados) {
    aic_str <- if (!is.na(r$aic))
      formatC(r$aic,    format = "f", digits = 2, decimal.mark = ",") else "—"
    dev_str <- if (!is.na(r$desvio))
      formatC(r$desvio, format = "f", digits = 2, decimal.mark = ",") else "—"
    cat("  ", pad_fn(r$ligacao, w["lig"], "left"),
               pad_fn(aic_str,  w["aic"], "center"),
               pad_fn(dev_str,  w["dev"], "center"), "\n", sep = "")
  }
  cat(sep_linha, "\n")
  cat(sprintf("  * Função de ligação selecionada: %s\n\n", ligacao_escolhida))
}


# ── EXTRATORES DE MODELOS EM PORTUGUÊS ────────────────────────────────────────

#' @title Coeficientes do Modelo
#' @description Extrai os coeficientes estimados de um modelo estatístico (regressão linear, MLG, ANOVA, etc.).
#' @param modelo Objeto de modelo ajustado.
#' @param ... Argumentos adicionais repassados para \code{stats::coef}.
#' @return Vetor numérico ou matriz com os coeficientes estimados.
#' @export
coeficientes <- function(modelo, ...) {
  if (inherits(modelo, c("estatR_logistica", "estatR_poisson", "estatR_nb", "estatR_zeroinfl"))) {
    return(stats::coef(modelo$modelo, ...))
  }
  stats::coef(modelo, ...)
}

#' @title Resíduos do Modelo
#' @description Extrai os resíduos de um modelo ajustado.
#' @param modelo Objeto de modelo ajustado.
#' @param ... Argumentos adicionais repassados para \code{stats::residuals}.
#' @return Vetor numérico com os resíduos.
#' @export
residuos <- function(modelo, ...) {
  if (inherits(modelo, c("estatR_logistica", "estatR_poisson", "estatR_nb", "estatR_zeroinfl"))) {
    return(stats::residuals(modelo$modelo, ...))
  }
  stats::residuals(modelo, ...)
}

#' @title Valores Ajustados do Modelo
#' @description Extrai os valores preditos / ajustados pelo modelo sobre a amostra.
#' @param modelo Objeto de modelo ajustado.
#' @param ... Argumentos adicionais repassados para \code{stats::fitted}.
#' @return Vetor numérico com os valores ajustados.
#' @export
#' @aliases ajustados valores_preditos
valores_ajustados <- function(modelo, ...) {
  if (inherits(modelo, c("estatR_logistica", "estatR_poisson", "estatR_nb", "estatR_zeroinfl"))) {
    return(stats::fitted(modelo$modelo, ...))
  }
  stats::fitted(modelo, ...)
}

#' @export
ajustados <- valores_ajustados

#' @export
valores_preditos <- valores_ajustados

#' @title Resumo do Modelo
#' @description Retorna ou exibe o resumo completo de um modelo estatístico.
#' @param modelo Objeto de modelo ajustado.
#' @param ... Argumentos adicionais repassados para \code{summary}.
#' @return Objeto de resumo estatístico.
#' @export
resumo <- function(modelo, ...) {
  if (inherits(modelo, c("estatR_logistica", "estatR_poisson", "estatR_nb", "estatR_zeroinfl"))) {
    return(summary(modelo$modelo, ...))
  }
  summary(modelo, ...)
}
