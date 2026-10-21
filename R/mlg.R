# ==============================================================================
# MÓDULO: MODELOS LINEARES GENERALIZADOS (MLG / GLM)
# ==============================================================================

#' @title Regressão por Modelos Lineares Generalizados (MLG / GLM)
#' @description Ajusta modelos lineares generalizados para diversas distribuições
#'   da família exponencial (Binomial/Logística, Poisson, Quasipoisson, Gama,
#'   Gaussiana, Binomial Negativa e Beta). Fornece resumo detalhado no console,
#'   tabela de coeficientes com Erro Padrão, intervalos de confiança de 95%,
#'   testes de hipótese (Wald e Razão de Verossimilhança / LRT), diagnóstico de ajuste,
#'   número de iterações do Escore de Fisher e suporte opcional a gráfico.
#'
#' @param dados Data frame contendo os dados (primeiro argumento para suporte a pipes).
#' @param formula Fórmula do modelo estatístico (ex: \code{y ~ x1 + x2}).
#' @param familia Família da distribuição. Aceita string (ex: \code{"poisson"},
#'   \code{"binomial"}, \code{"gama"}, \code{"quasipoisson"}, \code{"gaussian"},
#'   \code{"binomial_negativa"}, \code{"beta"}) ou função de família do R
#'   (ex: \code{poisson()}, \code{binomial(link = "probit")}).
#' @param ligacao Função de ligação (link). Se \code{NULL}, utiliza a ligação canônica da família.
#' @param teste Teste de hipótese para a significância dos termos:
#'   \code{"wald"} (padrão, rápido), \code{"lrt"} (Razão de Verossimilhança via Deviance),
#'   ou \code{"ambos"} (apresenta as duas estatísticas e p-valores lado a lado).
#' @param referencia Categoria de referência caso a variável resposta seja categórica (opcional).
#' @param grafico Lógico. Se \code{TRUE}, exibe o gráfico de diagnóstico ou efeito do modelo.
#'
#' @return Objeto ajustado da classe correspondente (\code{glm}, \code{negbin} ou \code{betareg})
#'   com atributo adicional de diagnóstico do estatR.
#' @export
regressao_mlg <- function(dados, formula,
                          familia = c("binomial", "poisson", "quasipoisson", "gama", "gaussian", "binomial_negativa", "beta"),
                          ligacao = NULL,
                          teste = c("wald", "lrt", "ambos"),
                          referencia = NULL,
                          grafico = FALSE) {

  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")
  if (!inherits(formula, "formula")) stop("O argumento 'formula' deve ser uma f\u00f3rmula v\u00e1lida (ex: y ~ x).")

  teste <- match.arg(teste)

  # Trata familia se for string
  nome_fam <- if (is.character(familia)) {
    match.arg(tolower(familia[1]), c("binomial", "poisson", "quasipoisson", "gama", "gamma", "gaussian", "gaussiana", "binomial_negativa", "beta"))
  } else if (inherits(familia, "family")) {
    familia$family
  } else {
    "custom"
  }

  y_name <- as.character(formula[[2]])
  if (!y_name %in% names(dados)) {
    stop(sprintf("Vari\u00e1vel resposta '%s' n\u00e3o encontrada nos dados.", y_name))
  }

  # Categoria de referência para respostas categóricas
  if (!is.null(referencia)) {
    if (!referencia %in% as.character(dados[[y_name]])) {
      stop(sprintf("A categoria de refer\u00eancia '%s' n\u00e3o existe na coluna '%s'.", referencia, y_name))
    }
    dados[[y_name]] <- stats::relevel(as.factor(dados[[y_name]]), ref = as.character(referencia))
  }

  # Despacho para famílias especializadas se ligacao for customizada ou específica
  if (nome_fam %in% c("binomial")) {
    link_use <- if (!is.null(ligacao)) ligacao else "logit"
    return(regressao_logistica(formula = formula, dados = dados, ligacao = link_use, grafico = grafico, referencia = referencia))
  }

  if (nome_fam %in% c("poisson")) {
    link_use <- if (!is.null(ligacao)) ligacao else "log"
    return(regressao_poisson(formula = formula, dados = dados, ligacao = link_use, grafico = grafico))
  }

  if (nome_fam %in% c("gama", "gamma")) {
    link_use <- if (!is.null(ligacao)) ligacao else "log"
    return(regressao_gama(formula = formula, dados = dados, ligacao = link_use, grafico = grafico))
  }

  if (nome_fam %in% c("binomial_negativa")) {
    return(regressao_binomial_negativa(formula = formula, dados = dados, grafico = grafico))
  }

  if (nome_fam %in% c("beta")) {
    return(regressao_beta(formula = formula, dados = dados, grafico = grafico))
  }

  # Caso padrão: Gaussian ou Quasipoisson ou GLM Geral
  obj_family <- if (is.character(familia)) {
    if (nome_fam %in% c("gaussian", "gaussiana")) {
      if (is.null(ligacao)) stats::gaussian() else stats::gaussian(link = ligacao)
    } else if (nome_fam == "quasipoisson") {
      if (is.null(ligacao)) stats::quasipoisson() else stats::quasipoisson(link = ligacao)
    } else {
      stats::gaussian()
    }
  } else {
    familia
  }

  mod <- stats::glm(formula, data = dados, family = obj_family)
  resumo <- summary(mod)
  coefs <- resumo$coefficients

  # Intervalos de Confiança 95%
  ci <- tryCatch(suppressMessages(stats::confint.default(mod)), error = function(e) {
    cbind(coefs[, 1] - 1.96 * coefs[, 2], coefs[, 1] + 1.96 * coefs[, 2])
  })

  # Cálculo do LRT se solicitado
  lrt_pvals <- NULL
  if (teste %in% c("lrt", "ambos")) {
    d1 <- tryCatch(stats::drop1(mod, test = "Chisq"), error = function(e) NULL)
    if (!is.null(d1) && "Pr(>Chi)" %in% names(d1)) {
      lrt_pvals <- d1[["Pr(>Chi)"]]
      names(lrt_pvals) <- rownames(d1)
    }
  }

  # Formatação e Exibição Console
  pad <- function(s, w, align = "left") {
    s <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right") return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }
  fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001")
    formatC(p, format = "f", digits = 3, decimal.mark = ",")
  }
  ast <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    if (p < 0.10) return("*")
    return("")
  }

  if (exists(".print_titulo", mode = "function")) {
    .print_titulo(sprintf("MODELO LINEAR GENERALIZADO (MLG) \u2014 FAM\u00cdLIA: %s", toupper(obj_family$family)))
  } else {
    cat(sprintf("\n── MODELO LINEAR GENERALIZADO (MLG) — FAMÍLIA: %s ──\n", toupper(obj_family$family)))
  }

  cat(sprintf("  Equa\u00e7\u00e3o:              %s\n", deparse(formula)))
  cat(sprintf("  Fun\u00e7\u00e3o de Liga\u00e7\u00e3o:       %s\n", obj_family$link))
  cat(sprintf("  N\u00ba de Observa\u00e7\u00f5es:     %d\n", nrow(mod$model)))
  cat(sprintf("  Itera\u00e7\u00f5es (Fisher):   %d\n\n", if (!is.null(mod$iter)) mod$iter else NA))

  if (exists(".print_topico", mode = "function")) .print_topico("AJUSTE DO MODELO")
  cat(sprintf("  Deviance Nula:       %.2f (df = %d)\n", mod$null.deviance, mod$df.null))
  cat(sprintf("  Deviance Residual:   %.2f (df = %d)\n", mod$deviance, mod$df.residual))
  cat(sprintf("  AIC:                 %.2f\n\n", stats::AIC(mod)))

  if (exists(".print_topico", mode = "function")) .print_topico("COEFICIENTES DO MODELO")

  rn <- rownames(coefs)
  rn[rn == "(Intercept)"] <- "(Intercepto)"

  w_c <- if (teste == "ambos") {
    c(var = 18, est = 12, ep = 12, stat = 10, p_wald = 15, p_lrt = 15, ic = 19)
  } else if (teste == "lrt") {
    c(var = 18, est = 12, ep = 12, p_lrt = 16, ic = 19)
  } else {
    c(var = 18, est = 12, ep = 12, stat = 10, p_wald = 16, ic = 19)
  }

  sep_line <- paste0("  ", strrep("\u2500", sum(w_c) + length(w_c) - 1))
  cat(sep_line, "\n  ")
  cat(pad("Vari\u00e1vel", w_c["var"], "left"))
  cat(pad("Estimativa", w_c["est"], "center"))
  cat(pad("Erro Pad.", w_c["ep"], "center"))
  if (teste %in% c("wald", "ambos")) cat(pad("Estat\u00edstica", w_c["stat"], "center"))
  if (teste %in% c("wald", "ambos")) cat(pad("p-valor (Wald)", w_c["p_wald"], "center"))
  if (teste %in% c("lrt", "ambos")) cat(pad("p-valor (LRT)", w_c["p_lrt"], "center"))
  cat(pad("IC (95%)", w_c["ic"], "center"), "\n")
  cat(sep_line, "\n")

  for (i in seq_len(nrow(coefs))) {
    v_nm <- rownames(coefs)[i]
    ic_txt <- sprintf("[%.2f; %.2f]", ci[i, 1], ci[i, 2])
    cat("  ")
    cat(pad(rn[i], w_c["var"], "left"))
    cat(pad(sprintf("%.3f", coefs[i, 1]), w_c["est"], "center"))
    cat(pad(sprintf("%.3f", coefs[i, 2]), w_c["ep"], "center"))
    if (teste %in% c("wald", "ambos")) cat(pad(sprintf("%.2f", coefs[i, 3]), w_c["stat"], "center"))
    if (teste %in% c("wald", "ambos")) {
      pw <- coefs[i, 4]
      cat(pad(paste(fmt_p(pw), ast(pw)), w_c["p_wald"], "center"))
    }
    if (teste %in% c("lrt", "ambos")) {
      pl <- if (!is.null(lrt_pvals) && v_nm %in% names(lrt_pvals)) lrt_pvals[[v_nm]] else NA
      cat(pad(paste(fmt_p(pl), ast(pl)), w_c["p_lrt"], "center"))
    }
    cat(pad(ic_txt, w_c["ic"], "center"), "\n")
  }
  cat(sep_line, "\n")
  cat("  Notas:\n")
  cat("  Signific\u00e2ncia: *** p < 0,01   ** p < 0,05   * p < 0,10\n\n")

  if (grafico && exists("grafico_coeficientes", mode = "function")) {
    tryCatch(print(grafico_coeficientes(mod)), error = function(e) NULL)
  }

  invisible(mod)
}

