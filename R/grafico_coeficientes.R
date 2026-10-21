# ==============================================================================
# MÓDULO: VISUALIZAÇÃO DE COEFICIENTES (FOREST PLOT / ESTILO SJPLOT)
# ==============================================================================

#' @title Gráfico de Coeficientes do Modelo (Forest Plot)
#' @description Gera um gráfico elegante de coeficientes com intervalos de confiança
#'   de 95% (Forest Plot), nativo em \code{ggplot2}, com grade dupla por padrão
#'   e integração ao \code{tema_estatR()}.
#'   Suporta modelos lineares (\code{lm}), generalizados (\code{glm}),
#'   sobrevivência (\code{coxph}, \code{estatR_cox}) e outros.
#'
#' @param modelo Objeto de modelo estatístico ajustado (\code{lm}, \code{glm}, \code{coxph}, etc.).
#' @param incluir_intercepto Lógico. Se \code{FALSE} (padrão), remove o intercepto para focar nos efeitos dos preditores.
#' @param exponencial Lógico ou \code{"auto"}. Se \code{TRUE}, plota os coeficientes exponenciados
#'   (Odds Ratio para Logística, Hazard Ratio para Cox, IRR para Poisson) e centraliza a linha nula em 1.
#'   Se \code{"auto"}, detecta automaticamente a família do modelo.
#' @param cor Cor dos pontos e intervalos de confiança (padrão: \code{"#1E3A8A"}).
#' @param tamanho_ponto Tamanho do ponto da estimativa (padrão: 3).
#' @param tema_estilo Estilo do \code{tema_estatR()} a ser aplicado (padrão: 2 com grade dupla).
#'
#' @return Um objeto \code{ggplot}.
#' @export
grafico_coeficientes <- function(modelo,
                                 incluir_intercepto = FALSE,
                                 exponencial = "auto",
                                 cor = "#1E3A8A",
                                 tamanho_ponto = 3,
                                 tema_estilo = 2) {

  if (!inherits(modelo, c("lm", "glm", "coxph", "estatR_cox", "betareg", "negbin"))) {
    stop("O modelo deve ser um objeto ajustado de regress\u00e3o (lm, glm, coxph, etc.).")
  }

  # Se for estatR_cox, extrai o modelo coxph interno
  if (inherits(modelo, "estatR_cox")) {
    modelo <- modelo$modelo
  }

  sm <- summary(modelo)
  coefs <- sm$coefficients

  # Determina se deve exponenciar
  deve_exp <- FALSE
  if (identical(exponencial, "auto")) {
    if (inherits(modelo, "coxph")) {
      deve_exp <- TRUE
    } else if (inherits(modelo, "glm") && !is.null(modelo$family)) {
      fam <- modelo$family$family
      link <- modelo$family$link
      if (fam %in% c("binomial", "quasibinomial") && link == "logit") deve_exp <- TRUE
      if (fam %in% c("poisson", "quasipoisson") && link == "log") deve_exp <- TRUE
      if (grepl("Negative Binomial", fam, ignore.case = TRUE)) deve_exp <- TRUE
    }
  } else {
    deve_exp <- isTRUE(exponencial)
  }

  # Coleta estimativas e ICs
  ci <- tryCatch(suppressMessages(stats::confint(modelo)), error = function(e) {
    # Fallback para confint.default (Wald)
    suppressMessages(stats::confint.default(modelo))
  })

  rn <- rownames(coefs)

  # Remove intercepto se solicitado
  if (!incluir_intercepto) {
    idx_int <- grep("^(Intercept|\\(Intercept\\)|Intercepto)", rn, ignore.case = TRUE)
    if (length(idx_int) > 0) {
      coefs <- coefs[-idx_int, , drop = FALSE]
      ci    <- ci[-idx_int, , drop = FALSE]
      rn    <- rn[-idx_int]
    }
  }

  if (nrow(coefs) == 0) {
    stop("Nenhum coeficiente restante para plotar ap\u00f3s remover o intercepto.")
  }

  # Extrai valores pontuais
  estimativas <- coefs[, 1]
  ic_inf <- ci[, 1]
  ic_sup <- ci[, 2]

  rotulo_eixo_x <- "Estimativa (\u03b2)"
  linha_nula <- 0

  if (deve_exp) {
    estimativas <- exp(estimativas)
    ic_inf      <- exp(ic_inf)
    ic_sup      <- exp(ic_sup)
    linha_nula  <- 1

    if (inherits(modelo, "coxph")) {
      rotulo_eixo_x <- "Raz\u00e3o de Riscos (Hazard Ratio - HR)"
    } else if (inherits(modelo, "glm") && modelo$family$family %in% c("binomial", "quasibinomial")) {
      rotulo_eixo_x <- "Raz\u00e3o de Chances (Odds Ratio - OR)"
    } else if (inherits(modelo, "glm") && modelo$family$family %in% c("poisson", "quasipoisson")) {
      rotulo_eixo_x <- "Raz\u00e3o de Taxas de Incid\u00eancia (IRR)"
    } else {
      rotulo_eixo_x <- "Estimativa Exponenciada (exp(\u03b2))"
    }
  }

  # Monta data.frame
  df_plot <- data.frame(
    Termo      = factor(rn, levels = rev(rn)),
    Estimativa = estimativas,
    IC_Inf     = ic_inf,
    IC_Sup     = ic_sup,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Estimativa, y = Termo)) +
    ggplot2::geom_vline(xintercept = linha_nula, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = IC_Inf, xmax = IC_Sup), height = 0.2, color = cor, linewidth = 0.8) +
    ggplot2::geom_point(size = tamanho_ponto, color = cor, fill = "white", shape = 21, stroke = 1.2) +
    ggplot2::labs(
      title = "Estimativas dos Coeficientes e Intervalos de Confian\u00e7a (95%)",
      x = rotulo_eixo_x,
      y = NULL
    )

  # Aplica tema estatR com grade dupla ativada
  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = tema_estilo, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal() + ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "gray90"))
  }

  p
}

