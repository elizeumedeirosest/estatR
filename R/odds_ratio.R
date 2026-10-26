# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: EXTRAÇÃO E INTERPRETAÇÃO DE ODDS RATIO (OR) E IRR
# ─────────────────────────────────────────────────────────────────────────────

#' @title Odds Ratio (Razão de Chances) com Interpretação e Gráfico Forest Plot
#' @description Extrai ou calcula a Razão de Chances (Odds Ratio - OR) a partir de
#'   modelos de regressão logística (\code{regressao_logistica}, \code{glm} binomial/quasibinomial,
#'   \code{polr} ou \code{multinom}), fornecendo intervalos de confiança, percentuais de aumento/redução,
#'   interpretação textual automatizada e gráfico tipo Forest Plot.
#'
#' @param modelo Objeto retornado por \code{regressao_logistica()}, \code{glm(..., family = binomial)},
#'   \code{polr()} ou \code{multinom()}.
#' @param confianca Nível de confiança bilateral para os intervalos (padrão 0.95).
#' @param digitos Número de casas decimais para exibição (padrão 3).
#' @param grafico Lógico. Se \code{TRUE}, gera o gráfico Forest Plot dos Odds Ratios.
#' @param ... Argumentos adicionais.
#'
#' @return Objeto da classe \code{c("estatR_odds_ratio", "data.frame")} contendo
#'   os coeficientes, erros-padrão, valores z, p-valores, OR, intervalos de confiança e interpretação.
#' @export
#' @aliases razao_chances
odds_ratio <- function(modelo, confianca = 0.95, digitos = 3, grafico = FALSE, ...) {
  mod_obj <- if (inherits(modelo, "estatR_logistica")) modelo$modelo else modelo

  if (!inherits(mod_obj, c("glm", "polr", "multinom"))) {
    stop("O modelo deve ser um objeto de regressao_logistica(), glm(family = binomial), polr() ou multinom().", call. = FALSE)
  }

  alpha  <- 1 - confianca
  z_crit <- stats::qnorm(1 - alpha / 2)

  if (inherits(mod_obj, "glm")) {
    sum_mod <- summary(mod_obj)
    coefs   <- sum_mod$coefficients
    termos  <- rownames(coefs)

    idx_pred <- which(termos != "(Intercept)")
    if (length(idx_pred) == 0) idx_pred <- 1:nrow(coefs)

    termos_pred <- termos[idx_pred]
    est         <- coefs[idx_pred, 1]
    se          <- coefs[idx_pred, 2]
    stat        <- coefs[idx_pred, 3]
    pval        <- coefs[idx_pred, 4]

    or_val  <- exp(est)
    ic_inf  <- exp(est - z_crit * se)
    ic_sup  <- exp(est + z_crit * se)
    pct_val <- (or_val - 1) * 100

    sig_vec <- ifelse(pval < alpha, "Sim", "Não")
    interp_vec <- character(length(termos_pred))

    for (k in seq_along(termos_pred)) {
      nm <- termos_pred[k]
      p_fmt <- if (pval[k] < 0.001) "< 0,001" else sprintf("= %.*f", digitos, pval[k])
      if (sig_vec[k] == "Sim") {
        if (or_val[k] > 1) {
          interp_vec[k] <- sprintf(
            "A cada aumento de 1 unid. em '%s', as chances do evento aumentam em %.1f%% (OR = %.*f; IC %g%%: [%.*f; %.*f]; p %s).",
            nm, pct_val[k], digitos, or_val[k], confianca * 100, digitos, ic_inf[k], digitos, ic_sup[k], p_fmt
          )
        } else {
          interp_vec[k] <- sprintf(
            "A cada aumento de 1 unid. em '%s', as chances do evento reduzem em %.1f%% (OR = %.*f; IC %g%%: [%.*f; %.*f]; p %s).",
            nm, abs(pct_val[k]), digitos, or_val[k], confianca * 100, digitos, ic_inf[k], digitos, ic_sup[k], p_fmt
          )
        }
      } else {
        interp_vec[k] <- sprintf(
          "O preditor '%s' não apresentou associação estatisticamente significativa com as chances do desfecho (p %s; OR = %.*f; IC %g%%: [%.*f; %.*f]).",
          nm, p_fmt, digitos, or_val[k], confianca * 100, digitos, ic_inf[k], digitos, ic_sup[k]
        )
      }
    }

    df_res <- data.frame(
      Termo         = termos_pred,
      Estimativa    = est,
      Erro_Padrao   = se,
      Z_valor       = stat,
      p_valor       = pval,
      OR            = or_val,
      IC_inf        = ic_inf,
      IC_sup        = ic_sup,
      Variacao_Pct  = pct_val,
      Sig           = sig_vec,
      Interpretacao = interp_vec,
      stringsAsFactors = FALSE
    )
  } else if (inherits(mod_obj, "polr")) {
    sum_mod     <- summary(mod_obj)
    coefs       <- sum_mod$coefficients
    termos_pred <- names(mod_obj$coefficients)
    coefs_pred  <- coefs[termos_pred, , drop = FALSE]

    est  <- coefs_pred[, "Value"]
    se   <- coefs_pred[, "Std. Error"]
    stat <- coefs_pred[, "t value"]
    pval <- 2 * (1 - stats::pnorm(abs(stat)))

    or_val  <- exp(est)
    ic_inf  <- exp(est - z_crit * se)
    ic_sup  <- exp(est + z_crit * se)
    pct_val <- (or_val - 1) * 100
    sig_vec <- ifelse(pval < alpha, "Sim", "Não")
    interp_vec <- character(length(termos_pred))

    for (k in seq_along(termos_pred)) {
      nm <- termos_pred[k]
      p_fmt <- if (pval[k] < 0.001) "< 0,001" else sprintf("= %.*f", digitos, pval[k])
      if (sig_vec[k] == "Sim") {
        dir_txt <- if (or_val[k] > 1) "aumentam" else "reduzem"
        interp_vec[k] <- sprintf(
          "A cada aumento de 1 unid. em '%s', as chances de categorias superiores %s em %.1f%% (OR = %.*f; IC %g%%: [%.*f; %.*f]; p %s).",
          nm, dir_txt, abs(pct_val[k]), digitos, or_val[k], confianca * 100, digitos, ic_inf[k], digitos, ic_sup[k], p_fmt
        )
      } else {
        interp_vec[k] <- sprintf(
          "O preditor '%s' não apresentou associação estatisticamente significativa na escala ordinal (p %s).",
          nm, p_fmt
        )
      }
    }

    df_res <- data.frame(
      Termo         = termos_pred,
      Estimativa    = est,
      Erro_Padrao   = se,
      Z_valor       = stat,
      p_valor       = pval,
      OR            = or_val,
      IC_inf        = ic_inf,
      IC_sup        = ic_sup,
      Variacao_Pct  = pct_val,
      Sig           = sig_vec,
      Interpretacao = interp_vec,
      stringsAsFactors = FALSE
    )
  } else {
    sum_mod  <- summary(mod_obj)
    coefs    <- sum_mod$coefficients
    ses      <- sum_mod$standard.errors
    df_lista <- list()

    for (cat_j in rownames(coefs)) {
      est  <- coefs[cat_j, ]
      se   <- ses[cat_j, ]
      stat <- est / se
      pval <- 2 * (1 - stats::pnorm(abs(stat)))
      idx_p <- which(names(est) != "(Intercept)")
      if (length(idx_p) == 0) idx_p <- seq_along(est)

      or_val  <- exp(est[idx_p])
      ic_inf  <- exp(est[idx_p] - z_crit * se[idx_p])
      ic_sup  <- exp(est[idx_p] + z_crit * se[idx_p])
      pct_val <- (or_val - 1) * 100
      sig_vec <- ifelse(pval[idx_p] < alpha, "Sim", "Não")
      t_nms   <- paste0(cat_j, ": ", names(est)[idx_p])

      df_lista[[cat_j]] <- data.frame(
        Termo         = t_nms,
        Estimativa    = est[idx_p],
        Erro_Padrao   = se[idx_p],
        Z_valor       = stat[idx_p],
        p_valor       = pval[idx_p],
        OR            = or_val,
        IC_inf        = ic_inf,
        IC_sup        = ic_sup,
        Variacao_Pct  = pct_val,
        Sig           = sig_vec,
        Interpretacao = paste0("Categoria ", cat_j, ": OR = ", round(or_val, digitos)),
        stringsAsFactors = FALSE
      )
    }
    df_res <- do.call(rbind, df_lista)
    rownames(df_res) <- NULL
  }

  attr(df_res, "confianca") <- confianca
  attr(df_res, "digitos")   <- digitos
  class(df_res) <- c("estatR_odds_ratio", "data.frame")

  if (grafico) {
    p <- grafico_odds_ratio(df_res)
    print(p)
  }

  df_res
}

#' @export
razao_chances <- odds_ratio


#' @title Razão de Taxa de Incidência (IRR) com Interpretação e Forest Plot
#' @description Calcula a Razão de Taxas de Incidência (Incidence Rate Ratio - IRR = exp(beta))
#'   para modelos de contagem (\code{regressao_poisson}, \code{regressao_binomial_negativa},
#'   \code{modelo_zero_inflado}, \code{glm(family = poisson)}, \code{glm.nb} ou \code{zeroinfl}),
#'   com intervalos de confiança, variação percentual, interpretação textual e gráfico Forest Plot.
#'
#' @param modelo Modelo ajustado de contagem (\code{regressao_poisson}, \code{regressao_binomial_negativa},
#'   \code{modelo_zero_inflado}, \code{glm}, \code{glm.nb} ou \code{zeroinfl}).
#' @param confianca Nível de confiança bilateral (padrão 0.95).
#' @param digitos Número de casas decimais para exibição (padrão 3).
#' @param grafico Lógico. Se \code{TRUE}, gera o gráfico Forest Plot dos IRRs.
#' @param ... Argumentos adicionais.
#'
#' @return Objeto da classe \code{c("estatR_irr", "data.frame")} com estatísticas,
#'   IRR, intervalos de confiança e interpretação textual.
#' @export
#' @aliases razao_taxas
irr <- function(modelo, confianca = 0.95, digitos = 3, grafico = FALSE, ...) {
  mod_obj <- if (inherits(modelo, c("estatR_poisson", "estatR_nb", "estatR_zeroinfl"))) modelo$modelo else modelo

  if (!inherits(mod_obj, c("glm", "negbin", "zeroinfl"))) {
    stop("O modelo deve ser um objeto de regressao_poisson(), regressao_binomial_negativa(), modelo_zero_inflado() ou glm(family = poisson).", call. = FALSE)
  }

  alpha  <- 1 - confianca
  z_crit <- stats::qnorm(1 - alpha / 2)

  if (inherits(mod_obj, "zeroinfl")) {
    sum_mod <- summary(mod_obj)
    coefs   <- sum_mod$coefficients$count
  } else {
    sum_mod <- summary(mod_obj)
    coefs   <- sum_mod$coefficients
  }

  termos   <- rownames(coefs)
  idx_pred <- which(termos != "(Intercept)")
  if (length(idx_pred) == 0) idx_pred <- 1:nrow(coefs)

  termos_pred <- termos[idx_pred]
  est         <- coefs[idx_pred, 1]
  se          <- coefs[idx_pred, 2]
  stat        <- coefs[idx_pred, 3]
  pval        <- coefs[idx_pred, 4]

  irr_val <- exp(est)
  ic_inf  <- exp(est - z_crit * se)
  ic_sup  <- exp(est + z_crit * se)
  pct_val <- (irr_val - 1) * 100

  sig_vec <- ifelse(pval < alpha, "Sim", "Não")
  interp_vec <- character(length(termos_pred))

  for (k in seq_along(termos_pred)) {
    nm <- termos_pred[k]
    p_fmt <- if (pval[k] < 0.001) "< 0,001" else sprintf("= %.*f", digitos, pval[k])
    if (sig_vec[k] == "Sim") {
      if (irr_val[k] > 1) {
        interp_vec[k] <- sprintf(
          "A cada aumento de 1 unid. em '%s', a taxa esperada do evento aumenta em %.1f%% (IRR = %.*f; IC %g%%: [%.*f; %.*f]; p %s).",
          nm, pct_val[k], digitos, irr_val[k], confianca * 100, digitos, ic_inf[k], digitos, ic_sup[k], p_fmt
        )
      } else {
        interp_vec[k] <- sprintf(
          "A cada aumento de 1 unid. em '%s', a taxa esperada do evento reduz em %.1f%% (IRR = %.*f; IC %g%%: [%.*f; %.*f]; p %s).",
          nm, abs(pct_val[k]), digitos, irr_val[k], confianca * 100, digitos, ic_inf[k], digitos, ic_sup[k], p_fmt
        )
      }
    } else {
      interp_vec[k] <- sprintf(
        "O preditor '%s' não apresentou associação estatisticamente significativa com a contagem esperada (p %s; IRR = %.*f; IC %g%%: [%.*f; %.*f]).",
        nm, p_fmt, digitos, irr_val[k], confianca * 100, digitos, ic_inf[k], digitos, ic_sup[k]
      )
    }
  }

  df_res <- data.frame(
    Termo         = termos_pred,
    Estimativa    = est,
    Erro_Padrao   = se,
    Z_valor       = stat,
    p_valor       = pval,
    IRR           = irr_val,
    IC_inf        = ic_inf,
    IC_sup        = ic_sup,
    Variacao_Pct  = pct_val,
    Sig           = sig_vec,
    Interpretacao = interp_vec,
    stringsAsFactors = FALSE
  )

  attr(df_res, "confianca") <- confianca
  attr(df_res, "digitos")   <- digitos
  class(df_res) <- c("estatR_irr", "data.frame")

  if (grafico) {
    p <- grafico_irr(df_res)
    print(p)
  }

  df_res
}

#' @export
razao_taxas <- irr


# ─────────────────────────────────────────────────────────────────────────────
# MÉTODOS DE IMPRESSÃO
# ─────────────────────────────────────────────────────────────────────────────

#' @export
print.estatR_odds_ratio <- function(x, ...) {
  digitos <- attr(x, "digitos")
  if (is.null(digitos)) digitos <- 3
  conf <- attr(x, "confianca")
  if (is.null(conf)) conf <- 0.95

  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("ODDS RATIO (RAZÃO DE CHANCES) — REGRESSÃO LOGÍSTICA")
  } else {
    cat("\n=== ODDS RATIO (RAZÃO DE CHANCES) ===\n")
  }
  cat(sprintf("  Nível de confiança dos intervalos: %g%%\n\n", conf * 100))

  if (exists(".print_topico", mode = "function")) {
    .print_topico("TABELA DE ODDS RATIO")
  } else {
    cat("--- TABELA DE ODDS RATIO ---\n")
  }
  cat("  ─────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("  %-20s  %7s  %7s  %9s  %9s  %16s  %4s\n",
              "Preditor", "Beta", "EP", "p-valor", "OR", sprintf("IC %g%%", conf * 100), "Sig"))
  cat("  ─────────────────────────────────────────────────────────────────────────────\n")
  for (i in 1:nrow(x)) {
    p_txt  <- if (x$p_valor[i] < 0.001) "< 0.001" else sprintf("%.4f", x$p_valor[i])
    ic_str <- sprintf(paste0("[%.*f; %.*f]"), digitos, x$IC_inf[i], digitos, x$IC_sup[i])
    cat(sprintf("  %-20s  %7.*f  %7.*f  %9s  %9.*f  %16s  %4s\n",
                x$Termo[i], digitos, x$Estimativa[i], digitos, x$Erro_Padrao[i],
                p_txt, digitos, x$OR[i], ic_str, x$Sig[i]))
  }
  cat("  ─────────────────────────────────────────────────────────────────────────────\n")
  cat("  * OR > 1: aumenta as chances | OR < 1: reduz as chances | OR = 1: efeito nulo.\n\n")

  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO AUTOMATIZADA")
  } else {
    cat("--- INTERPRETAÇÃO AUTOMATIZADA ---\n")
  }
  for (i in 1:nrow(x)) {
    cat(sprintf("  • %s\n", x$Interpretacao[i]))
  }
  cat("\n")
  invisible(x)
}

#' @export
print.estatR_irr <- function(x, ...) {
  digitos <- attr(x, "digitos")
  if (is.null(digitos)) digitos <- 3
  conf <- attr(x, "confianca")
  if (is.null(conf)) conf <- 0.95

  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("RAZÃO DE TAXAS DE INCIDÊNCIA (IRR) — MODELOS DE CONTAGEM")
  } else {
    cat("\n=== RAZÃO DE TAXAS DE INCIDÊNCIA (IRR) ===\n")
  }
  cat(sprintf("  Nível de confiança dos intervalos: %g%%\n\n", conf * 100))

  if (exists(".print_topico", mode = "function")) {
    .print_topico("TABELA DE INCIDENCE RATE RATIO (IRR = exp(beta))")
  } else {
    cat("--- TABELA DE INCIDENCE RATE RATIO ---\n")
  }
  cat("  ─────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("  %-20s  %7s  %7s  %9s  %9s  %16s  %4s\n",
              "Preditor", "Beta", "EP", "p-valor", "IRR", sprintf("IC %g%%", conf * 100), "Sig"))
  cat("  ─────────────────────────────────────────────────────────────────────────────\n")
  for (i in 1:nrow(x)) {
    p_txt  <- if (x$p_valor[i] < 0.001) "< 0.001" else sprintf("%.4f", x$p_valor[i])
    ic_str <- sprintf(paste0("[%.*f; %.*f]"), digitos, x$IC_inf[i], digitos, x$IC_sup[i])
    cat(sprintf("  %-20s  %7.*f  %7.*f  %9s  %9.*f  %16s  %4s\n",
                x$Termo[i], digitos, x$Estimativa[i], digitos, x$Erro_Padrao[i],
                p_txt, digitos, x$IRR[i], ic_str, x$Sig[i]))
  }
  cat("  ─────────────────────────────────────────────────────────────────────────────\n")
  cat("  * IRR > 1: aumenta a taxa | IRR < 1: reduz a taxa | IRR = 1: efeito nulo.\n\n")

  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO AUTOMATIZADA")
  } else {
    cat("--- INTERPRETAÇÃO AUTOMATIZADA ---\n")
  }
  for (i in 1:nrow(x)) {
    cat(sprintf("  • %s\n", x$Interpretacao[i]))
  }
  cat("\n")
  invisible(x)
}


# ─────────────────────────────────────────────────────────────────────────────
# MÉTODOS GRÁFICOS (FOREST PLOTS)
# ─────────────────────────────────────────────────────────────────────────────

#' @title Gráfico Forest Plot de Odds Ratio
#' @description Gera gráfico Forest Plot para os Odds Ratios estimados com intervalos de confiança.
#' @param x Objeto retornado por \code{odds_ratio()}.
#' @param cor_sig Cor para termos significativos (padrão \code{"#0099FF"}).
#' @param cor_ns Cor para termos não significativos (padrão \code{"gray50"}).
#' @return Objeto \code{ggplot}.
#' @export
grafico_odds_ratio <- function(x, cor_sig = "#0099FF", cor_ns = "gray50") {
  if (!inherits(x, "estatR_odds_ratio")) {
    stop("O argumento 'x' deve ser um objeto retornado por odds_ratio().", call. = FALSE)
  }

  df_plot <- as.data.frame(x)
  df_plot$Termo <- factor(df_plot$Termo, levels = rev(df_plot$Termo))
  conf <- attr(x, "confianca")
  if (is.null(conf)) conf <- 0.95

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = OR, y = Termo, color = Sig)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = IC_inf, xmax = IC_sup), height = 0.25, linewidth = 0.85) +
    ggplot2::geom_point(size = 3.2) +
    ggplot2::scale_color_manual(
      values = c("Não" = cor_ns, "Sim" = cor_sig),
      labels = c("Não" = "Não significativo (cruza 1)", "Sim" = "Significativo (p < 0.05)")
    ) +
    ggplot2::labs(
      title = "Odds Ratio (Razão de Chances) — Forest Plot",
      subtitle = sprintf("Linha tracejada vermelha em 1 indica ausência de efeito (IC %g%%)", conf * 100),
      x = "Odds Ratio (OR)",
      y = "Preditor",
      color = "Efeito"
    )

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p <- p + ggplot2::theme(legend.position = "bottom")
  p
}

#' @export
plot.estatR_odds_ratio <- function(x, ...) {
  print(grafico_odds_ratio(x, ...))
}

#' @title Gráfico Forest Plot de Razão de Taxa de Incidência (IRR)
#' @description Gera gráfico Forest Plot para os IRRs estimados com intervalos de confiança.
#' @param x Objeto retornado por \code{irr()}.
#' @param cor_sig Cor para termos significativos (padrão \code{"#0099FF"}).
#' @param cor_ns Cor para termos não significativos (padrão \code{"gray50"}).
#' @return Objeto \code{ggplot}.
#' @export
grafico_irr <- function(x, cor_sig = "#0099FF", cor_ns = "gray50") {
  if (!inherits(x, "estatR_irr")) {
    stop("O argumento 'x' deve ser um objeto retornado por irr().", call. = FALSE)
  }

  df_plot <- as.data.frame(x)
  df_plot$Termo <- factor(df_plot$Termo, levels = rev(df_plot$Termo))
  conf <- attr(x, "confianca")
  if (is.null(conf)) conf <- 0.95

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = IRR, y = Termo, color = Sig)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = IC_inf, xmax = IC_sup), height = 0.25, linewidth = 0.85) +
    ggplot2::geom_point(size = 3.2) +
    ggplot2::scale_color_manual(
      values = c("Não" = cor_ns, "Sim" = cor_sig),
      labels = c("Não" = "Não significativo (cruza 1)", "Sim" = "Significativo (p < 0.05)")
    ) +
    ggplot2::labs(
      title = "Razão de Taxas de Incidência (IRR) — Forest Plot",
      subtitle = sprintf("Linha tracejada vermelha em 1 indica ausência de efeito (IC %g%%)", conf * 100),
      x = "Incidence Rate Ratio (IRR)",
      y = "Preditor",
      color = "Efeito"
    )

  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p <- p + ggplot2::theme(legend.position = "bottom")
  p
}

#' @export
plot.estatR_irr <- function(x, ...) {
  print(grafico_irr(x, ...))
}
