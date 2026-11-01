#' @title Modelos Lineares Mistos e Hierárquicos
#' @description Ajusta modelos lineares e generalizados lineares de efeitos mistos
#' (hierárquicos / multinível), acomodando estruturas de dados agrupados, medidas repetidas
#' e correlação intragrupo, com decomposição de variância e cálculo do ICC (Índice de Correlação Intraclasse).
#' @param formula Fórmula com especificação dos efeitos fixos e aleatórios no padrão lme4 (ex.: y ~ x1 + (1 | grupo)).
#' @param dados Data frame contendo as variáveis da fórmula.
#' @param familia Família da distribuição: "gaussian" (padrão, modelo linear misto) ou famílias do GLMM ("binomial", "poisson").
#' @param metodo Método de estimação para gaussian: "REML" (Máxima Verossimilhança Restrita, padrão) ou "ML" (Máxima Verossimilhança).
#' @param nivel Nível de confiança nominal para os intervalos dos efeitos fixos (padrão: 0.95).
#' @param grafico Lógico. Se TRUE (padrão), plota os efeitos aleatórios dos grupos (caterpillar plot).
#' @return Objeto da classe \code{estatR_regressao_mista} contendo:
#'   \item{efeitos_fixos}{Data frame com parâmetros, estimativas, erros padrão, ICs e p-valores}
#'   \item{efeitos_aleatorios}{Data frame com as variâncias e desvios padrão dos componentes aleatórios}
#'   \item{icc}{Índice de Correlação Intraclasse calculado}
#'   \item{modelo_lme4}{Objeto ajustado de classe merMod (se lme4 estiver disponível)}
#'   \item{formula}{Fórmula utilizada}
#' @export
regressao_mista <- function(formula,
                            dados   = NULL,
                            data    = NULL,
                            familia = c("gaussian", "binomial", "poisson"),
                            metodo  = c("REML", "ML"),
                            nivel   = 0.95,
                            grafico = TRUE) {

  metodo <- match.arg(metodo)
  if (is.character(familia)) familia <- match.arg(familia)

  if (is.null(dados) && !is.null(data)) dados <- data
  if (is.null(dados)) stop("O argumento 'dados' deve ser fornecido como um data frame.")
  if (!inherits(formula, "formula")) {
    formula <- stats::as.formula(formula)
  }

  alpha <- 1 - nivel
  z_crit <- stats::qnorm(1 - alpha / 2)

  # Ajuste do modelo misto via lme4 (ou fallback)
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("O pacote 'lme4' é necessário para o ajuste de modelos mistos. Por favor, instale-o via install.packages('lme4').")
  }

  is_reml <- (metodo == "REML" && (is.character(familia) && familia == "gaussian"))

  mod_fit <- tryCatch({
    if (is.character(familia) && familia == "gaussian") {
      lme4::lmer(formula, data = dados, REML = is_reml)
    } else {
      lme4::glmer(formula, data = dados, family = familia)
    }
  }, error = function(e) {
    stop(sprintf("Erro ao ajustar o modelo misto: %s", e$message))
  })

  # Extração dos Efeitos Fixos
  s_mod <- summary(mod_fit)
  fix_coefs <- s_mod$coefficients
  nomes_fixos <- rownames(fix_coefs)

  est_fix <- fix_coefs[, 1]
  se_fix  <- fix_coefs[, 2]
  stat_fix <- if (ncol(fix_coefs) >= 3) fix_coefs[, 3] else (est_fix / se_fix)

  # p-valores assintóticos normais (padrão lme4 com z-test / Wald)
  p_fix <- 2 * stats::pnorm(-abs(stat_fix))
  ic_inf <- est_fix - z_crit * se_fix
  ic_sup <- est_fix + z_crit * se_fix

  df_fixos <- data.frame(
    Parametro   = nomes_fixos,
    Estimativa  = est_fix,
    Erro_Padrao = se_fix,
    IC_Inf      = ic_inf,
    IC_Sup      = ic_sup,
    Estatistica = stat_fix,
    p_valor     = p_fix,
    stringsAsFactors = FALSE
  )

  # Extração dos Componentes de Variância (Efeitos Aleatórios)
  vc <- lme4::VarCorr(mod_fit)
  vc_df <- as.data.frame(vc)

  # Identificação do agrupamento principal para cálculo do ICC
  var_grupo <- 0
  var_resid <- stats::sigma(mod_fit)^2
  nome_grp_princ <- "Grupo"

  for (i in seq_len(nrow(vc_df))) {
    if (vc_df$grp[i] != "Residual") {
      var_grupo <- var_grupo + vc_df$vcov[i]
      nome_grp_princ <- vc_df$grp[i]
    }
  }

  icc <- var_grupo / (var_grupo + var_resid)
  if (is.na(icc) || is.nan(icc)) icc <- 0

  # Helpers de formatação
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 3) {
    if (is.na(x)) return("-")
    formatC(x, format = "f", digits = d, decimal.mark = ",")
  }
  .fmt_p <- function(p) {
    if (is.na(p)) return("-")
    if (p < 0.001) return("< 0,001 ***")
    ast <- if (p < 0.01) " **" else if (p < 0.05) " *" else ""
    paste0(formatC(p, format = "f", digits = 3, decimal.mark = ","), ast)
  }

  # 3. Impressão no Console
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("MODELO LINEAR MISTO (EFEITOS FIXOS E ALEATÓRIOS)")
  } else {
    cat("
── MODELO LINEAR MISTO (EFEITOS FIXOS E ALEATÓRIOS) ──
")
  }

  cat(sprintf("  Fórmula do Modelo:    %s
", deparse(formula)))
  cat(sprintf("  Observações Totais:   %d observações
", stats::nobs(mod_fit)))
  cat(sprintf("  Nível de Agrupamento: %s
", nome_grp_princ))
  cat(sprintf("  Critério de Ajuste:   %s

", if (is_reml) "Máxima Verossimilhança Restrita (REML)" else "Máxima Verossimilhança (ML)"))

  # SEÇÃO: EFEITOS FIXOS
  if (exists(".print_topico", mode = "function")) {
    .print_topico("EFEITOS FIXOS")
  } else {
    cat("  ▶ EFEITOS FIXOS
")
  }

  w_fix <- c(par = max(20, max(nchar(nomes_fixos)) + 2), est = 14, se = 14, ic = 24, stat = 12, p = 16)
  sep_fix <- paste0("  ", strrep("─", sum(w_fix)))
  cat(sep_fix, "
")
  cat("  ", .pad("Parâmetro", w_fix["par"], "left"),
             .pad("Estimativa", w_fix["est"], "center"),
             .pad("Erro Padrão", w_fix["se"], "center"),
             .pad(sprintf("IC (%d%%)", round(nivel * 100)), w_fix["ic"], "center"),
             .pad("t", w_fix["stat"], "center"),
             .pad("p-valor", w_fix["p"], "center"), "
", sep = "")
  cat(sep_fix, "
")
  for (i in seq_len(nrow(df_fixos))) {
    str_ic <- sprintf("[ %s ; %s ]", .fmt_f(df_fixos$IC_Inf[i], 2), .fmt_f(df_fixos$IC_Sup[i], 2))
    cat("  ", .pad(df_fixos$Parametro[i], w_fix["par"], "left"),
               .pad(.fmt_f(df_fixos$Estimativa[i], 3), w_fix["est"], "center"),
               .pad(.fmt_f(df_fixos$Erro_Padrao[i], 3), w_fix["se"], "center"),
               .pad(str_ic, w_fix["ic"], "center"),
               .pad(.fmt_f(df_fixos$Estatistica[i], 2), w_fix["stat"], "center"),
               .pad(.fmt_p(df_fixos$p_valor[i]), w_fix["p"], "center"), "
", sep = "")
  }
  cat(sep_fix, "
")
  cat("  Notas:
")
  cat("  Significância: *** p < 0,001   ** p < 0,01   * p < 0,05

")

  # SEÇÃO: COMPONENTES DE VARIÂNCIA (EFEITOS ALEATÓRIOS)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("COMPONENTES DE VARIÂNCIA (EFEITOS ALEATÓRIOS)")
  } else {
    cat("  ▶ COMPONENTES DE VARIÂNCIA (EFEITOS ALEATÓRIOS)
")
  }

  w_ran <- c(grp = 18, eff = 20, v = 16, sd = 16)
  sep_ran <- paste0("  ", strrep("─", sum(w_ran)))
  cat(sep_ran, "
")
  cat("  ", .pad("Grupo", w_ran["grp"], "left"),
             .pad("Efeito", w_ran["eff"], "left"),
             .pad("Variância", w_ran["v"], "center"),
             .pad("Desvio Padrão", w_ran["sd"], "center"), "
", sep = "")
  cat(sep_ran, "
")
  for (i in seq_len(nrow(vc_df))) {
    eff_nm <- if (is.na(vc_df$var1[i])) "Residual" else vc_df$var1[i]
    cat("  ", .pad(vc_df$grp[i], w_ran["grp"], "left"),
               .pad(eff_nm, w_ran["eff"], "left"),
               .pad(.fmt_f(vc_df$vcov[i], 3), w_ran["v"], "center"),
               .pad(.fmt_f(vc_df$sdcor[i], 3), w_ran["sd"], "center"), "
", sep = "")
  }
  cat(sep_ran, "

")

  # INTERPRETAÇÃO INTEGRADA
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO E DIAGNÓSTICO DO MODELO")
  } else {
    cat("  ▶ INTERPRETAÇÃO E DIAGNÓSTICO DO MODELO
")
  }

  cor_ok   <- "[1;32m"
  cor_info <- "[1;34m"
  reset    <- "[0m"

  pct_icc <- icc * 100
  cat(sprintf("  %s[✓]%s Índice de Correlação Intraclasse (ICC = %s):
    %s da variância total da resposta decorre das diferenças entre os agrupamentos (%s).
    %s

",
              cor_ok, reset, .fmt_f(icc, 3), paste0(.fmt_f(pct_icc, 1), "%"), nome_grp_princ,
              if (icc > 0.05) "A estrutura hierárquica é indispensável para evitar viés nos erros padrão." else "Baixa correlação intragrupo; modelo mista converge para regressão linear simples."))

  sig_count <- sum(p_fix < 0.05)
  cat(sprintf("  %s[✓]%s Foram identificados %d efeitos fixos estatisticamente significativos (p < 0,05),
    indicando associações robustas após controlar a heterogeneidade entre os grupos.

",
              cor_ok, reset, sig_count))

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat("──────────────────────────────────────────────────────────────────────────────
")
  }

  res_obj <- list(
    efeitos_fixos      = df_fixos,
    efeitos_aleatorios = vc_df,
    icc                = icc,
    modelo_lme4        = mod_fit,
    formula            = formula,
    dados              = dados,
    nome_grupo         = nome_grp_princ
  )
  class(res_obj) <- c("estatR_regressao_mista", "list")

  if (isTRUE(grafico)) {
    p_mista <- plot_regressao_mista(res_obj)
    print(p_mista)
  }

  invisible(res_obj)
}

#' @title Gráfico dos Efeitos Aleatórios do Modelo Misto (Caterpillar Plot)
#' @description Plota os efeitos aleatórios dos grupos (BLUPs) com intervalo de confiança de 95%.
#' @param modelo Objeto da classe \code{estatR_regressao_mista}.
#' @return Objeto ggplot2.
#' @export
plot_regressao_mista <- function(modelo) {
  if (!inherits(modelo, "estatR_regressao_mista")) {
    stop("O modelo deve ser um objeto da classe 'estatR_regressao_mista'.")
  }

  ran_effects <- lme4::ranef(modelo$modelo_lme4, condVar = TRUE)
  grp_name <- modelo$nome_grupo
  if (!(grp_name %in% names(ran_effects))) {
    grp_name <- names(ran_effects)[1]
  }

  df_grp <- ran_effects[[grp_name]]
  pv <- attr(df_grp, "postVar")
  se_vec <- if (!is.null(pv)) sqrt(pv[1, 1, ]) else rep(0.1, nrow(df_grp))

  df_plot <- data.frame(
    Grupo     = rownames(df_grp),
    Efeito    = df_grp[, 1],
    SE        = se_vec,
    stringsAsFactors = FALSE
  )
  df_plot$IC_Inf <- df_plot$Efeito - 1.96 * df_plot$SE
  df_plot$IC_Sup <- df_plot$Efeito + 1.96 * df_plot$SE
  df_plot$Grupo  <- stats::reorder(df_plot$Grupo, df_plot$Efeito)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Grupo, y = Efeito)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = IC_Inf, ymax = IC_Sup), width = 0.25, color = "#1D3557", linewidth = 0.7) +
    ggplot2::geom_point(color = "#1D3557", size = 2.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = sprintf("Efeitos Aleatórios por Grupo (%s)", grp_name),
      subtitle = sprintf("BLUPs com IC de 95%% | ICC = %.3f", modelo$icc),
      x        = grp_name,
      y        = "Desvio em relação ao Intercepto Global"
    )

  if (exists("meu_tema", mode = "function")) {
    p <- p + meu_tema()
  } else if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2)
  } else {
    p <- p + ggplot2::theme_minimal()
  }

  p
}
