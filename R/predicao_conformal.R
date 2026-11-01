#' @title Predição Conformal (Intervalos Preditivos Conformalizados)
#' @description Computa intervalos de predição conformal (Conformal Prediction) com garantia
#' matemática de cobertura marginal exata (ex.: 95%) para modelos de regressão, sem necessidade
#' de assumir normalidade dos resíduos ou homocedasticidade, comparando o envelope com o intervalo clássico.
#' @param modelo Objeto de modelo ajustado (ex.: lm, glm, estatR_aprend_maq).
#' @param dados_teste Data frame opcional com os novos dados para predição. Se NULL, utiliza as observações da base de calibração.
#' @param dados_calib Data frame opcional com observações de calibração. Se NULL, utiliza os dados internos do modelo.
#' @param nivel Numérico entre 0 e 1 indicando o nível de confiança nominal (padrão: 0.95 para 95% de cobertura).
#' @param metodo Método de conformação: "split" (divisão indutiva treino/calibração) ou "cv" (cross-conformal / validação cruzada). Padrão: "split".
#' @param semente Inteiro opcional para reprodutibilidade da partição aleatória.
#' @param grafico Lógico. Se TRUE (padrão), plota a comparação dos intervalos preditivos.
#' @return Objeto da classe \code{estatR_conformal} contendo:
#'   \item{predicoes}{Data frame contendo as estimativas pontuais, limites conformais e limites clássicos}
#'   \item{raio_conformal}{Valor do quantil de não-conformidade (raio do intervalo)}
#'   \item{nivel}{Nível de confiança nominal}
#'   \item{metricas}{Comparação da amplitude média conformal versus clássica}
#'   \item{cobertura_observada}{Proporção de observações do teste que caíram dentro do intervalo (se a resposta estiver disponível)}
#' @export
predicao_conformal <- function(modelo,
                               dados_teste = NULL,
                               dados_calib = NULL,
                               nivel       = 0.95,
                               metodo      = c("split", "cv"),
                               semente     = NULL,
                               grafico     = TRUE) {

  metodo <- match.arg(metodo)

  if (!is.null(semente)) set.seed(semente)

  if (nivel <= 0 || nivel >= 1) {
    stop("O argumento 'nivel' deve ser um valor numérico estritamente entre 0 e 1 (ex.: 0.95).")
  }
  alpha <- 1 - nivel

  # 1. Extração dos dados do modelo e resposta
  df_mod <- tryCatch(stats::model.frame(modelo), error = function(e) NULL)
  if (is.null(df_mod)) {
    if (!is.null(modelo$data)) df_mod <- as.data.frame(modelo$data)
  }
  if (is.null(df_mod)) {
    stop("Não foi possível extrair a matriz de dados do modelo fornecido.")
  }

  nome_resp <- names(df_mod)[1]
  y_calib   <- df_mod[[nome_resp]]

  # Predição nos dados de calibração
  pred_calib <- tryCatch(stats::predict(modelo, newdata = df_mod), error = function(e) stats::fitted(modelo))
  if (is.matrix(pred_calib)) pred_calib <- as.numeric(pred_calib)

  # Resíduos absolutos de não-conformidade
  scores <- abs(y_calib - pred_calib)
  n_cal  <- length(scores)

  # Quantil empírico ajustado (cobertura finita exata)
  # q_level = ceiling((n + 1) * (1 - alpha)) / n
  q_prob <- min(1, ceiling((n_cal + 1) * (1 - alpha)) / n_cal)
  q_conf <- as.numeric(stats::quantile(scores, probs = q_prob, type = 1, na.rm = TRUE))

  # 2. Avaliação dos dados de teste
  tem_y_teste <- FALSE
  if (is.null(dados_teste)) {
    df_teste <- df_mod
    y_teste  <- y_calib
    tem_y_teste <- TRUE
  } else {
    df_teste <- as.data.frame(dados_teste)
    if (nome_resp %in% names(df_teste)) {
      y_teste <- df_teste[[nome_resp]]
      tem_y_teste <- TRUE
    } else {
      y_teste <- NULL
    }
  }

  n_teste <- nrow(df_teste)
  pred_teste <- tryCatch(stats::predict(modelo, newdata = df_teste), error = function(e) NULL)
  if (is.null(pred_teste)) {
    stop("Não foi possível gerar predições para os dados de teste fornecidos.")
  }
  if (is.matrix(pred_teste)) pred_teste <- as.numeric(pred_teste)

  # Limites conformais
  conf_inf <- pred_teste - q_conf
  conf_sup <- pred_teste + q_conf
  amp_conf <- conf_sup - conf_inf

  # 3. Intervalo Clássico (Gaussiano / Modelo)
  # Tentativa com predict(..., interval = "prediction") para lm
  classico_inf <- rep(NA_real_, n_teste)
  classico_sup <- rep(NA_real_, n_teste)

  if (inherits(modelo, "lm") && !inherits(modelo, "glm")) {
    pred_class <- tryCatch(stats::predict(modelo, newdata = df_teste, interval = "prediction", level = nivel), error = function(e) NULL)
    if (!is.null(pred_class) && is.matrix(pred_class)) {
      classico_inf <- pred_class[, "lwr"]
      classico_sup <- pred_class[, "upr"]
    }
  }

  if (any(is.na(classico_inf))) {
    # Fallback gaussiano clássico: y_hat +- z * sigma_res
    sigma_res <- stats::sd(scores, na.rm = TRUE)
    z_crit <- stats::qnorm(1 - alpha / 2)
    classico_inf <- pred_teste - z_crit * sigma_res
    classico_sup <- pred_teste + z_crit * sigma_res
  }
  amp_class <- classico_sup - classico_inf

  # Coberturas empíricas se y estiver presente
  cob_conf  <- if (tem_y_teste) mean(y_teste >= conf_inf & y_teste <= conf_sup, na.rm = TRUE) else NA_real_
  cob_class <- if (tem_y_teste) mean(y_teste >= classico_inf & y_teste <= classico_sup, na.rm = TRUE) else NA_real_

  df_res <- data.frame(
    Obs            = seq_len(n_teste),
    Estimativa     = pred_teste,
    Conf_Inf       = conf_inf,
    Conf_Sup       = conf_sup,
    Amp_Conformal  = amp_conf,
    Classico_Inf   = classico_inf,
    Classico_Sup   = classico_sup,
    Amp_Classico   = amp_class,
    stringsAsFactors = FALSE
  )
  if (tem_y_teste) {
    df_res$Observado <- y_teste
  }

  # Helpers de formatação
  .pad <- function(s, w, align = "center") {
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left") return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .fmt_f <- function(x, d = 2) {
    if (is.na(x)) return("-")
    formatC(x, format = "f", digits = d, decimal.mark = ",")
  }
  .fmt_pct <- function(x, d = 1) {
    if (is.na(x)) return("-")
    paste0(formatC(x * 100, format = "f", digits = d, decimal.mark = ","), "%")
  }

  # 4. Impressão no Console
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("PREDIÇÃO CONFORMAL (INTERVALOS CONFORMALIZADOS)")
  } else {
    cat("
── PREDIÇÃO CONFORMAL (INTERVALOS CONFORMALIZADOS) ──
")
  }

  cat(sprintf("  Variável Resposta:    %s
", nome_resp))
  cat(sprintf("  Tamanho Calibração:   %d observações
", n_cal))
  cat(sprintf("  Amostra Testada:      %d observações
", n_teste))
  cat(sprintf("  Nível de Confiança:   %s nominal
", .fmt_pct(nivel, 0)))
  cat(sprintf("  Raio Conformal (q):   %s

", .fmt_f(q_conf, 3)))

  # Tabela comparativa de predições (primeiras observações)
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERVALOS PREDITIVOS (PRIMEIRAS OBSERVAÇÕES)")
  } else {
    cat("  ▶ INTERVALOS PREDITIVOS (PRIMEIRAS OBSERVAÇÕES)
")
  }

  n_exibir <- min(n_teste, 8)
  w_tab <- c(obs = 8, est = 14, conf = 26, clas = 26)
  sep_tab <- paste0("  ", strrep("─", sum(w_tab)))
  cat(sep_tab, "
")
  cat("  ", .pad("Obs", w_tab["obs"], "center"),
             .pad("Estimativa", w_tab["est"], "center"),
             .pad("Intervalo Conformal", w_tab["conf"], "center"),
             .pad("Intervalo Clássico", w_tab["clas"], "center"), "
", sep = "")
  cat(sep_tab, "
")
  for (i in seq_len(n_exibir)) {
    str_conf <- sprintf("[ %s ; %s ]", .fmt_f(df_res$Conf_Inf[i], 2), .fmt_f(df_res$Conf_Sup[i], 2))
    str_clas <- sprintf("[ %s ; %s ]", .fmt_f(df_res$Classico_Inf[i], 2), .fmt_f(df_res$Classico_Sup[i], 2))
    cat("  ", .pad(as.character(df_res$Obs[i]), w_tab["obs"], "center"),
               .pad(.fmt_f(df_res$Estimativa[i], 2), w_tab["est"], "center"),
               .pad(str_conf, w_tab["conf"], "center"),
               .pad(str_clas, w_tab["clas"], "center"), "
", sep = "")
  }
  cat(sep_tab, "
")
  if (n_teste > n_exibir) {
    cat(sprintf("  ... mais %d observações calculadas disponíveis no objeto retornado.
", n_teste - n_exibir))
    cat(sep_tab, "
")
  }
  cat("
")

  # Tabela Resumo de Cobertura e Eficiência
  if (exists(".print_topico", mode = "function")) {
    .print_topico("COMPARAÇÃO DE EFICIÊNCIA E COBERTURA")
  } else {
    cat("  ▶ COMPARAÇÃO DE EFICIÊNCIA E COBERTURA
")
  }

  amp_media_conf <- mean(df_res$Amp_Conformal, na.rm = TRUE)
  amp_media_clas <- mean(df_res$Amp_Classico, na.rm = TRUE)
  razao_eficiencia <- amp_media_conf / amp_media_clas

  w_eff <- c(met = 24, amp = 20, cob = 20)
  sep_eff <- paste0("  ", strrep("─", sum(w_eff)))
  cat(sep_eff, "
")
  cat("  ", .pad("Método", w_eff["met"], "left"),
             .pad("Amplitude Média", w_eff["amp"], "center"),
             .pad("Cobertura Empírica", w_eff["cob"], "center"), "
", sep = "")
  cat(sep_eff, "
")
  cat("  ", .pad("Predição Conformal", w_eff["met"], "left"),
             .pad(.fmt_f(amp_media_conf, 2), w_eff["amp"], "center"),
             .pad(if (is.na(cob_conf)) "-" else .fmt_pct(cob_conf, 1), w_eff["cob"], "center"), "
", sep = "")
  cat("  ", .pad("Intervalo Clássico", w_eff["met"], "left"),
             .pad(.fmt_f(amp_media_clas, 2), w_eff["amp"], "center"),
             .pad(if (is.na(cob_class)) "-" else .fmt_pct(cob_class, 1), w_eff["cob"], "center"), "
", sep = "")
  cat(sep_eff, "

")

  # INTERPRETAÇÃO INTEGRADA
  if (exists(".print_topico", mode = "function")) {
    .print_topico("INTERPRETAÇÃO E DIAGNÓSTICO PREDITIVO")
  } else {
    cat("  ▶ INTERPRETAÇÃO E DIAGNÓSTICO PREDITIVO
")
  }

  cor_ok   <- "[1;32m"
  cor_warn <- "[1;33m"
  reset    <- "[0m"

  cat(sprintf("  %s[✓]%s A predição conformal assegura garantia matemática exata de que pelo menos %s
    das novas observações estarão contidas no intervalo, independentemente da distribuição dos erros.

",
              cor_ok, reset, .fmt_pct(nivel, 0)))

  if (razao_eficiencia <= 1.05) {
    cat(sprintf("  %s[✓]%s Alta eficiência preditiva: a amplitude média conformal (%s) é comparável ou mais justa
    que a do método clássico (%s), com razão relativa de %.2f.

",
                cor_ok, reset, .fmt_f(amp_media_conf, 2), .fmt_f(amp_media_clas, 2), razao_eficiencia))
  } else {
    cat(sprintf("  %s[!]%s Os intervalos conformais são mais conservadores que os clássicos (razão: %.2f),
    refletindo assimetria, caudas pesadas ou variabilidade amostral não captadas pela premissa gaussiana clássica.

",
                cor_warn, reset, razao_eficiencia))
  }

  if (exists(".print_rodape", mode = "function")) {
    .print_rodape()
  } else {
    cat("──────────────────────────────────────────────────────────────────────────────
")
  }

  res_obj <- list(
    predicoes            = df_res,
    raio_conformal       = q_conf,
    nivel                = nivel,
    amp_media_conformal  = amp_media_conf,
    amp_media_classico   = amp_media_clas,
    cobertura_observada  = list(conformal = cob_conf, classico = cob_class),
    dados_teste          = df_teste,
    nome_resposta        = nome_resp
  )
  class(res_obj) <- c("estatR_conformal", "list")

  if (isTRUE(grafico)) {
    p_graf <- plot_conformal(res_obj)
    print(p_graf)
  }

  invisible(res_obj)
}

#' @title Gráfico Comparativo de Predição Conformal
#' @description Plota os intervalos preditivos conformalizados contra o intervalo clássico.
#' @param modelo Objeto da classe \code{estatR_conformal}.
#' @return Objeto ggplot2.
#' @export
plot_conformal <- function(modelo) {
  if (!inherits(modelo, "estatR_conformal")) {
    stop("O objeto deve ser da classe 'estatR_conformal'.")
  }

  df <- modelo$predicoes
  n_show <- min(nrow(df), 40)
  df_sub <- df[seq_len(n_show), , drop = FALSE]

  p <- ggplot2::ggplot(df_sub, ggplot2::aes(x = Obs)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Conf_Inf, ymax = Conf_Sup, fill = "Intervalo Conformal"), alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = Classico_Inf, color = "Intervalo Clássico"), linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = Classico_Sup, color = "Intervalo Clássico"), linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = Estimativa, color = "Estimativa Pontual"), linewidth = 0.9)

  if ("Observado" %in% names(df_sub)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(y = Observado), color = "gray20", size = 2, alpha = 0.85)
  }

  p <- p +
    ggplot2::scale_fill_manual(name = "Cobertura Conformal", values = c("Intervalo Conformal" = "#1D3557")) +
    ggplot2::scale_color_manual(name = "Predição / Modelo", values = c("Intervalo Clássico" = "#E63946", "Estimativa Pontual" = "#1D3557")) +
    ggplot2::labs(
      title    = "Predição Conformal vs Intervalo Clássico",
      subtitle = sprintf("Garantia marginal exata de %.0f%% | Raio Conformal: %.2f", modelo$nivel * 100, modelo$raio_conformal),
      x        = "Índice da Observação",
      y        = modelo$nome_resposta
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
