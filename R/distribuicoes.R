# ── Modulo: distribuicao() ────────────────────────────────────────────────────
# Visao geral pedagogica de uma distribuicao: parametros, momentos,
# quantis, regra empirica (Normal), grafico pedagogico ou painel de simulacao.

#' @title Distribuição Estatística
#' @description Exibe uma visão completa e pedagógica de uma distribuição de probabilidade:
#' parâmetros, momentos teóricos, quantis-chave e gráfico ilustrativo.
#' Suporta: "normal", "binomial", "poisson", "exponencial", "uniforme", "t", "qui_quadrado", "f", "gama", "weibull".
#'
#' @param tipo Tipo da distribuição.
#' @param media,dp Parâmetros da Normal.
#' @param n,p Parâmetros da Binomial.
#' @param lambda Parâmetro da Poisson.
#' @param taxa Parâmetro da Exponencial.
#' @param min,max Parâmetros da Uniforme.
#' @param gl Graus de liberdade (t de Student e Qui-Quadrado).
#' @param gl1,gl2 Graus de liberdade (F de Snedecor).
#' @param forma,escala Parâmetros da Gama.
#' @param grafico Lógico. Se TRUE, exibe o gráfico pedagógico da distribuição.
#' @param simulacao Lógico. Se TRUE, exibe painel com histogramas simulados.
#' @param semente Semente para reprodutibilidade da simulação.
#' @export
distribuicao <- function(tipo = c("normal", "binomial", "poisson", "exponencial",
                                  "uniforme", "t", "qui_quadrado", "f", "gama", "weibull"),
                         media = 0, dp = 1, n = 10, p = 0.5, lambda = 1, taxa = 1,
                         min = 0, max = 1, gl = 10, gl1 = 5, gl2 = 10,
                         forma = 2, escala = 1,
                         grafico = TRUE, simulacao = FALSE, semente = 42, ...) {

  tipo <- match.arg(tipo)
  dots <- list(...)
  if ('df' %in% names(dots)) gl <- dots$df
  if ('df1' %in% names(dots)) gl1 <- dots$df1
  if ('df2' %in% names(dots)) gl2 <- dots$df2
  if ('shape' %in% names(dots)) forma <- dots$shape
  if ('scale' %in% names(dots)) escala <- dots$scale
  if ('rate' %in% names(dots)) taxa <- dots$rate
  if ('size' %in% names(dots)) n <- dots$size
  if ('prob' %in% names(dots)) p <- dots$prob

  # ── Helpers internos ──────────────────────────────────────────────────────
  .pad <- function(s, w, align = "left") {
    if (length(s) > 1) return(sapply(s, .pad, w = w, align = align))
    s <- trimws(as.character(s))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    if (align == "left")  return(paste0(s, strrep(" ", pad)))
    if (align == "right") return(paste0(strrep(" ", pad), s))
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  .f2 <- function(x) formatC(x, format = "f", digits = 2, decimal.mark = ",")
  .f3 <- function(x) { if (is.na(x) || is.nan(x)) return("—"); if (is.infinite(x)) return("Inf"); if (abs(x) >= 10) formatC(x, format="f", digits=2, decimal.mark=",") else if (abs(x) >= 1) formatC(x, format="f", digits=3, decimal.mark=",") else formatC(x, format="f", digits=4, decimal.mark=",") }

  # ── Calculos por tipo ─────────────────────────────────────────────────────
  tem_regra <- FALSE
  pcts <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
  
  if (tipo == "normal") {
    if (dp <= 0) stop("O desvio-padrão (dp) deve ser > 0.")
    teo_media  <- media; teo_var <- dp^2; teo_dp <- dp
    teo_assim  <- 0; teo_curt <- 0
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O NORMAL  \u2014  \u03bc = %g  |  \u03c3 = %g", media, dp)
    params_txt <- sprintf("  M\u00e9dia (\u03bc): %g    Desvio Padr\u00e3o (\u03c3): %g\n", media, dp)
    quantis <- qnorm(pcts, mean = media, sd = dp)
    regra68  <- c(media - dp,  media + dp)
    regra95  <- c(media - 2*dp, media + 2*dp)
    regra997 <- c(media - 3*dp, media + 3*dp)
    tem_regra <- TRUE
    pdf_formula <- sprintf("f(x) = (1 / (\u03c3\u221a2\u03c0)) \u00d7 exp(-\u00bd((x-\u03bc)/\u03c3)\u00b2)")

  } else if (tipo == "binomial") {
    if (n <= 0 || n != floor(n)) stop("'n' deve ser inteiro positivo.")
    if (p < 0 || p > 1)          stop("'p' deve estar entre 0 e 1.")
    teo_media  <- n * p; teo_var <- n * p * (1 - p); teo_dp <- sqrt(teo_var)
    teo_assim  <- (1 - 2*p) / teo_dp
    teo_curt   <- (1 - 6*p*(1-p)) / teo_var
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O BINOMIAL  \u2014  n = %d  |  p = %g", n, p)
    params_txt <- sprintf("  Ensaios (n): %d    Probabilidade (p): %g\n", n, p)
    quantis <- qbinom(pcts, n, p)
    pdf_formula <- sprintf("P(X=k) = C(n,k) \u00d7 p\u1d4f \u00d7 (1-p)^(n-k)")

  } else if (tipo == "poisson") {
    if (lambda <= 0) stop("'lambda' deve ser > 0.")
    teo_media  <- lambda; teo_var <- lambda; teo_dp <- sqrt(lambda)
    teo_assim  <- 1 / sqrt(lambda); teo_curt <- 1 / lambda
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O POISSON  \u2014  \u03bb = %g", lambda)
    params_txt <- sprintf("  Taxa (\u03bb): %g\n", lambda)
    quantis <- qpois(pcts, lambda)
    pdf_formula <- sprintf("P(X=k) = (e^-\u03bb \u00d7 \u03bb\u1d4f) / k!")

  } else if (tipo == "exponencial") {
    if (taxa <= 0) stop("'taxa' deve ser > 0.")
    teo_media  <- 1 / taxa; teo_var <- 1 / taxa^2; teo_dp <- 1 / taxa
    teo_assim  <- 2; teo_curt <- 6
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O EXPONENCIAL  \u2014  \u03bb = %g", taxa)
    params_txt <- sprintf("  Taxa (\u03bb): %g    M\u00e9dia (1/\u03bb): %g\n", taxa, 1/taxa)
    quantis <- qexp(pcts, rate = taxa)
    pdf_formula <- sprintf("f(x) = \u03bb \u00d7 e^(-\u03bb\u00d7x),  x \u2265 0")

  } else if (tipo == "uniforme") {
    if (min >= max) stop("'min' deve ser menor que 'max'.")
    teo_media <- (min + max) / 2; teo_var <- (max - min)^2 / 12; teo_dp <- sqrt(teo_var)
    teo_assim <- 0; teo_curt <- -1.2
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O UNIFORME  \u2014  U(%g, %g)", min, max)
    params_txt <- sprintf("  M\u00ednimo (a): %g    M\u00e1ximo (b): %g\n", min, max)
    quantis <- qunif(pcts, min, max)
    pdf_formula <- sprintf("f(x) = 1 / (b - a),  a \u2264 x \u2264 b")

  } else if (tipo == "t") {
    if (gl <= 0) stop("Graus de liberdade 'gl' deve ser > 0.")
    teo_media <- if (gl > 1) 0 else NaN
    teo_var   <- if (gl > 2) gl / (gl - 2) else if (gl > 1) Inf else NaN
    teo_dp    <- sqrt(teo_var)
    teo_assim <- if (gl > 3) 0 else NaN
    teo_curt  <- if (gl > 4) 6 / (gl - 4) else NaN
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O t DE STUDENT  \u2014  gl = %g", gl)
    params_txt <- sprintf("  Graus de Liberdade (\u03bd): %g\n", gl)
    quantis <- qt(pcts, df = gl)
    pdf_formula <- "f(x) \u221d (1 + x\u00b2/\u03bd)^(-(\u03bd+1)/2)"

  } else if (tipo == "qui_quadrado") {
    if (gl <= 0) stop("Graus de liberdade 'gl' deve ser > 0.")
    teo_media <- gl; teo_var <- 2 * gl; teo_dp <- sqrt(teo_var)
    teo_assim <- sqrt(8 / gl); teo_curt <- 12 / gl
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O QUI-QUADRADO (\u03c7\u00b2)  \u2014  gl = %g", gl)
    params_txt <- sprintf("  Graus de Liberdade (k): %g\n", gl)
    quantis <- qchisq(pcts, df = gl)
    pdf_formula <- "f(x) \u221d x^(k/2 - 1) \u00d7 e^(-x/2),  x \u2265 0"

  } else if (tipo == "f") {
    if (gl1 <= 0 || gl2 <= 0) stop("Graus de liberdade 'gl1' e 'gl2' devem ser > 0.")
    teo_media <- if (gl2 > 2) gl2 / (gl2 - 2) else NaN
    teo_var   <- if (gl2 > 4) 2 * gl2^2 * (gl1 + gl2 - 2) / (gl1 * (gl2 - 2)^2 * (gl2 - 4)) else NaN
    teo_dp    <- sqrt(teo_var)
    teo_assim <- if (gl2 > 6) (2*gl1 + gl2 - 2)*sqrt(8*(gl2-4))/( (gl2-6)*sqrt(gl1*(gl1+gl2-2)) ) else NaN
    teo_curt  <- NaN # Omitido para F complexa
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O F DE SNEDECOR  \u2014  gl1 = %g, gl2 = %g", gl1, gl2)
    params_txt <- sprintf("  Graus Liberdade Num (d1): %g    Graus Liberdade Den (d2): %g\n", gl1, gl2)
    quantis <- qf(pcts, df1 = gl1, df2 = gl2)
    pdf_formula <- "f(x) \u221d ( (d1 \u00d7 x / d2)^(d1/2) ) / (1 + d1 \u00d7 x / d2)^((d1+d2)/2),  x \u2265 0"

  } else if (tipo == "gama") {
    if (forma <= 0 || escala <= 0) stop("Parâmetros 'forma' e 'escala' devem ser > 0.")
    teo_media <- forma * escala; teo_var <- forma * escala^2; teo_dp <- sqrt(teo_var)
    teo_assim <- 2 / sqrt(forma); teo_curt <- 6 / forma
    titulo_dist <- sprintf("DISTRIBUI\u00c7\u00c3O GAMA  \u2014  forma (k) = %g, escala (\u03b8) = %g", forma, escala)
    params_txt <- sprintf("  Forma (k ou \u03b1): %g    Escala (\u03b8 ou 1/\u03b2): %g\n", forma, escala)
    quantis <- qgamma(pcts, shape = forma, scale = escala)
    pdf_formula <- "f(x) = (x^(k-1) \u00d7 e^(-x/\u03b8)) / (\u03b8\u1d4f \u00d7 \u0393(k)),  x \u2265 0"

  } else if (tipo == "weibull") {
    if (forma <= 0 || escala <= 0) stop("Parâmetros 'forma' e 'escala' devem ser > 0.")
    teo_media <- escala * gamma(1 + 1 / forma)
    teo_var   <- escala^2 * (gamma(1 + 2 / forma) - (gamma(1 + 1 / forma))^2)
    teo_dp    <- sqrt(teo_var)
    teo_assim <- (gamma(1 + 3/forma)*escala^3 - 3*teo_media*teo_var - teo_media^3) / (teo_dp^3)
    teo_curt  <- (gamma(1 + 4/forma)*escala^4 - 4*teo_media*gamma(1 + 3/forma)*escala^3 + 6*teo_media^2*teo_var + 3*teo_media^4) / (teo_var^2) - 3
    titulo_dist <- sprintf("DISTRIBUIÇÃO WEIBULL  —  forma (k) = %g, escala (λ) = %g", forma, escala)
    params_txt <- sprintf("  Forma (k): %g    Escala (λ): %g\n", forma, escala)
    quantis <- stats::qweibull(pcts, shape = forma, scale = escala)
    pdf_formula <- "f(x) = (k/λ) * (x/λ)^(k-1) * e^(-(x/λ)^k),  x ≥ 0"
  }

  # ── Saida Console ─────────────────────────────────────────────────────────
  if (exists(".print_titulo", mode = "function")) .print_titulo(titulo_dist) else cat(sprintf("\n\u2500\u2500 %s \u2500\u2500\n", titulo_dist))

  if (exists(".print_topico", mode = "function")) .print_topico("PAR\u00c2METROS")
  cat(params_txt)
  cat(sprintf("  F\u00f3rmula:  %s\n\n", pdf_formula))

  if (exists(".print_topico", mode = "function")) .print_topico("MOMENTOS TE\u00d3RICOS")
  wm <- 22; wv <- 12
  sep_m <- paste0("  ", strrep("\u2500", wm + wv + 1))
  cat(sep_m, "\n")
  cat("  ", .pad("Medida", wm, "left"), .pad("Valor", wv, "center"), "\n", sep = "")
  cat(sep_m, "\n")
  cat("  ", .pad("M\u00e9dia", wm, "left"),      .pad(.f3(teo_media), wv, "center"), "\n", sep = "")
  cat("  ", .pad("Vari\u00e2ncia", wm, "left"),  .pad(.f3(teo_var),   wv, "center"), "\n", sep = "")
  cat("  ", .pad("Desvio Padr\u00e3o", wm, "left"), .pad(.f3(teo_dp), wv, "center"), "\n", sep = "")
  cat("  ", .pad("Assimetria", wm, "left"),    .pad(.f3(teo_assim), wv, "center"), "\n", sep = "")
  cat("  ", .pad("Curtose (excesso)", wm, "left"), .pad(.f3(teo_curt), wv, "center"), "\n", sep = "")
  cat(sep_m, "\n\n")

  if (tem_regra) {
    if (exists(".print_topico", mode = "function")) .print_topico("REGRA EMP\u00cdRICA (68 \u2014 95 \u2014 99,7%)")
    cat(sprintf("  \u03bc \u00b1 1\u03c3  \u2192  [%s ; %s]  \u2192  68,27%% dos dados\n", .f2(regra68[1]),  .f2(regra68[2])))
    cat(sprintf("  \u03bc \u00b1 2\u03c3  \u2192  [%s ; %s]  \u2192  95,45%% dos dados\n", .f2(regra95[1]),  .f2(regra95[2])))
    cat(sprintf("  \u03bc \u00b1 3\u03c3  \u2192  [%s ; %s]  \u2192  99,73%% dos dados\n\n", .f2(regra997[1]), .f2(regra997[2])))
  }

  if (exists(".print_topico", mode = "function")) .print_topico("QUANTIS")
  p_labels <- c("P01","P05","P10","P25","P50","P75","P90","P95","P99")
  wq <- 8
  cat("  ", paste(.pad(p_labels, wq, "center"), collapse = ""), "\n", sep = "")
  cat("  ", strrep("\u2500", wq * length(p_labels)), "\n", sep = "")
  cat("  ", paste(.pad(.f2(quantis), wq, "center"), collapse = ""), "\n\n", sep = "")

  if (exists(".print_rodape", mode = "function")) .print_rodape()

  # ── Funcoes Auxiliares de Grafico ─────────────────────────────────────────
  .gera <- function(nn) {
    set.seed(semente + nn)
    switch(tipo,
      normal      = rnorm(nn, mean = media, sd = dp),
      binomial    = rbinom(nn, size = n, prob = p),
      poisson     = rpois(nn, lambda = lambda),
      exponencial = rexp(nn, rate = taxa),
      uniforme    = runif(nn, min = min, max = max),
      t           = rt(nn, df = gl),
      qui_quadrado= rchisq(nn, df = gl),
      f           = rf(nn, df1 = gl1, df2 = gl2),
      gama        = rgamma(nn, shape = forma, scale = escala),
      weibull     = rweibull(nn, shape = forma, scale = escala)
    )
  }
  
  .curva_teo <- function(gg, nn) {
    switch(tipo,
      normal = gg + ggplot2::stat_function(fun = dnorm, args = list(mean = media, sd = dp), color = "#D90429", linewidth = 1),
      exponencial = gg + ggplot2::stat_function(fun = dexp, args = list(rate = taxa), color = "#D90429", linewidth = 1),
      uniforme = gg + ggplot2::stat_function(fun = dunif, args = list(min = min, max = max), color = "#D90429", linewidth = 1),
      t = gg + ggplot2::stat_function(fun = dt, args = list(df = gl), color = "#D90429", linewidth = 1),
      qui_quadrado = gg + ggplot2::stat_function(fun = dchisq, args = list(df = gl), color = "#D90429", linewidth = 1),
      f = gg + ggplot2::stat_function(fun = df, args = list(df1 = gl1, df2 = gl2), color = "#D90429", linewidth = 1),
      gama = gg + ggplot2::stat_function(fun = dgamma, args = list(shape = forma, scale = escala), color = "#D90429", linewidth = 1),
      weibull = gg + ggplot2::stat_function(fun = dweibull, args = list(shape = forma, scale = escala), color = "#D90429", linewidth = 1),
      binomial = {
        am <- .gera(nn)
        xs <- min(am):max(am)
        df_teo <- data.frame(x = xs, y = dbinom(xs, size = n, prob = p))
        gg + ggplot2::geom_segment(data = df_teo, ggplot2::aes(x = x, xend = x, y = 0, yend = y), color = "#D90429", linewidth = 1) +
             ggplot2::geom_point(data = df_teo, ggplot2::aes(x = x, y = y), color = "#D90429", size = 2)
      },
      poisson = {
        am <- .gera(nn)
        xs <- min(am):max(am)
        df_teo <- data.frame(x = xs, y = dpois(xs, lambda = lambda))
        gg + ggplot2::geom_segment(data = df_teo, ggplot2::aes(x = x, xend = x, y = 0, yend = y), color = "#D90429", linewidth = 1) +
             ggplot2::geom_point(data = df_teo, ggplot2::aes(x = x, y = y), color = "#D90429", size = 2)
      }
    )
  }

  .tema <- function(gg, estilo = 1) {
    t_base <- if (exists("meu_tema", mode = "function")) {
      meu_tema(estilo = 1, grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      tema_estatR(estilo = 1, grade = "dupla")
    } else if (requireNamespace('estatR', quietly = TRUE)) {
      fn <- get('tema_estatR', envir = asNamespace('estatR'))
      fn(estilo = estilo, grade = 'dupla')
    } else {
      ggplot2::theme_minimal()
    }
    gg + t_base + ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold', size = 13))
  }

  # ── Modo Simulacao ────────────────────────────────────────────────────────
  if (simulacao && grafico) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      message("[Aviso] Instale o pacote 'patchwork' para o painel de simulacao.")
      return(invisible(NULL))
    }
    ns <- c(10, 50, 100, 250, 500, 1000)
    plots <- lapply(ns, function(nn) {
      df <- data.frame(X = .gera(nn))
      eh_continua <- !(tipo %in% c("binomial", "poisson"))
      if (eh_continua) {
        gg <- ggplot2::ggplot(df, ggplot2::aes(x = X)) +
          ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                  bins = max(8, min(30, round(sqrt(nn)))),
                                  fill = "#AAAAAA", color = "#555555", alpha = 0.7)
      } else {
        gg <- ggplot2::ggplot(df, ggplot2::aes(x = X)) +
          ggplot2::geom_bar(ggplot2::aes(y = ggplot2::after_stat(prop)), fill = "#AAAAAA", color = "#555555", alpha = 0.7)
      }
      gg <- .curva_teo(gg, nn) + ggplot2::labs(title = sprintf("n = %d", nn), x = "x", y = if (eh_continua) "Densidade" else "Proporção")
      .tema(gg)
    })

    painel <- (plots[[1]] | plots[[2]] | plots[[3]]) / (plots[[4]] | plots[[5]] | plots[[6]]) +
      patchwork::plot_annotation(title = sprintf("Simulação \u2014 %s", tools::toTitleCase(gsub("_", "-", tipo))),
                                 subtitle = "Convergência da distribuição amostral para a teórica (vermelha)")
    print(painel)
    return(invisible(NULL))
  }

  # ── Modo Grafico Estatico ─────────────────────────────────────────────────
  if (grafico && !simulacao) {
    if (tipo %in% c("normal", "t")) {
      eh_n <- tipo == "normal"
      media_plt <- if (eh_n) media else 0
      dp_plt <- if (eh_n) dp else teo_dp
      if (is.nan(dp_plt) || is.infinite(dp_plt)) dp_plt <- 2 # fallback grafico
      
      x_seq <- seq(media_plt - 4.2*dp_plt, media_plt + 4.2*dp_plt, length.out = 800)
      df_c <- data.frame(x = x_seq, y = if(eh_n) dnorm(x_seq, media, dp) else dt(x_seq, gl))
      
      .poly_region <- function(a, b) {
        sub <- df_c[df_c$x >= a & df_c$x <= b, ]
        rbind(data.frame(x = sub$x[1], y = 0), sub, data.frame(x = sub$x[nrow(sub)], y = 0))
      }
      
      p <- ggplot2::ggplot(df_c, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_polygon(data = .poly_region(media_plt - 3*dp_plt, media_plt + 3*dp_plt), ggplot2::aes(x = x, y = y), fill = "#F5DEB3", alpha = 1) +
        ggplot2::geom_polygon(data = .poly_region(media_plt - 2*dp_plt, media_plt + 2*dp_plt), ggplot2::aes(x = x, y = y), fill = "#D4A843", alpha = 1) +
        ggplot2::geom_polygon(data = .poly_region(media_plt - dp_plt, media_plt + dp_plt), ggplot2::aes(x = x, y = y), fill = "#B8860B", alpha = 1) +
        ggplot2::geom_line(color = "#8B6914", linewidth = 1.4)
        
      if (eh_n) {
        breaks_x <- c(media - 3*dp, media - 2*dp, media - dp, media, media + dp, media + 2*dp, media + 3*dp)
        labs_x <- c(expression(mu - 3*sigma), expression(mu - 2*sigma), expression(mu - sigma), expression(mu), expression(mu + sigma), expression(mu + 2*sigma), expression(mu + 3*sigma))
        p <- p + ggplot2::scale_x_continuous(breaks = breaks_x, labels = labs_x) +
          ggplot2::annotate("text", x = media - 0.5*dp, y = max(df_c$y) * 0.35, label = "34,13%", size = 5.5, color = "white", fontface = "bold") +
          ggplot2::annotate("text", x = media + 0.5*dp, y = max(df_c$y) * 0.35, label = "34,13%", size = 5.5, color = "white", fontface = "bold") +
          ggplot2::annotate("text", x = media - 1.5*dp, y = max(df_c$y) * 0.12, label = "13,59%", size = 4.5, color = "#4A3000", fontface = "bold") +
          ggplot2::annotate("text", x = media + 1.5*dp, y = max(df_c$y) * 0.12, label = "13,59%", size = 4.5, color = "#4A3000", fontface = "bold") +
          ggplot2::annotate("text", x = media - 2.5*dp, y = max(df_c$y) * 0.02, label = "2,14%", size = 3.8, color = "#4A3000", vjust = 1.5) +
          ggplot2::annotate("text", x = media + 2.5*dp, y = max(df_c$y) * 0.02, label = "2,14%", size = 3.8, color = "#4A3000", vjust = 1.5) +
          ggplot2::geom_vline(xintercept = media, color = '#8B6914', linetype = 'dashed', linewidth = 0.8) +
          ggplot2::coord_cartesian(clip = 'off') +
          ggplot2::labs(title = sprintf("Distribuição Normal  —  N(μ = %g, σ = %g)", media, dp), subtitle = "Regra Empírica: 68% — 95% — 99,7%", x = "x", y = "f(x)")
      } else {
        p <- p + ggplot2::labs(title = sprintf("Distribui\u00e7\u00e3o t de Student  \u2014  t(\u03bd = %g)", gl), subtitle = sprintf("Formato de sino, caudas mais pesadas que a Normal  |  gl = %g", gl), x = "x", y = "f(x)")
      }
      p <- .tema(p)
      
    } else if (tipo %in% c("exponencial", "qui_quadrado", "f", "gama", "weibull")) {
      if (tipo == "exponencial") {
        x_max <- qexp(0.999, taxa); f_den <- function(x) dexp(x, taxa); med_x <- 1/taxa
        tit <- sprintf("Distribui\u00e7\u00e3o Exponencial  \u2014  Exp(\u03bb = %g)", taxa)
      } else if (tipo == "qui_quadrado") {
        x_max <- qchisq(0.999, gl); f_den <- function(x) dchisq(x, gl); med_x <- gl
        tit <- sprintf("Distribui\u00e7\u00e3o Qui-Quadrado  \u2014  \u03c7\u00b2(k = %g)", gl)
      } else if (tipo == "f") {
        x_max <- qf(0.95, gl1, gl2); f_den <- function(x) df(x, gl1, gl2); med_x <- teo_media
        tit <- sprintf("Distribuição F de Snedecor  —  F(gl1 = %g, gl2 = %g)", gl1, gl2)
      } else if (tipo == "weibull") {
        x_max <- stats::qweibull(0.999, shape = forma, scale = escala)
        f_den <- function(x) stats::dweibull(x, shape = forma, scale = escala)
        med_x <- teo_media
        tit <- sprintf("Distribui\u00e7\u00e3o Weibull  \u2014  Weibull(k = %g, \u03bb = %g)", forma, escala)
      } else {
        x_max <- qgamma(0.999, shape = forma, scale = escala); f_den <- function(x) dgamma(x, shape=forma, scale=escala); med_x <- teo_media
        tit <- sprintf("Distribui\u00e7\u00e3o Gama  \u2014  Gama(\u03b1 = %g, \u03b8 = %g)", forma, escala)
      }
      
      x_seq <- seq(0.001, x_max * 1.05, length.out = 500)
      df_c <- data.frame(x = x_seq, y = f_den(x_seq))
      
      p <- ggplot2::ggplot(df_c, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_area(fill = "#AACCFF", alpha = 0.5) +
        ggplot2::geom_line(color = "#1A5DB5", linewidth = 1.2)
        
      if (!is.nan(med_x)) {
        p <- p + ggplot2::geom_vline(xintercept = med_x, color = "#D90429", linetype = "dashed", linewidth = 0.9) +
          ggplot2::annotate("text", x = med_x + x_max * 0.03, y = max(df_c$y) * 0.7, label = sprintf("\u03bc = %g", med_x), color = "#D90429", size = 6, hjust = 0, fontface = "bold")
      }
      p <- p + ggplot2::labs(title = tit, subtitle = if(!is.nan(med_x)) sprintf("Média = %g", med_x) else "Média indefinida", x = "x", y = "f(x)")
      p <- .tema(p)
      
    } else if (tipo == "uniforme") {
      df_c <- data.frame(x = c(min - 0.5, min, min, max, max, max + 0.5), y = c(0, 0, 1/(max-min), 1/(max-min), 0, 0))
      p <- ggplot2::ggplot(df_c, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_area(fill = "#AACCFF", alpha = 0.5) +
        ggplot2::geom_line(color = "#1A5DB5", linewidth = 1.2) +
        ggplot2::labs(title = sprintf("Distribuição Uniforme  \u2014  U(%g, %g)", min, max), subtitle = sprintf("Área = 1, Altura = %g", 1/(max-min)), x = "x", y = "f(x)")
      p <- .tema(p)
      
    } else if (tipo %in% c("binomial", "poisson")) {
      if (tipo == "binomial") {
        xs <- 0:n; ys <- dbinom(xs, n, p); med_d <- n * p
        tit <- sprintf("Distribuição Binomial  —  Bin(n = %d, p = %g)", n, p)
      } else {
        xs <- 0:qpois(0.999, lambda); ys <- dpois(xs, lambda); med_d <- lambda
        tit <- sprintf("Distribuição Poisson  —  Poisson(λ = %g)", lambda)
      }
      df_d <- data.frame(x = xs, y = ys)
      p <- ggplot2::ggplot(df_d, ggplot2::aes(x = factor(x), y = y)) +
        ggplot2::geom_col(fill = ifelse(df_d$x == round(med_d), "#1A5DB5", "#AACCFF"), color = "#555555", alpha = 0.85) +
        ggplot2::geom_point(ggplot2::aes(x = as.numeric(factor(x))), color = "#D90429", size = 2.5) +
        ggplot2::geom_line(ggplot2::aes(x = as.numeric(factor(x))), color = "#D90429", linewidth = 1, group = 1) +
        ggplot2::labs(title = tit, subtitle = sprintf("Média = %g", med_d), x = "x", y = "P(X = x)")
      p <- .tema(p)
    }
    suppressMessages(print(p))
  }
  invisible(list(media = teo_media, variancia = teo_var, dp = teo_dp, quantis = stats::setNames(quantis, p_labels)))
}