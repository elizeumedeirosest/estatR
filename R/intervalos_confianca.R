# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: INTERVALOS DE CONFIANÇA
# ─────────────────────────────────────────────────────────────────────────────

# ── Helpers internos ──────────────────────────────────────────────────────────
.ic_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
  paste0(s, strrep(" ", sp))
}

.ic_sep <- function(w) paste0("  ", strrep("\u2500", w))

.ic_fmt <- function(x, decimais = 3) {
  formatC(x, format = "f", digits = decimais, decimal.mark = ",")
}

.plot_ic <- function(est, linf, lsup, titulo, subtitulo, x_label) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  if (!exists("meu_tema")) return(NULL)
  
  df_plot <- data.frame(x = est, y = 1, linf = linf, lsup = lsup)
  
  amp <- lsup - linf
  
  old_w <- getOption("warn"); options(warn = -1)
  p <- ggplot2::ggplot(df_plot) +
    ggplot2::geom_segment(ggplot2::aes(x = linf, xend = lsup, y = 1, yend = 1),
                          color = "#555555", linewidth = 1.2) +
    ggplot2::geom_segment(ggplot2::aes(x = linf, xend = linf, y = 0.95, yend = 1.05),
                          color = "#555555", linewidth = 1) +
    ggplot2::geom_segment(ggplot2::aes(x = lsup, xend = lsup, y = 0.95, yend = 1.05),
                          color = "#555555", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = y),
                        color = "#D90429", size = 6, shape = 18) +
    ggplot2::annotate("text", x = linf, y = 1.10, label = .ic_fmt(linf, 2), 
                      size = 4.5, fontface = "bold", color = "#333333") +
    ggplot2::annotate("text", x = lsup, y = 1.10, label = .ic_fmt(lsup, 2), 
                      size = 4.5, fontface = "bold", color = "#333333") +
    ggplot2::annotate("text", x = est,  y = 0.90, label = .ic_fmt(est, 2), 
                      size = 4.5, fontface = "bold", color = "#D90429") +
    ggplot2::scale_y_continuous(limits = c(0.7, 1.3)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
    ggplot2::labs(title = titulo, subtitle = subtitulo, x = x_label, y = NULL, caption = "estatR") +
    meu_tema(grade = "dupla") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank())
  
  suppressMessages(print(p))
  options(warn = old_w)
}


# ─────────────────────────────────────────────────────────────────────────────
#' @title Intervalo de Confiança para a Média
#' @description Calcula o intervalo de confiança para a média de uma amostra.
#' @param x Vetor numérico.
#' @param confianca Nível de confiança desejado (padrão 0.95).
#' @param decimais Casas decimais para exibição (padrão 3).
#' @param grafico Lógico. Se TRUE, plota o intervalo com a estimativa pontual.
#' @return Retorna invisivelmente um data frame com os resultados.
#' @export
ic_media <- function(x, confianca = 0.95, decimais = 3, grafico = TRUE) {
  if (!is.numeric(x)) stop("'x' deve ser num\u00e9rico.")
  x_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x")) var_nome <- "x"
  
  x_c <- x[!is.na(x)]
  n   <- length(x_c)
  if (n < 2) stop("A amostra deve ter pelo menos 2 observa\u00e7\u00f5es v\u00e1lidas.")
  
  media <- mean(x_c)
  dp    <- sd(x_c)
  ep    <- dp / sqrt(n)
  gl    <- n - 1
  
  alpha  <- 1 - confianca
  tc     <- qt(1 - alpha / 2, df = gl)
  margem <- tc * ep
  
  linf <- media - margem
  lsup <- media + margem
  
  w_sep <- 67
  cat(sprintf("\n\u2500\u2500 INTERVALO DE CONFIAN\u00c7A PARA A M\u00c9DIA %s\n", strrep("\u2500", w_sep - 39)))
  cat(sprintf("  Vari\u00e1vel: %s   |   N v\u00e1lido: %d   |   Confian\u00e7a: %.0f%%\n", var_nome, n, confianca * 100))
  cat("  Distribui\u00e7\u00e3o utilizada: t de Student (vari\u00e2ncia desconhecida)\n\n")
  
  cat("  \u25b6 ESTIMATIVAS\n")
  wc <- c(ep=20, err=15, margem=16, gl=12)
  hdr <- paste0(.ic_pad("Estimativa Pontual", wc["ep"], "center"),
                .ic_pad("Erro Padr\u00e3o", wc["err"], "center"),
                .ic_pad("Margem de Erro", wc["margem"], "center"),
                .ic_pad("Graus Lib.", wc["gl"], "center"))
  cat(.ic_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.ic_sep(nchar(hdr)), "\n")
  cat("  ",
      .ic_pad(.ic_fmt(media, decimais), wc["ep"], "center"),
      .ic_pad(.ic_fmt(ep, decimais), wc["err"], "center"),
      .ic_pad(.ic_fmt(margem, decimais), wc["margem"], "center"),
      .ic_pad(gl, wc["gl"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr)), "\n\n")
  
  cat(sprintf("  \u25b6 INTERVALO DE CONFIAN\u00c7A (%.0f%%)\n", confianca * 100))
  cat(.ic_sep(nchar(hdr)), "\n")
  cat(sprintf("  %s\n", .ic_pad(sprintf("[ %s  ;  %s ]", .ic_fmt(linf, decimais), .ic_fmt(lsup, decimais)), nchar(hdr), "center")))
  cat(.ic_sep(nchar(hdr)), "\n\n")
  
  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  cat(sprintf("  Com %.0f%% de confian\u00e7a, estima-se que a verdadeira m\u00e9dia\n", confianca * 100))
  cat(sprintf("  populacional de '%s' est\u00e1 entre %s e %s.\n", var_nome, .ic_fmt(linf, decimais), .ic_fmt(lsup, decimais)))
  cat(strrep("\u2500", w_sep), "\n\n")
  
  if (grafico) {
    tryCatch({
      .plot_ic(media, linf, lsup, 
               sprintf("Intervalo de Confian\u00e7a (%.0f%%) para a M\u00e9dia", confianca * 100), 
               sprintf("Vari\u00e1vel: %s", var_nome), var_nome)
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }
  
  invisible(data.frame(media = media, erro_padrao = ep, margem = margem, li = linf, ls = lsup))
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Intervalo de Confiança para a Proporção
#' @description Calcula o intervalo de confiança para uma proporção (aproximação normal/Wald).
#' @param x Número de sucessos.
#' @param n Tamanho da amostra.
#' @param confianca Nível de confiança desejado (padrão 0.95).
#' @param decimais Casas decimais para exibição (padrão 3).
#' @param grafico Lógico. Se TRUE, plota o intervalo com a estimativa pontual.
#' @return Retorna invisivelmente um data frame com os resultados.
#' @export
ic_proporcao <- function(x, n, confianca = 0.95, decimais = 3, grafico = TRUE) {
  if (x > n) stop("O n\u00famero de sucessos (x) n\u00e3o pode ser maior que o tamanho da amostra (n).")
  
  p_hat <- x / n
  ep    <- sqrt(p_hat * (1 - p_hat) / n)
  alpha <- 1 - confianca
  zc    <- qnorm(1 - alpha / 2)
  margem <- zc * ep
  
  linf <- max(0, p_hat - margem)
  lsup <- min(1, p_hat + margem)
  
  w_sep <- 67
  cat(sprintf("\n\u2500\u2500 INTERVALO DE CONFIAN\u00c7A PARA A PROPOR\u00c7\u00c3O %s\n", strrep("\u2500", w_sep - 42)))
  cat(sprintf("  Sucessos: %d   |   N total: %d   |   Confian\u00e7a: %.0f%%\n", x, n, confianca * 100))
  cat("  M\u00e9todo: Aproxima\u00e7\u00e3o Normal (Wald)\n\n")
  
  cat("  \u25b6 ESTIMATIVAS\n")
  wc <- c(ep=22, err=15, margem=16, z=10)
  hdr <- paste0(.ic_pad("Propor\u00e7\u00e3o Amostral", wc["ep"], "center"),
                .ic_pad("Erro Padr\u00e3o", wc["err"], "center"),
                .ic_pad("Margem de Erro", wc["margem"], "center"),
                .ic_pad("Valor Z", wc["z"], "center"))
  cat(.ic_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.ic_sep(nchar(hdr)), "\n")
  cat("  ",
      .ic_pad(sprintf("%s (%.1f%%)", .ic_fmt(p_hat, decimais), p_hat*100), wc["ep"], "center"),
      .ic_pad(.ic_fmt(ep, decimais), wc["err"], "center"),
      .ic_pad(sprintf("%s (%.1f%%)", .ic_fmt(margem, decimais), margem*100), wc["margem"], "center"),
      .ic_pad(.ic_fmt(zc, 2), wc["z"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr)), "\n\n")
  
  cat(sprintf("  \u25b6 INTERVALO DE CONFIAN\u00c7A (%.0f%%)\n", confianca * 100))
  cat(.ic_sep(nchar(hdr)), "\n")
  str_ic <- sprintf("[ %s%%  ;  %s%% ]", .ic_fmt(linf*100, 1), .ic_fmt(lsup*100, 1))
  cat(sprintf("  %s\n", .ic_pad(str_ic, nchar(hdr), "center")))
  cat(.ic_sep(nchar(hdr)), "\n\n")
  
  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  cat(sprintf("  Com %.0f%% de confian\u00e7a, estima-se que a verdadeira propor\u00e7\u00e3o\n", confianca * 100))
  cat(sprintf("  populacional est\u00e1 entre %s%% e %s%%.\n", .ic_fmt(linf*100, 1), .ic_fmt(lsup*100, 1)))
  cat(strrep("\u2500", w_sep), "\n\n")
  
  if (grafico) {
    tryCatch({
      .plot_ic(p_hat, linf, lsup, 
               sprintf("Intervalo de Confian\u00e7a (%.0f%%) para a Propor\u00e7\u00e3o", confianca * 100), 
               sprintf("Estimativa: %s%% \u00b1 %s%%", .ic_fmt(p_hat*100, 1), .ic_fmt(margem*100, 1)), 
               "Propor\u00e7\u00e3o")
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }
  
  invisible(data.frame(proporcao = p_hat, erro_padrao = ep, margem = margem, li = linf, ls = lsup))
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Intervalo de Confiança para a Variância
#' @description Calcula o intervalo de confiança para a variância e desvio padrão.
#' @param x Vetor numérico.
#' @param confianca Nível de confiança desejado (padrão 0.95).
#' @param decimais Casas decimais para exibição (padrão 3).
#' @return Retorna invisivelmente um data frame com os resultados.
#' @export
ic_variancia <- function(x, confianca = 0.95, decimais = 3) {
  if (!is.numeric(x)) stop("'x' deve ser num\u00e9rico.")
  x_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x")) var_nome <- "x"
  
  x_c <- x[!is.na(x)]
  n   <- length(x_c)
  if (n < 2) stop("A amostra deve ter pelo menos 2 observa\u00e7\u00f5es v\u00e1lidas.")
  
  v  <- var(x_c)
  dp <- sd(x_c)
  gl <- n - 1
  
  alpha <- 1 - confianca
  chi_inf <- qchisq(1 - alpha / 2, df = gl)
  chi_sup <- qchisq(alpha / 2, df = gl)
  
  linf_v <- (gl * v) / chi_inf
  lsup_v <- (gl * v) / chi_sup
  
  linf_dp <- sqrt(linf_v)
  lsup_dp <- sqrt(lsup_v)
  
  w_sep <- 67
  cat(sprintf("\n\u2500\u2500 INTERVALO DE CONFIAN\u00c7A PARA A VARI\u00c2NCIA %s\n", strrep("\u2500", w_sep - 42)))
  cat(sprintf("  Vari\u00e1vel: %s   |   N v\u00e1lido: %d   |   Confian\u00e7a: %.0f%%\n", var_nome, n, confianca * 100))
  cat("  Distribui\u00e7\u00e3o utilizada: Qui-Quadrado (\u03c7\u00b2)\n\n")
  
  cat("  \u25b6 ESTIMATIVAS PONTUAIS\n")
  wc <- c(met=20, est=20)
  hdr <- paste0(.ic_pad("M\u00e9trica", wc["met"], "center"), .ic_pad("Valor", wc["est"], "center"))
  cat(.ic_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.ic_sep(nchar(hdr)), "\n")
  cat("  ", .ic_pad("Vari\u00e2ncia", wc["met"], "center"), .ic_pad(.ic_fmt(v, decimais), wc["est"], "center"), "\n", sep = "")
  cat("  ", .ic_pad("Desvio Padr\u00e3o", wc["met"], "center"), .ic_pad(.ic_fmt(dp, decimais), wc["est"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr)), "\n\n")
  
  cat(sprintf("  \u25b6 INTERVALOS DE CONFIAN\u00c7A (%.0f%%)\n", confianca * 100))
  w_ic <- c(met=20, int=30)
  hdr_ic <- paste0(.ic_pad("Par\u00e2metro", w_ic["met"], "center"), .ic_pad("Intervalo [ LI ; LS ]", w_ic["int"], "center"))
  cat(.ic_sep(nchar(hdr_ic)), "\n")
  cat("  ", hdr_ic, "\n", sep = "")
  cat(.ic_sep(nchar(hdr_ic)), "\n")
  cat("  ", .ic_pad("Vari\u00e2ncia (\u03c3\u00b2)", w_ic["met"], "center"), 
      .ic_pad(sprintf("[ %s ; %s ]", .ic_fmt(linf_v, decimais), .ic_fmt(lsup_v, decimais)), w_ic["int"], "center"), "\n", sep = "")
  cat("  ", .ic_pad("Desvio Padr\u00e3o (\u03c3)", w_ic["met"], "center"), 
      .ic_pad(sprintf("[ %s ; %s ]", .ic_fmt(linf_dp, decimais), .ic_fmt(lsup_dp, decimais)), w_ic["int"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr_ic)), "\n\n")
  
  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  cat(sprintf("  Com %.0f%% de confian\u00e7a, estima-se que a variabilidade da\n", confianca * 100))
  cat(sprintf("  popula\u00e7\u00e3o, medida pelo desvio padr\u00e3o, est\u00e1 entre %s e %s.\n", .ic_fmt(linf_dp, decimais), .ic_fmt(lsup_dp, decimais)))
  cat(strrep("\u2500", w_sep), "\n\n")
  
  invisible(data.frame(variancia = v, desvio_padrao = dp, li_var = linf_v, ls_var = lsup_v, li_dp = linf_dp, ls_dp = lsup_dp))
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Intervalo de Confiança para a Diferença de Médias
#' @description Calcula o intervalo de confiança para a diferença entre médias de dois grupos.
#' @param x Vetor numérico ou data frame.
#' @param grupo Vetor com os grupos (2 níveis) ou nome da coluna categórica se x for data frame.
#' @param variancia_igual Lógico. Assume variâncias iguais? (Padrão FALSE = Welch).
#' @param confianca Nível de confiança desejado (padrão 0.95).
#' @param decimais Casas decimais para exibição (padrão 3).
#' @param grafico Lógico. Se TRUE, plota o intervalo com a estimativa pontual.
#' @return Retorna invisivelmente um data frame com os resultados.
#' @export
ic_diferenca_medias <- function(x, grupo, variancia_igual = FALSE, confianca = 0.95, decimais = 3, grafico = TRUE) {
  grupo_sub  <- substitute(grupo)
  grupo_nome <- sub(".*\\$", "", deparse(grupo_sub))
  x_sub      <- substitute(x)
  x_nome     <- sub(".*\\$", "", deparse(x_sub))
  
  if (is.data.frame(x)) {
    grupo_vec <- if (grupo_nome %in% names(x)) x[[grupo_nome]] else eval(grupo_sub, parent.frame())
    x_vec     <- x[[x_nome]]
    var_nome  <- x_nome
  } else {
    grupo_vec <- eval(grupo_sub, parent.frame())
    x_vec     <- x
    var_nome  <- x_nome
  }
  
  grupo_vec <- as.factor(grupo_vec)
  niveis    <- levels(grupo_vec)
  if (length(niveis) != 2) stop("O agrupamento deve conter exatamente 2 n\u00edveis.")
  
  res <- t.test(x_vec ~ grupo_vec, var.equal = variancia_igual, conf.level = confianca)
  
  media1 <- res$estimate[1]
  media2 <- res$estimate[2]
  dif    <- media1 - media2
  linf   <- res$conf.int[1]
  lsup   <- res$conf.int[2]
  gl     <- res$parameter
  
  tipo_txt <- if (variancia_igual) "Student (vari\u00e2ncias iguais)" else "Welch (vari\u00e2ncias n\u00e3o iguais)"
  
  w_sep <- 67
  cat(sprintf("\n\u2500\u2500 IC PARA A DIFEREN\u00c7A DE M\u00c9DIAS %s\n", strrep("\u2500", w_sep - 32)))
  cat(sprintf("  Vari\u00e1vel: %s   |   Grupos: %s vs %s\n", var_nome, niveis[1], niveis[2]))
  cat(sprintf("  M\u00e9todo: Teste t de %s\n\n", tipo_txt))
  
  cat("  \u25b6 ESTIMATIVAS PONTUAIS\n")
  wc <- c(grp=15, med=15)
  hdr_g <- paste0(.ic_pad("Grupo", wc["grp"], "center"), .ic_pad("M\u00e9dia", wc["med"], "center"))
  cat(.ic_sep(nchar(hdr_g)), "\n")
  cat("  ", hdr_g, "\n", sep = "")
  cat(.ic_sep(nchar(hdr_g)), "\n")
  cat("  ", .ic_pad(niveis[1], wc["grp"], "center"), .ic_pad(.ic_fmt(media1, decimais), wc["med"], "center"), "\n", sep = "")
  cat("  ", .ic_pad(niveis[2], wc["grp"], "center"), .ic_pad(.ic_fmt(media2, decimais), wc["med"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr_g)), "\n\n")
  
  cat(sprintf("  \u25b6 INTERVALO DE CONFIAN\u00c7A PARA A DIFEREN\u00c7A (%.0f%%)\n", confianca * 100))
  w_ic <- c(dif=18, int=30)
  hdr_ic <- paste0(.ic_pad("Diferen\u00e7a (\u0394)", w_ic["dif"], "center"), .ic_pad("Intervalo [ LI ; LS ]", w_ic["int"], "center"))
  cat(.ic_sep(nchar(hdr_ic)), "\n")
  cat("  ", hdr_ic, "\n", sep = "")
  cat(.ic_sep(nchar(hdr_ic)), "\n")
  cat("  ", .ic_pad(.ic_fmt(dif, decimais), w_ic["dif"], "center"), 
      .ic_pad(sprintf("[ %s ; %s ]", .ic_fmt(linf, decimais), .ic_fmt(lsup, decimais)), w_ic["int"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr_ic)), "\n\n")
  
  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
  if (linf * lsup > 0) {
    cat(sprintf("  Como o intervalo n\u00e3o cont\u00e9m o zero, h\u00e1 diferen\u00e7a significativa\n"))
    cat(sprintf("  entre as m\u00e9dias dos grupos (ao n\u00edvel de %.0f%% de confian\u00e7a).\n", confianca*100))
  } else {
    cat(sprintf("  Como o intervalo cont\u00e9m o zero, a diferen\u00e7a entre as m\u00e9dias\n"))
    cat(sprintf("  n\u00e3o \u00e9 estatisticamente significativa (ao n\u00edvel de %.0f%% de confian\u00e7a).\n", confianca*100))
  }
  cat(strrep("\u2500", w_sep), "\n\n")
  
  if (grafico) {
    tryCatch({
      .plot_ic(dif, linf, lsup, 
               sprintf("Intervalo de Confian\u00e7a (%.0f%%) para a Diferen\u00e7a", confianca * 100), 
               sprintf("%s - %s", niveis[1], niveis[2]), 
               "\u0394 M\u00e9dia")
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }
  
  invisible(data.frame(diferenca = dif, li = linf, ls = lsup))
}
