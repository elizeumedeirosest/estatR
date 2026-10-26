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
    ggplot2::labs(title = titulo, subtitle = subtitulo, x = x_label, y = NULL) +
    meu_tema(estilo = 1, grade = "dupla") +
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
  if (!is.numeric(x)) stop("'x' deve ser numérico.")
  x_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x")) var_nome <- "x"
  
  x_c <- x[!is.na(x)]
  n   <- length(x_c)
  if (n < 2) stop("A amostra deve ter pelo menos 2 observações válidas.")
  
  media <- mean(x_c)
  dp    <- sd(x_c)
  ep    <- dp / sqrt(n)
  gl    <- n - 1
  
  alpha  <- 1 - confianca
  tc     <- qt(1 - alpha / 2, df = gl)
  margem <- tc * ep
  
  linf <- media - margem
  lsup <- media + margem
  
  .print_titulo("INTERVALO DE CONFIANÇA PARA A MÉDIA")
  cat(sprintf("  Variável: %s   |   N válido: %d   |   Confiança: %.0f%%\n", var_nome, n, confianca * 100))
  cat("  Distribuição utilizada: t de Student (variância desconhecida)\n\n")
  
  .print_topico("ESTIMATIVAS")
  wc <- c(ep=20, err=15, margem=16, gl=12)
  hdr <- paste0(.ic_pad("Estimativa Pontual", wc["ep"], "center"),
                .ic_pad("Erro Padrão", wc["err"], "center"),
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
  
  .print_topico(sprintf("INTERVALO DE CONFIANÇA (%.0f%%)", confianca * 100))
  cat(.ic_sep(nchar(hdr)), "\n")
  cat(sprintf("  %s\n", .ic_pad(sprintf("[ %s  ;  %s ]", .ic_fmt(linf, decimais), .ic_fmt(lsup, decimais)), nchar(hdr), "center")))
  cat(.ic_sep(nchar(hdr)), "\n\n")
  
  .print_topico("INTERPRETAÇÃO")
  cat(sprintf("  Com %.0f%% de confiança, estima-se que a verdadeira média\n", confianca * 100))
  cat(sprintf("  populacional de '%s' está entre %s e %s.\n\n", var_nome, .ic_fmt(linf, decimais), .ic_fmt(lsup, decimais)))
  cat(sprintf("  Segundo a interpretação frequentista, podemos dizer que se repetirmos\n"))
  cat(sprintf("  o experimento diversas vezes, %.0f%% dos intervalos de confiança\n", confianca * 100))
  cat(sprintf("  construídos a partir dessas amostras irão conter o valor verdadeiro\n"))
  cat(sprintf("  do parâmetro populacional.\n"))
  .print_rodape()
  
  if (grafico) {
    tryCatch({
      .plot_ic(media, linf, lsup, 
               sprintf("Intervalo de Confiança (%.0f%%) para a Média", confianca * 100), 
               sprintf("Variável: %s", var_nome), var_nome)
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
#' @param variavel Nome descritivo da variável para exibição (opcional).
#' @return Retorna invisivelmente um data frame com os resultados.
#' @export
ic_proporcao <- function(x, n, confianca = 0.95, decimais = 3, grafico = TRUE, variavel = NULL) {
  if (x > n) stop("O número de sucessos (x) não pode ser maior que o tamanho da amostra (n).")
  var_label <- if (!is.null(variavel)) as.character(variavel) else "Proporção Amostral"
  
  p_hat <- x / n
  ep    <- sqrt(p_hat * (1 - p_hat) / n)
  alpha <- 1 - confianca
  zc    <- qnorm(1 - alpha / 2)
  margem <- zc * ep
  
  linf <- max(0, p_hat - margem)
  lsup <- min(1, p_hat + margem)
  
  .print_titulo("INTERVALO DE CONFIANÇA PARA A PROPORÇÃO")
  cat(sprintf("  Sucessos: %d   |   N total: %d   |   Confiança: %.0f%%\n", x, n, confianca * 100))
  cat("  Método: Aproximação Normal (Wald)\n\n")
  
  .print_topico("ESTIMATIVAS")
  wc <- c(ep = 24, err = 16, margem = 20, z = 12)
  hdr <- paste0(.ic_pad("Proporção Amostral", wc["ep"], "center"),
                .ic_pad("Erro Padrão", wc["err"], "center"),
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
  
  .print_topico(sprintf("INTERVALO DE CONFIANÇA (%.0f%%)", confianca * 100))
  cat(.ic_sep(nchar(hdr)), "\n")
  str_ic <- sprintf("[ %s%%  ;  %s%% ]", .ic_fmt(linf*100, 1), .ic_fmt(lsup*100, 1))
  cat(sprintf("  %s\n", .ic_pad(str_ic, nchar(hdr), "center")))
  cat(.ic_sep(nchar(hdr)), "\n\n")
  
  .print_topico("INTERPRETAÇÃO")
  cat(sprintf("  Com %.0f%% de confiança, estima-se que a verdadeira proporção\n", confianca * 100))
  cat(sprintf("  populacional está entre %s%% e %s%%.\n\n", .ic_fmt(linf*100, 1), .ic_fmt(lsup*100, 1)))
  cat(sprintf("  Segundo a interpretação frequentista, podemos dizer que se repetirmos\n"))
  cat(sprintf("  o experimento diversas vezes, %.0f%% dos intervalos de confiança\n", confianca * 100))
  cat(sprintf("  construídos a partir dessas amostras irão conter o valor verdadeiro\n"))
  cat(sprintf("  do parâmetro populacional.\n"))
  .print_rodape()
  
  if (grafico) {
    tryCatch({
      .plot_ic(p_hat, linf, lsup, 
               sprintf("Intervalo de Confiança (%.0f%%) para a Proporção", confianca * 100), 
               sprintf("Variável: %s (x = %d, n = %d)", var_label, x, n), 
               "Proporção")
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
#' @param grafico Lógico. Se TRUE, plota o intervalo do desvio padrão.
#' @return Retorna invisivelmente um data frame com os resultados.
#' @export
ic_variancia <- function(x, confianca = 0.95, decimais = 3, grafico = TRUE) {
  if (!is.numeric(x)) stop("'x' deve ser numérico.")
  x_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", x_expr)
  if (identical(var_nome, "x")) var_nome <- "x"
  
  x_c <- x[!is.na(x)]
  n   <- length(x_c)
  if (n < 2) stop("A amostra deve ter pelo menos 2 observações válidas.")
  
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
  
  .print_titulo("INTERVALO DE CONFIANÇA PARA A VARIÂNCIA")
  cat(sprintf("  Variável: %s   |   N válido: %d   |   Confiança: %.0f%%\n", var_nome, n, confianca * 100))
  cat("  Distribuição utilizada: Qui-Quadrado (χ²)\n\n")
  
  .print_topico("ESTIMATIVAS PONTUAIS")
  wc <- c(met=20, est=20)
  hdr <- paste0(.ic_pad("Métrica", wc["met"], "center"), .ic_pad("Valor", wc["est"], "center"))
  cat(.ic_sep(nchar(hdr)), "\n")
  cat("  ", hdr, "\n", sep = "")
  cat(.ic_sep(nchar(hdr)), "\n")
  cat("  ", .ic_pad("Variância", wc["met"], "center"), .ic_pad(.ic_fmt(v, decimais), wc["est"], "center"), "\n", sep = "")
  cat("  ", .ic_pad("Desvio Padrão", wc["met"], "center"), .ic_pad(.ic_fmt(dp, decimais), wc["est"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr)), "\n\n")
  
  .print_topico(sprintf("INTERVALOS DE CONFIANÇA (%.0f%%)", confianca * 100))
  w_ic <- c(met=20, int=30)
  hdr_ic <- paste0(.ic_pad("Parâmetro", w_ic["met"], "center"), .ic_pad("Intervalo [ LI ; LS ]", w_ic["int"], "center"))
  cat(.ic_sep(nchar(hdr_ic)), "\n")
  cat("  ", hdr_ic, "\n", sep = "")
  cat(.ic_sep(nchar(hdr_ic)), "\n")
  cat("  ", .ic_pad("Variância (σ²)", w_ic["met"], "center"), 
      .ic_pad(sprintf("[ %s ; %s ]", .ic_fmt(linf_v, decimais), .ic_fmt(lsup_v, decimais)), w_ic["int"], "center"), "\n", sep = "")
  cat("  ", .ic_pad("Desvio Padrão (σ)", w_ic["met"], "center"), 
      .ic_pad(sprintf("[ %s ; %s ]", .ic_fmt(linf_dp, decimais), .ic_fmt(lsup_dp, decimais)), w_ic["int"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr_ic)), "\n\n")
  
  .print_topico("INTERPRETAÇÃO")
  cat(sprintf("  Com %.0f%% de confiança, estima-se que a variabilidade da\n", confianca * 100))
  cat(sprintf("  população, medida pelo desvio padrão, está entre %s e %s.\n\n", .ic_fmt(linf_dp, decimais), .ic_fmt(lsup_dp, decimais)))
  cat(sprintf("  Segundo a interpretação frequentista, podemos dizer que se repetirmos\n"))
  cat(sprintf("  o experimento diversas vezes, %.0f%% dos intervalos de confiança\n", confianca * 100))
  cat(sprintf("  construídos a partir dessas amostras irão conter o valor verdadeiro\n"))
  cat(sprintf("  do parâmetro populacional.\n"))
  .print_rodape()
  
  if (grafico) {
    tryCatch({
      .plot_ic(dp, linf_dp, lsup_dp, 
               sprintf("Intervalo de Confiança (%.0f%%) para o Desvio Padrão", confianca * 100), 
               sprintf("Variável: %s", var_nome), "Desvio Padrão")
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }
  
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
  if (length(niveis) != 2) stop("O agrupamento deve conter exatamente 2 níveis.")
  
  res <- t.test(x_vec ~ grupo_vec, var.equal = variancia_igual, conf.level = confianca)
  
  media1 <- res$estimate[1]
  media2 <- res$estimate[2]
  dif    <- media1 - media2
  linf   <- res$conf.int[1]
  lsup   <- res$conf.int[2]
  gl     <- res$parameter
  
  tipo_txt <- if (variancia_igual) "Student (variâncias iguais)" else "Welch (variâncias não iguais)"
  
  .print_titulo("IC PARA A DIFERENÇA DE MÉDIAS")
  cat(sprintf("  Variável: %s   |   Grupos: %s vs %s\n", var_nome, niveis[1], niveis[2]))
  cat(sprintf("  Método: Teste t de %s\n\n", tipo_txt))
  
  .print_topico("ESTIMATIVAS PONTUAIS")
  wc <- c(grp=15, med=15)
  hdr_g <- paste0(.ic_pad("Grupo", wc["grp"], "center"), .ic_pad("Média", wc["med"], "center"))
  cat(.ic_sep(nchar(hdr_g)), "\n")
  cat("  ", hdr_g, "\n", sep = "")
  cat(.ic_sep(nchar(hdr_g)), "\n")
  cat("  ", .ic_pad(niveis[1], wc["grp"], "center"), .ic_pad(.ic_fmt(media1, decimais), wc["med"], "center"), "\n", sep = "")
  cat("  ", .ic_pad(niveis[2], wc["grp"], "center"), .ic_pad(.ic_fmt(media2, decimais), wc["med"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr_g)), "\n\n")
  
  .print_topico(sprintf("INTERVALO DE CONFIANÇA PARA A DIFERENÇA (%.0f%%)\n", confianca * 100))
  w_ic <- c(dif=18, int=30)
  hdr_ic <- paste0(.ic_pad("Diferença (Δ)", w_ic["dif"], "center"), .ic_pad("Intervalo [ LI ; LS ]", w_ic["int"], "center"))
  cat(.ic_sep(nchar(hdr_ic)), "\n")
  cat("  ", hdr_ic, "\n", sep = "")
  cat(.ic_sep(nchar(hdr_ic)), "\n")
  cat("  ", .ic_pad(.ic_fmt(dif, decimais), w_ic["dif"], "center"), 
      .ic_pad(sprintf("[ %s ; %s ]", .ic_fmt(linf, decimais), .ic_fmt(lsup, decimais)), w_ic["int"], "center"), "\n", sep = "")
  cat(.ic_sep(nchar(hdr_ic)), "\n\n")
  
  .print_topico("INTERPRETAÇÃO")
  if (linf * lsup > 0) {
    cat(sprintf("  Como o intervalo não contém o zero, há diferença significativa\n"))
    cat(sprintf("  entre as médias dos grupos (ao nível de %.0f%% de confiança).\n\n", confianca*100))
  } else {
    cat(sprintf("  Como o intervalo contém o zero, a diferença entre as médias\n"))
    cat(sprintf("  não é estatisticamente significativa (ao nível de %.0f%% de confiança).\n\n", confianca*100))
  }
  cat(sprintf("  Segundo a interpretação frequentista, podemos dizer que se repetirmos\n"))
  cat(sprintf("  o experimento diversas vezes, %.0f%% dos intervalos de confiança\n", confianca * 100))
  cat(sprintf("  construídos a partir dessas amostras irão conter o valor verdadeiro\n"))
  cat(sprintf("  do parâmetro populacional.\n"))
  .print_rodape()
  
  if (grafico) {
    tryCatch({
      .plot_ic(dif, linf, lsup, 
               sprintf("Intervalo de Confiança (%.0f%%) para a Diferença", confianca * 100), 
               sprintf("Variável: %s (%s vs %s)", var_nome, niveis[1], niveis[2]), 
               "Δ Média")
    }, error = function(e) message("[Aviso] Falha ao plotar: ", e$message))
  }
  
  invisible(data.frame(diferenca = dif, li = linf, ls = lsup))
}
