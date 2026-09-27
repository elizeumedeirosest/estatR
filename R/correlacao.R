#' @title Teste de Correlação (Bivariada)
#' @description Realiza o teste de correlação entre duas variáveis contínuas.
#' @param dados Data frame
#' @param var1 Variável 1 (sem aspas)
#' @param var2 Variável 2 (sem aspas)
#' @param metodo Método ("pearson", "spearman", "kendall")
#' @param grafico Lógico. Se TRUE, exibe o gráfico de dispersão.
#' @export
teste_correlacao <- function(dados, var1, var2, metodo = c("pearson", "spearman", "kendall"), grafico = TRUE) {
  metodo <- match.arg(metodo)
  
  v1_expr <- rlang::enexpr(var1)
  v2_expr <- rlang::enexpr(var2)
  v1_str  <- rlang::as_name(v1_expr)
  v2_str  <- rlang::as_name(v2_expr)
  
  if (!(v1_str %in% names(dados))) stop(sprintf("Variável '%s' não encontrada.", v1_str))
  if (!(v2_str %in% names(dados))) stop(sprintf("Variável '%s' não encontrada.", v2_str))
  
  x <- dados[[v1_str]]
  y <- dados[[v2_str]]
  
  # Remove NAs par-a-par
  idx <- complete.cases(x, y)
  x <- x[idx]
  y <- y[idx]
  n <- length(x)
  
  res <- cor.test(x, y, method = metodo, exact = FALSE)
  
  r <- res$estimate
  p <- res$p.value
  
  formata_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("<0,001")
    return(formatC(p, format = "f", digits = 3, decimal.mark = ","))
  }
  
  centraliza <- function(val, w) {
    s <- trimws(as.character(val))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  
  linha <- paste0("\n", strrep("\u2500", 54), "\n\n")
  
  cat("\nTESTE DE CORRELAÇÃO\n\n")
  cat(sprintf("Variável 1: %s\n", v1_str))
  cat(sprintf("Variável 2: %s\n", v2_str))
  cat(sprintf("Método:     %s\n", tools::toTitleCase(metodo)))
  cat(sprintf("N válido:   %d\n", n))
  
  cat(linha)
  cat("RESULTADOS\n\n")
  
  if (metodo == "pearson") {
    cat(sprintf("%-18s %s %s\n", "Estatística", centraliza("IC (95%)", 20), centraliza("p-valor", 10)))
    ic_str <- sprintf("[%.3f ; %.3f]", res$conf.int[1], res$conf.int[2])
    cat(sprintf("%-18s %s %s\n", 
                sprintf("%s = %.3f", "r", r),
                centraliza(ic_str, 20),
                centraliza(formata_p(p), 10)))
  } else {
    cat(sprintf("%-18s %s\n", "Estatística", centraliza("p-valor", 10)))
    cat(sprintf("%-18s %s\n", 
                sprintf("%s = %.3f", switch(metodo, spearman="rho", kendall="tau"), r),
                centraliza(formata_p(p), 10)))
  }
  
  cat(linha)
  cat("INTERPRETAÇÃO\n\n")
  
  forca <- if (abs(r) < 0.3) "fraca" else if (abs(r) < 0.7) "moderada" else "forte"
  direcao <- if (r > 0) "positiva" else "negativa"
  
  if (p < 0.05) {
    cat(sprintf("Foi detectada uma associação linear %s e %s\n(p %s) entre '%s' e '%s'.\n\n",
                forca, direcao, ifelse(p < 0.001, "< 0,001", sprintf("= %s", formata_p(p))), v1_str, v2_str))
  } else {
    cat(sprintf("Não há evidências suficientes para afirmar que existe\numa correlação significativa entre '%s' e '%s' (p = %s).\n\n",
                v1_str, v2_str, formata_p(p)))
  }
  
  if (grafico) {
    if (exists("grafico_de_dispersao") && exists("meu_tema")) {
      texto_legenda <- switch(metodo, pearson="R", spearman="\u03c1", kendall="\u03c4")
      anotacao <- sprintf("%s = %.2f\np %s", texto_legenda, r, ifelse(p < 0.001, "< 0,001", sprintf("= %.3f", p)))
      
      chamada_grafico <- bquote(
        grafico_de_dispersao(data = dados, x = .(as.name(v1_str)), y = .(as.name(v2_str)),
                             reta = TRUE, ic = TRUE, outlier = FALSE, paleta = 1)
      )
      
      tryCatch({
        old_w <- getOption("warn")
        options(warn = -1)
        
        # Posição inteligente para não trombar com os pontos (igual na regressão)
        if (r < 0) {
          pos_x <- Inf; pos_y <- Inf; anc_h <- 1.1; anc_v <- 1.5
        } else {
          pos_x <- -Inf; pos_y <- Inf; anc_h <- -0.1; anc_v <- 1.5
        }
        
        p_plot <- ggplot2::ggplot() + eval(chamada_grafico) +
          ggplot2::annotate("text", x = pos_x, y = pos_y, label = anotacao,
                            hjust = anc_h, vjust = anc_v, size = 6.5, fontface = "bold", color = "#333333") +
          ggplot2::labs(title = sprintf("Correlação: %s vs %s", v2_str, v1_str),
                        caption = "estatR") +
          meu_tema(grade = "dupla")
        
        suppressMessages(print(p_plot))
        options(warn = old_w)
      }, error = function(e) {
        message("[Aviso] Não foi possível plotar o gráfico: ", e$message)
      })
    }
  }
  
  invisible(res)
}

#' @title Matriz de Correlação
#' @description Plota o correlograma de todas as variáveis numéricas do data frame.
#' @param dados Data frame (usa todas as variáveis numéricas automaticamente)
#' @param variaveis Vetor opcional com nomes das variáveis. Se NULL, usa todas numéricas.
#' @param metodo Método ("pearson", "spearman", "kendall")
#' @export
matriz_correlacao <- function(dados, variaveis = NULL, metodo = c("pearson", "spearman", "kendall")) {
  metodo <- match.arg(metodo)
  
  if (!is.null(variaveis)) {
    dados <- dados[, variaveis, drop = FALSE]
  }
  
  df_num <- dados[, sapply(dados, is.numeric), drop = FALSE]
  
  if (ncol(df_num) < 2) {
    stop("São necessárias pelo menos 2 variáveis numéricas para correlação.")
  }
  
  mat_r <- cor(df_num, method = metodo, use = "pairwise.complete.obs")
  
  if (exists("grafico_correlacao") && exists("meu_tema")) {
    tryCatch({
      legenda_titulo <- switch(metodo, pearson="R", spearman="\u03c1", kendall="\u03c4")
      
      p_corr <- ggplot2::ggplot() + 
                grafico_correlacao(data = df_num, 
                                   metodo = metodo,
                                   sig = TRUE, 
                                   mostrar_val = TRUE, 
                                   paleta = 1) + 
                meu_tema(estilo = 3) + 
                ggplot2::labs(title = "Correlograma", caption = "estatR", fill = legenda_titulo)
      
      suppressMessages(print(p_corr))
    }, error = function(e) {
      message("[Aviso] Não foi possível gerar o correlograma: ", e$message)
    })
  } else {
    message("[Aviso] A função 'grafico_correlacao' não foi encontrada. Verifique se metaR está carregado.")
  }
  
  invisible(mat_r)
}
