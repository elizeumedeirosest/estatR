#' @title Teste de Correlação (Bivariada)
#' @description Realiza o teste de correlação entre duas variáveis.
#' @param dados Data frame
#' @param var1 Variável 1 (sem aspas)
#' @param var2 Variável 2 (sem aspas)
#' @param metodo Método ("pearson", "spearman", "kendall")
#' @param grafico Lógico. Se TRUE, exibe o gráfico de dispersão.
#' @export
teste_correlacao <- function(dados, var1, var2, metodo = c("pearson", "spearman", "kendall"), grafico = TRUE) {
  metodo <- match.arg(metodo)
  
  v1_expr <- rlang::enquo(var1)
  v2_expr <- rlang::enquo(var2)
  
  if (rlang::quo_is_missing(v1_expr) || rlang::quo_is_missing(v2_expr)) {
    stop("Forne\u00e7a 'var1' e 'var2'. Para a matriz completa, use matriz_correlacao(dados).")
  }
  
  v1_str <- rlang::as_name(v1_expr)
  v2_str <- rlang::as_name(v2_expr)
  
  if (!(v1_str %in% names(dados))) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", v1_str))
  if (!(v2_str %in% names(dados))) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", v2_str))
  
  x <- dados[[v1_str]]; y <- dados[[v2_str]]
  idx <- complete.cases(x, y)
  x <- x[idx]; y <- y[idx]; n <- length(x)
  
  res <- cor.test(x, y, method = metodo, exact = FALSE)
  r <- res$estimate; p <- res$p.value
  
  formata_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("< 0.001")
    return(formatC(p, format = "f", digits = 3, decimal.mark = ","))
  }
  asterisk <- function(p) {
    if (is.na(p) || p >= 0.10) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    return("*")
  }
  centraliza <- function(val, w) {
    s <- trimws(as.character(val))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  
  cat("\n\u2500\u2500 TESTE DE CORRELA\u00c7\u00c3O ", strrep("\u2500", 41), "\n")
  cat(sprintf("  Pares: %s vs %s\n", v1_str, v2_str))
  cat(sprintf("  M\u00e9todo: %s   |   N v\u00e1lido: %d\n\n", tools::toTitleCase(metodo), n))
  
  .print_topico("RESULTADOS")
  cat("  ", strrep("\u2500", 63), "\n", sep = "")
  
  if (metodo == "pearson") {
    cat(sprintf("  %-18s %s %s\n", "Estat\u00edstica", centraliza("IC (95%)", 20), centraliza("p-valor", 15)))
    ic_str <- sprintf("[%s; %s]", 
                      formatC(res$conf.int[1], format="f", digits=3, decimal.mark=","),
                      formatC(res$conf.int[2], format="f", digits=3, decimal.mark=","))
    cat(sprintf("  %-18s %s %s\n", 
                sprintf("r = %.3f", r), centraliza(ic_str, 20), centraliza(paste(formata_p(p), asterisk(p)), 15)))
  } else {
    cat(sprintf("  %-18s %s\n", "Estat\u00edstica", centraliza("p-valor", 20)))
    simbolo <- switch(metodo, spearman="rho (\u03c1)", kendall="tau (\u03c4)")
    cat(sprintf("  %-18s %s\n", 
                sprintf("%s = %.3f", simbolo, r), centraliza(paste(formata_p(p), asterisk(p)), 20)))
  }
  cat("  ", strrep("\u2500", 63), "\n", sep = "")
  if (p < 0.10) cat("  Signific\u00e2ncia: *** p < 0.01   ** p < 0.05   * p < 0.10\n")
  cat("\n")
  
  .print_topico("INTERPRETA\u00c7\u00c3O")
  forca <- if (abs(r) < 0.3) "fraca" else if (abs(r) < 0.7) "moderada" else "forte"
  direcao <- if (r > 0) "positiva" else "negativa"
  
  if (p < 0.05) {
    cat(sprintf("  Foi detectada uma associa\u00e7\u00e3o linear %s e %s\n", forca, direcao))
    cat(sprintf("  (p %s) entre '%s' e '%s'.\n\n", ifelse(p < 0.001, "< 0.001", sprintf("= %s", formata_p(p))), v1_str, v2_str))
  } else {
    cat(sprintf("  N\u00e3o h\u00e1 evid\u00eancias suficientes para afirmar que existe\n"))
    cat(sprintf("  uma correla\u00e7\u00e3o significativa entre '%s' e '%s' (p = %s).\n\n", v1_str, v2_str, formata_p(p)))
  }
  
  if (grafico && exists("grafico_de_dispersao") && exists("meu_tema")) {
    tryCatch({
      texto_leg <- switch(metodo, pearson="R", spearman="\u03c1", kendall="\u03c4")
      anotacao <- sprintf("%s = %.2f\np %s%s", texto_leg, r, 
                          ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)),
                          ifelse(nchar(asterisk(p)) > 0, paste0(" ", asterisk(p)), ""))
      
      chamada_grafico <- bquote(
        grafico_de_dispersao(data = dados, x = .(as.name(v1_str)), y = .(as.name(v2_str)),
                             reta = TRUE, ic = TRUE, outlier = FALSE, paleta = 1)
      )
      
      old_w <- getOption("warn"); options(warn = -1)
      if (r < 0) { pos_x <- Inf; pos_y <- Inf; anc_h <- 1.1; anc_v <- 1.5 } else { pos_x <- -Inf; pos_y <- Inf; anc_h <- -0.1; anc_v <- 1.5 }
      
      p_plot <- ggplot2::ggplot() + eval(chamada_grafico) +
        ggplot2::annotate("text", x = pos_x, y = pos_y, label = anotacao,
                          hjust = anc_h, vjust = anc_v, size = 6.5, fontface = "bold", color = "#333333") +
        ggplot2::labs(title = sprintf("Correla\u00e7\u00e3o: %s vs %s", v2_str, v1_str), caption = "estatR") +
        meu_tema()
      suppressMessages(print(p_plot))
      options(warn = old_w)
    }, error = function(e) {
      message("[Aviso] N\u00e3o foi poss\u00edvel plotar o gr\u00e1fico: ", e$message)
    })
  }
  invisible(res)
}

#' @title Matriz de Correlação (Tabela e Correlograma)
#' @description Gera e imprime no console a matriz de correlação completa (com p-valores em asteriscos)
#'   e plota o correlograma.
#' @param dados Data frame (usa todas as variáveis numéricas automaticamente)
#' @param variaveis Vetor opcional com nomes específicos de variáveis para incluir.
#' @param metodo Método ("pearson", "spearman", "kendall")
#' @param top Opcional. Se definido (ex: 5), exibe no gráfico apenas as variáveis envolvidas
#'   nos 'top' pares com maior correlação absoluta.
#' @export
matriz_correlacao <- function(dados, variaveis = NULL, metodo = c("pearson", "spearman", "kendall"), top = NULL) {
  metodo <- match.arg(metodo)
  
  if (!is.null(variaveis)) {
    dados <- dados[, variaveis, drop = FALSE]
  }
  df_num <- dados[, sapply(dados, is.numeric), drop = FALSE]
  
  if (ncol(df_num) < 2) stop("S\u00e3o necess\u00e1rias pelo menos 2 vari\u00e1veis num\u00e9ricas para correla\u00e7\u00e3o.")
  
  # ── Cálculo da Matriz de Correlação e p-valores ───────────────────────────
  k <- ncol(df_num)
  vars <- names(df_num)
  
  mat_r <- cor(df_num, method = metodo, use = "pairwise.complete.obs")
  mat_p <- matrix(NA, k, k)
  
  for (i in 1:(k-1)) {
    for (j in (i+1):k) {
      x_val <- df_num[[i]]; y_val <- df_num[[j]]
      idx <- complete.cases(x_val, y_val)
      if (sum(idx) > 2) {
        ct <- suppressWarnings(cor.test(x_val[idx], y_val[idx], method = metodo, exact = FALSE))
        mat_p[i, j] <- ct$p.value
        mat_p[j, i] <- ct$p.value
      }
    }
  }
  
  asterisk <- function(p) {
    if (is.na(p) || p >= 0.10) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    return("*")
  }
  
  pad <- function(s, w, align = "left") {
    s <- as.character(s); sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right")  return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }
  
  # ── Impressão no Console (Tabela Arrumada) ────────────────────────────────
  .print_titulo("MATRIZ DE CORRELA\u00c7\u00c3O %s")
  cat(sprintf("  M\u00e9todo: %s\n\n", tools::toTitleCase(metodo)))
  
  # Preparando a tabela de texto
  # Para não quebrar a tela, se houver muitas variáveis, avisamos e truncamos ou imprimimos tudo.
  # Geralmente, 10 variáveis cabem bem se usarmos 8 caracteres por coluna.
  vars_rotulos <- sapply(vars, function(x) substr(x, 1, 8))
  w_var <- max(nchar(vars_rotulos))
  w_col <- 10 # Largura das colunas numéricas
  
  hdr_linha <- paste0("  ", pad("Vari\u00e1vel", w_var, "left"), " |")
  for (v in vars_rotulos) {
    hdr_linha <- paste0(hdr_linha, pad(v, w_col, "center"))
  }
  
  w_total <- nchar(hdr_linha) - 2
  cat("  ", strrep("\u2500", w_total), "\n", sep = "")
  cat(hdr_linha, "\n")
  cat("  ", strrep("\u2500", w_total), "\n", sep = "")
  
  for (i in 1:k) {
    linha <- paste0("  ", pad(vars_rotulos[i], w_var, "left"), " |")
    for (j in 1:k) {
      if (i == j) {
        val_txt <- "1.00"
      } else {
        r_txt <- formatC(mat_r[i, j], format="f", digits=2, decimal.mark=".") # Usando ponto para alinhar visualmente
        ast <- asterisk(mat_p[i, j])
        val_txt <- paste0(r_txt, ast)
      }
      linha <- paste0(linha, pad(val_txt, w_col, "center"))
    }
    cat(linha, "\n")
  }
  cat("  ", strrep("\u2500", w_total), "\n", sep = "")
  cat("  Signific\u00e2ncia: *** p < 0.01   ** p < 0.05   * p < 0.10\n\n")
  
  # ── Gráfico (Correlograma com filtro top opcional) ────────────────────────
  if (!is.null(top) && is.numeric(top) && top > 0) {
    mat_tmp <- mat_r
    mat_tmp[lower.tri(mat_tmp, diag = TRUE)] <- NA
    
    df_pairs <- as.data.frame(as.table(mat_tmp))
    df_pairs <- df_pairs[!is.na(df_pairs$Freq), ]
    df_pairs <- df_pairs[order(abs(df_pairs$Freq), decreasing = TRUE), ]
    
    n_pares <- min(top, nrow(df_pairs))
    pares_selecionados <- df_pairs[1:n_pares, ]
    
    vars_top <- unique(c(as.character(pares_selecionados$Var1), as.character(pares_selecionados$Var2)))
    df_num <- df_num[, vars_top, drop = FALSE]
    
    message(sprintf("[estatR] Gr\u00e1fico reduzido: exibindo as %d vari\u00e1veis mais fortemente correlacionadas.", length(vars_top)))
  }
  
  if (exists("grafico_correlacao") && exists("meu_tema")) {
    tryCatch({
      legenda_titulo <- switch(metodo, pearson="R", spearman="\u03c1", kendall="\u03c4")
      
      p_corr <- ggplot2::ggplot() + 
                grafico_correlacao(data = df_num, metodo = metodo, sig = TRUE, 
                                   mostrar_val = TRUE, paleta = 1) + 
                meu_tema(estilo = 3) + 
                ggplot2::labs(title = "Correlograma", caption = "estatR", fill = legenda_titulo)
      
      suppressMessages(print(p_corr))
    }, error = function(e) {
      message("[Aviso] N\u00e3o foi poss\u00edvel gerar o correlograma: ", e$message)
    })
  } else {
    message("[Aviso] A fun\u00e7\u00e3o 'grafico_correlacao' n\u00e3o foi encontrada.")
  }
  
  invisible(mat_r)
}
