#' @title Teste de Correlação (Bivariada ou Multivariada)
#' @description Realiza o teste de correlação. Se duas variáveis forem fornecidas, faz o teste
#'   detalhado entre elas. Se nenhuma variável for fornecida, calcula a correlação entre todas
#'   as variáveis numéricas do banco de dados e exibe os pares mais correlacionados.
#' @param dados Data frame
#' @param var1 Variável 1 (opcional, sem aspas)
#' @param var2 Variável 2 (opcional, sem aspas)
#' @param metodo Método ("pearson", "spearman", "kendall")
#' @param top Quando \code{var1} e \code{var2} não são fornecidos, especifica o número de pares
#'   mais correlacionados a serem exibidos no console (padrão 10).
#' @param grafico Lógico. Se TRUE e for um teste bivariado, exibe o gráfico de dispersão.
#' @export
teste_correlacao <- function(dados, var1 = NULL, var2 = NULL, metodo = c("pearson", "spearman", "kendall"), 
                             top = 10, grafico = TRUE) {
  metodo <- match.arg(metodo)
  
  v1_expr <- rlang::enquo(var1)
  v2_expr <- rlang::enquo(var2)
  
  v1_nulo <- rlang::quo_is_null(v1_expr)
  v2_nulo <- rlang::quo_is_null(v2_expr)
  
  # ── Helper de formatação ────────────────────────────────────────────────
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
  
  linha <- paste0("\n", strrep("\u2500", 65), "\n\n")
  
  # ── ROTA MULTIVARIADA (Toda a base) ───────────────────────────────────────
  if (v1_nulo && v2_nulo) {
    df_num <- dados[, sapply(dados, is.numeric), drop = FALSE]
    if (ncol(df_num) < 2) stop("O banco precisa ter pelo menos 2 vari\u00e1veis num\u00e9ricas.")
    
    vars <- names(df_num)
    k <- length(vars)
    res_list <- list()
    
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
        x_val <- df_num[[i]]
        y_val <- df_num[[j]]
        idx <- complete.cases(x_val, y_val)
        
        if (sum(idx) > 2) {
          ct <- suppressWarnings(cor.test(x_val[idx], y_val[idx], method = metodo, exact = FALSE))
          res_list[[length(res_list) + 1]] <- data.frame(
            Var1 = vars[i], Var2 = vars[j],
            r = ct$estimate, p = ct$p.value,
            row.names = NULL, stringsAsFactors = FALSE
          )
        }
      }
    }
    
    res_df <- do.call(rbind, res_list)
    res_df <- res_df[order(abs(res_df$r), decreasing = TRUE), ]
    
    n_show <- min(top, nrow(res_df))
    df_show <- res_df[1:n_show, ]
    
    cat("\n\u2500\u2500 MATRIZ DE CORRELA\u00c7\u00c3O (RANKING) ", strrep("\u2500", 33), "\n")
    cat(sprintf("  M\u00e9todo: %s   |   Total de pares testados: %d\n\n", tools::toTitleCase(metodo), nrow(res_df)))
    cat(sprintf("  \u25b6 TOP %d PARES MAIS CORRELACIONADOS\n", n_show))
    
    hdr <- sprintf("  %-15s %-15s %-12s %s", "Vari\u00e1vel 1", "Vari\u00e1vel 2", "Correla\u00e7\u00e3o", "p-valor")
    cat("  ", strrep("\u2500", 63), "\n", sep = "")
    cat(hdr, "\n")
    cat("  ", strrep("\u2500", 63), "\n", sep = "")
    
    for (i in 1:nrow(df_show)) {
      v1_txt <- substr(df_show$Var1[i], 1, 14)
      v2_txt <- substr(df_show$Var2[i], 1, 14)
      r_txt  <- formatC(df_show$r[i], format="f", digits=3, decimal.mark=",")
      p_txt  <- paste(formata_p(df_show$p[i]), asterisk(df_show$p[i]))
      
      cat(sprintf("  %-15s %-15s %-12s %s\n", v1_txt, v2_txt, r_txt, p_txt))
    }
    cat("  ", strrep("\u2500", 63), "\n", sep = "")
    cat("  Signific\u00e2ncia: *** p < 0.01   ** p < 0.05   * p < 0.10\n\n")
    
    return(invisible(res_df))
  }
  
  # ── ROTA BIVARIADA (Duas variáveis) ───────────────────────────────────────
  v1_str <- rlang::as_name(v1_expr)
  v2_str <- rlang::as_name(v2_expr)
  
  if (!(v1_str %in% names(dados))) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", v1_str))
  if (!(v2_str %in% names(dados))) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", v2_str))
  
  x <- dados[[v1_str]]; y <- dados[[v2_str]]
  idx <- complete.cases(x, y)
  x <- x[idx]; y <- y[idx]; n <- length(x)
  
  res <- cor.test(x, y, method = metodo, exact = FALSE)
  r <- res$estimate; p <- res$p.value
  
  centraliza <- function(val, w) {
    s <- trimws(as.character(val))
    pad <- w - nchar(s)
    if (pad <= 0) return(s)
    paste0(strrep(" ", floor(pad/2)), s, strrep(" ", ceiling(pad/2)))
  }
  
  cat("\n\u2500\u2500 TESTE DE CORRELA\u00c7\u00c3O ", strrep("\u2500", 41), "\n")
  cat(sprintf("  Pares: %s vs %s\n", v1_str, v2_str))
  cat(sprintf("  M\u00e9todo: %s   |   N v\u00e1lido: %d\n\n", tools::toTitleCase(metodo), n))
  
  cat("  \u25b6 RESULTADOS\n")
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
  
  cat("  \u25b6 INTERPRETA\u00c7\u00c3O\n")
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

#' @title Matriz de Correlação (Correlograma)
#' @description Plota o correlograma de todas as variáveis numéricas do data frame.
#'   Permite filtrar pelas variáveis mais correlacionadas para simplificar visualizações densas.
#' @param dados Data frame (usa todas as variáveis numéricas automaticamente)
#' @param variaveis Vetor opcional com nomes específicos de variáveis para incluir.
#' @param metodo Método ("pearson", "spearman", "kendall")
#' @param top Opcional. Se definido (ex: 5), exibe no gráfico apenas as variáveis envolvidas
#'   nos 'top' pares com maior correlação absoluta, ajudando a limpar gráficos muito grandes.
#' @export
matriz_correlacao <- function(dados, variaveis = NULL, metodo = c("pearson", "spearman", "kendall"), top = NULL) {
  metodo <- match.arg(metodo)
  
  if (!is.null(variaveis)) {
    dados <- dados[, variaveis, drop = FALSE]
  }
  df_num <- dados[, sapply(dados, is.numeric), drop = FALSE]
  
  if (ncol(df_num) < 2) stop("S\u00e3o necess\u00e1rias pelo menos 2 vari\u00e1veis num\u00e9ricas para correla\u00e7\u00e3o.")
  
  # Filtragem do Top N mais correlacionados
  if (!is.null(top) && is.numeric(top) && top > 0) {
    mat_tmp <- cor(df_num, method = metodo, use = "pairwise.complete.obs")
    mat_tmp[lower.tri(mat_tmp, diag = TRUE)] <- NA
    
    df_pairs <- as.data.frame(as.table(mat_tmp))
    df_pairs <- df_pairs[!is.na(df_pairs$Freq), ]
    df_pairs <- df_pairs[order(abs(df_pairs$Freq), decreasing = TRUE), ]
    
    n_pares <- min(top, nrow(df_pairs))
    pares_selecionados <- df_pairs[1:n_pares, ]
    
    vars_top <- unique(c(as.character(pares_selecionados$Var1), as.character(pares_selecionados$Var2)))
    df_num <- df_num[, vars_top, drop = FALSE]
    
    message(sprintf("[estatR] Exibindo matriz reduzida: %d vari\u00e1veis mais fortemente correlacionadas.", length(vars_top)))
  }
  
  mat_r <- cor(df_num, method = metodo, use = "pairwise.complete.obs")
  
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
