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
        ggplot2::labs(title = sprintf("Correla\u00e7\u00e3o: %s vs %s", v2_str, v1_str)) +
        meu_tema()
      suppressMessages(print(p_plot))
      options(warn = old_w)
    }, error = function(e) {
      message("[Aviso] N\u00e3o foi poss\u00edvel plotar o gr\u00e1fico: ", e$message)
    })
  }
  invisible(res)
}

#' @title Matriz de Correlação com Ranking dos Maiores Pares
#' @description Gera e imprime no console a matriz de correlação completa (com p-valores em asteriscos),
#'   exibe tabelas de ranking das maiores correlações positivas e negativas com classificação qualitativa
#'   centralizada, e plota o correlograma em \code{ggplot2}.
#'
#' @param dados Data frame (seleciona todas as variáveis numéricas automaticamente).
#' @param variaveis Vetor opcional com nomes específicos de variáveis para incluir.
#' @param metodo Método de correlação: \code{"pearson"} (padrão), \code{"spearman"} ou \code{"kendall"}.
#' @param top Número de pares a exibir nos rankings de maiores correlações (padrão: 5).
#'   Se \code{0} ou \code{FALSE}, desativa a exibição das tabelas de ranking.
#' @param filtrar_matriz Lógico. Se \code{TRUE}, restringe a matriz de correlação quadrada e o gráfico
#'   apenas às variáveis presentes no ranking TOP selecionado. Padrão \code{FALSE}.
#' @param criterio Critério de ordenação do ranking: \code{"maiores"} (padrão, associações mais fortes)
#'   ou \code{"menores"} (associações mais fracas / próximas de zero).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe o correlograma.
#'
#' @return Retorna invisivelmente uma lista com:
#'   \item{matriz_r}{Matriz com os coeficientes de correlação}
#'   \item{matriz_p}{Matriz com os p-valores correspondentes}
#'   \item{top_positivas}{Data frame com o ranking das correlações positivas}
#'   \item{top_negativas}{Data frame com o ranking das correlações negativas}
#' @export
matriz_correlacao <- function(dados,
                              variaveis      = NULL,
                              metodo         = c("pearson", "spearman", "kendall"),
                              top            = 5,
                              filtrar_matriz = FALSE,
                              criterio       = c("maiores", "menores"),
                              grafico        = TRUE) {
  metodo   <- match.arg(metodo)
  criterio <- match.arg(criterio)
  
  if (!is.null(variaveis)) {
    dados <- dados[, variaveis, drop = FALSE]
  }
  df_num <- dados[, sapply(dados, is.numeric), drop = FALSE]
  
  if (ncol(df_num) < 2) stop("S\u00e3o necess\u00e1rias pelo menos 2 vari\u00e1veis num\u00e9ricas para correla\u00e7\u00e3o.")
  
  # ── Cálculo da Matriz de Correlação e p-valores Completa ───────────────────
  k_total   <- ncol(df_num)
  vars_orig <- names(df_num)
  
  mat_r_completa <- stats::cor(df_num, method = metodo, use = "pairwise.complete.obs")
  mat_p_completa <- matrix(NA, k_total, k_total, dimnames = list(vars_orig, vars_orig))
  
  for (i in 1:(k_total - 1)) {
    for (j in (i + 1):k_total) {
      x_val <- df_num[[i]]; y_val <- df_num[[j]]
      idx   <- stats::complete.cases(x_val, y_val)
      if (sum(idx) > 2) {
        ct <- suppressWarnings(stats::cor.test(x_val[idx], y_val[idx], method = metodo, exact = FALSE))
        mat_p_completa[i, j] <- ct$p.value
        mat_p_completa[j, i] <- ct$p.value
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
    s  <- as.character(s)
    sp <- w - nchar(s)
    if (sp <= 0) return(s)
    if (align == "right")  return(paste0(strrep(" ", sp), s))
    if (align == "center") return(paste0(strrep(" ", floor(sp/2)), s, strrep(" ", ceiling(sp/2))))
    paste0(s, strrep(" ", sp))
  }
  
  # Classificação qualitativa da correlação
  classificar_cor <- function(val) {
    a <- abs(val)
    if (is.na(a)) return("-")
    if (a >= 0.90) return("Alt\u00edssima")
    if (a >= 0.70) return("Alta")
    if (a >= 0.40) return("Moderada")
    if (a >= 0.20) return("Fraca")
    return("Muito Fraca")
  }
  
  # ── Extração de todos os pares únicos ─────────────────────────────────────
  pares_list <- list()
  for (i in 1:(k_total - 1)) {
    for (j in (i + 1):k_total) {
      r_v <- mat_r_completa[i, j]
      p_v <- mat_p_completa[i, j]
      if (!is.na(r_v)) {
        pares_list[[length(pares_list) + 1]] <- data.frame(
          Var1          = vars_orig[i],
          Var2          = vars_orig[j],
          Par           = sprintf("%s  \u2194  %s", vars_orig[i], vars_orig[j]),
          Correlacao    = r_v,
          p_valor       = p_v,
          Significancia = asterisk(p_v),
          Classificacao = classificar_cor(r_v),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  df_todos_pares <- if (length(pares_list) > 0) do.call(rbind, pares_list) else data.frame()
  
  # ── Rankings TOP (Positivas e Negativas) ──────────────────────────────────
  top_pos <- data.frame()
  top_neg <- data.frame()
  vars_top_selecionadas <- character(0)
  
  n_top <- if (is.numeric(top) && top > 0) as.integer(top) else 0
  
  if (nrow(df_todos_pares) > 0 && n_top > 0) {
    # Positivas
    df_p <- df_todos_pares[df_todos_pares$Correlacao > 0, ]
    if (nrow(df_p) > 0) {
      ord_p <- if (criterio == "maiores") order(df_p$Correlacao, decreasing = TRUE) else order(df_p$Correlacao, decreasing = FALSE)
      df_p <- df_p[ord_p, ]
      top_pos <- head(df_p, n_top)
    }
    
    # Negativas
    df_n <- df_todos_pares[df_todos_pares$Correlacao < 0, ]
    if (nrow(df_n) > 0) {
      ord_n <- if (criterio == "maiores") order(abs(df_n$Correlacao), decreasing = TRUE) else order(abs(df_n$Correlacao), decreasing = FALSE)
      df_n <- df_n[ord_n, ]
      top_neg <- head(df_n, n_top)
    }
    
    vars_top_selecionadas <- unique(c(top_pos$Var1, top_pos$Var2, top_neg$Var1, top_neg$Var2))
  }
  
  # ── Filtragem da Matriz se solicitado ─────────────────────────────────────
  if (isTRUE(filtrar_matriz) && length(vars_top_selecionadas) >= 2) {
    df_num <- df_num[, vars_top_selecionadas, drop = FALSE]
    mat_r  <- mat_r_completa[vars_top_selecionadas, vars_top_selecionadas, drop = FALSE]
    mat_p  <- mat_p_completa[vars_top_selecionadas, vars_top_selecionadas, drop = FALSE]
  } else {
    mat_r <- mat_r_completa
    mat_p <- mat_p_completa
  }
  
  k    <- ncol(df_num)
  vars <- names(df_num)
  
  # ── Impressão no Console ──────────────────────────────────────────────────
  if (exists(".print_titulo", mode = "function")) {
    .print_titulo("MATRIZ DE CORRELA\u00c7\u00c3O")
  } else {
    cat("\n── MATRIZ DE CORRELAÇÃO ──\n")
  }
  
  cat(sprintf("  M\u00e9todo: %s\n", tools::toTitleCase(metodo)))
  cat(sprintf("  Vari\u00e1veis Analisadas: %d%s\n\n", k,
              if (isTRUE(filtrar_matriz)) " (Filtrada para as vari\u00e1veis do TOP)" else ""))
  
  vars_rotulos <- sapply(vars, function(x) substr(x, 1, 8))
  w_var <- max(nchar(vars_rotulos), 8)
  w_col <- 11
  
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
        r_txt <- formatC(mat_r[i, j], format="f", digits=2, decimal.mark=".")
        ast   <- asterisk(mat_p[i, j])
        val_txt <- paste0(r_txt, if (nchar(ast) > 0) paste0(" ", ast) else "")
      }
      linha <- paste0(linha, pad(val_txt, w_col, "center"))
    }
    cat(linha, "\n")
  }
  cat("  ", strrep("\u2500", w_total), "\n\n")
  
  # ── Impressão das Tabelas TOP (Positivas e Negativas) ─────────────────────
  if (n_top > 0) {
    w_t <- c(par = 28, cor = 16, clas = 18)
    sep_top <- paste0("  ", strrep("\u2500", sum(w_t) + length(w_t) - 1))
    
    # 1. Positivas
    if (nrow(top_pos) > 0) {
      if (exists(".print_topico", mode = "function")) {
        .print_topico(sprintf("MAIORES CORRELA\u00c7\u00d5ES POSITIVAS (TOP %d)", nrow(top_pos)))
      } else {
        cat(sprintf("\n▶ MAIORES CORRELAÇÕES POSITIVAS (TOP %d)\n", nrow(top_pos)))
      }
      cat(sep_top, "\n  ")
      cat(pad("Par de Vari\u00e1veis", w_t["par"], "center"))
      cat(pad("Correla\u00e7\u00e3o", w_t["cor"], "center"))
      cat(pad("Classifica\u00e7\u00e3o", w_t["clas"], "center"), "\n")
      cat(sep_top, "\n")
      
      for (m in seq_len(nrow(top_pos))) {
        cor_txt <- paste0(formatC(top_pos$Correlacao[m], format="f", digits=2, decimal.mark="."),
                          if (nchar(top_pos$Significancia[m]) > 0) paste0(" ", top_pos$Significancia[m]) else "")
        cat("  ")
        cat(pad(top_pos$Par[m], w_t["par"], "center"))
        cat(pad(cor_txt, w_t["cor"], "center"))
        cat(pad(top_pos$Classificacao[m], w_t["clas"], "center"), "\n")
      }
      cat(sep_top, "\n\n")
    }
    
    # 2. Negativas
    if (nrow(top_neg) > 0) {
      if (exists(".print_topico", mode = "function")) {
        .print_topico(sprintf("MAIORES CORRELA\u00c7\u00d5ES NEGATIVAS (TOP %d)", nrow(top_neg)))
      } else {
        cat(sprintf("\n▶ MAIORES CORRELAÇÕES NEGATIVAS (TOP %d)\n", nrow(top_neg)))
      }
      cat(sep_top, "\n  ")
      cat(pad("Par de Vari\u00e1veis", w_t["par"], "center"))
      cat(pad("Correla\u00e7\u00e3o", w_t["cor"], "center"))
      cat(pad("Classifica\u00e7\u00e3o", w_t["clas"], "center"), "\n")
      cat(sep_top, "\n")
      
      for (m in seq_len(nrow(top_neg))) {
        cor_txt <- paste0(formatC(top_neg$Correlacao[m], format="f", digits=2, decimal.mark="."),
                          if (nchar(top_neg$Significancia[m]) > 0) paste0(" ", top_neg$Significancia[m]) else "")
        cat("  ")
        cat(pad(top_neg$Par[m], w_t["par"], "center"))
        cat(pad(cor_txt, w_t["cor"], "center"))
        cat(pad(top_neg$Classificacao[m], w_t["clas"], "center"), "\n")
      }
      cat(sep_top, "\n\n")
    }
    
    cat("  Notas:\n")
    cat("  Signific\u00e2ncia: *** p < 0,01   ** p < 0,05   * p < 0,10\n")
    cat("  Classifica\u00e7\u00e3o (|r|): 0,90\u20131,00 (Alt\u00edssima); 0,70\u20130,89 (Alta);\n")
    cat("                       0,40\u20130,69 (Moderada); 0,20\u20130,39 (Fraca); 0,00\u20130,19 (Muito Fraca).\n\n")
  }
  
  # ── Gráfico (Correlograma) ────────────────────────────────────────────────
  if (isTRUE(grafico)) {
    if (exists("grafico_correlacao") && exists("meu_tema")) {
      tryCatch({
        legenda_titulo <- switch(metodo, pearson="R", spearman="\u03c1", kendall="\u03c4")
        p_corr <- ggplot2::ggplot() + 
                  grafico_correlacao(data = df_num, metodo = metodo, sig = TRUE, 
                                     mostrar_val = TRUE, paleta = 1) + 
                  meu_tema(estilo = 3) + 
                  ggplot2::labs(title = "Correlograma", fill = legenda_titulo)
        suppressMessages(print(p_corr))
      }, error = function(e) {
        message("[Aviso] N\u00e3o foi poss\u00edvel gerar o correlograma: ", e$message)
      })
    }
  }
  
  invisible(list(
    matriz_r      = mat_r,
    matriz_p      = mat_p,
    top_positivas = top_pos,
    top_negativas = top_neg
  ))
}

