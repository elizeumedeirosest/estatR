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




# ── COEFICIENTE DE CORRELAÇÃO INTRACLASSE (ICC) ──────────────────────────────

#' @title Coeficiente de Correlação Intraclasse (ICC)
#' @description Calcula os Coeficientes de Correlação Intraclasse (ICC) para avaliar a
#'   concordância, reprodutibilidade e confiabilidade entre múltiplos avaliadores, juízes
#'   ou métodos de mensuração contínua, conforme Shrout & Fleiss (1979) e McGraw & Wong (1996).
#'
#' @param dados Matriz ou \code{data.frame} em formato largo onde cada coluna representa um
#'   avaliador/juiz e cada linha representa um sujeito/objeto avaliado.
#' @param tipo Tipo de modelo de ICC a reportar: \code{"todos"} (padrão, exibe os 6 modelos clássicos),
#'   \code{"icc1"} (1 fator aleatório), \code{"icc2"} (2 fatores aleatórios / concordância absoluta)
#'   ou \code{"icc3"} (2 fatores mistos / consistência).
#' @param confianca Nível de confiança para os intervalos (padrão 0.95).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe gráfico de perfis das medições por avaliador.
#'
#' @return Objeto da classe \code{estatR_icc} contendo:
#' \itemize{
#'   \item \code{tabela}: Tabela completa com os modelos, estimativas pontuais, IC 95%, estatística F e p-valor.
#'   \item \code{anova}: Quadro de ANOVA da decomposição da variância (entre sujeitos, entre juízes e erro).
#'   \item \code{grafico}: Objeto gráfico do ggplot2.
#' }
#' @export
coeficiente_cor_intra_classe <- function(dados,
                                         tipo      = c("todos", "icc1", "icc2", "icc3"),
                                         confianca = 0.95,
                                         grafico   = TRUE) {

  tipo <- match.arg(tipo)

  mat <- as.matrix(dados)
  if (!is.numeric(mat)) stop("Os 'dados' devem ser numéricos (sujeitos nas linhas, juízes nas colunas).", call. = FALSE)

  # Remover linhas com valores ausentes
  idx_completos <- stats::complete.cases(mat)
  mat <- mat[idx_completos, , drop = FALSE]

  n <- nrow(mat) # número de sujeitos
  k <- ncol(mat) # número de avaliadores

  if (n < 3) stop("São necessários pelo menos 3 sujeitos avaliados.", call. = FALSE)
  if (k < 2) stop("São necessários pelo menos 2 avaliadores (colunas).", call. = FALSE)

  # ── Decomposição ANOVA Two-Way sem repetição ─────────────────────────────────
  grand_mean <- mean(mat)
  row_means  <- rowMeans(mat)
  col_means  <- colMeans(mat)

  ss_total <- sum((mat - grand_mean)^2)
  ss_rows  <- k * sum((row_means - grand_mean)^2) # Entre sujeitos
  ss_cols  <- n * sum((col_means - grand_mean)^2) # Entre avaliadores
  ss_error <- ss_total - ss_rows - ss_cols        # Resíduo / Interação
  ss_within <- ss_total - ss_rows                 # Dentro dos sujeitos

  df_rows  <- n - 1
  df_cols  <- k - 1
  df_error <- (n - 1) * (k - 1)
  df_within <- n * (k - 1)

  ms_r <- ss_rows / df_rows
  ms_c <- ss_cols / df_cols
  ms_e <- max(ss_error / df_error, 1e-12)
  ms_w <- ss_within / df_within

  # ── Cálculo dos 6 Modelos de ICC (Shrout & Fleiss, 1979) ──────────────────────
  # ICC(1,1) e ICC(1,k) - One-way random
  icc1_1 <- (ms_r - ms_w) / (ms_r + (k - 1) * ms_w)
  icc1_k <- (ms_r - ms_w) / ms_r
  f1 <- ms_r / ms_w
  p1 <- 1 - stats::pf(f1, df_rows, df_within)

  f1_l <- f1 / stats::qf(1 - (1 - confianca)/2, df_rows, df_within)
  f1_u <- f1 * stats::qf(1 - (1 - confianca)/2, df_within, df_rows)
  ic_icc1_1 <- c((f1_l - 1)/(f1_l + k - 1), (f1_u - 1)/(f1_u + k - 1))
  ic_icc1_k <- c(1 - 1/f1_l, 1 - 1/f1_u)

  # ICC(2,1) e ICC(2,k) - Two-way random (Concordância Absoluta)
  icc2_1 <- (ms_r - ms_e) / (ms_r + (k - 1) * ms_e + (k / n) * (ms_c - ms_e))
  icc2_k <- (ms_r - ms_e) / (ms_r + (ms_c - ms_e) / n)
  f2 <- ms_r / ms_e
  p2 <- 1 - stats::pf(f2, df_rows, df_error)

  # Aproximação exata F de McGraw & Wong para ICC2
  f2_l <- f2 / stats::qf(1 - (1 - confianca)/2, df_rows, df_error)
  f2_u <- f2 * stats::qf(1 - (1 - confianca)/2, df_error, df_rows)
  ic_icc2_1 <- c(max(-1, (f2_l - 1)/(f2_l + k - 1)), min(1, (f2_u - 1)/(f2_u + k - 1)))
  ic_icc2_k <- c(max(-1, 1 - 1/f2_l), min(1, 1 - 1/f2_u))

  # ICC(3,1) e ICC(3,k) - Two-way mixed (Consistência)
  icc3_1 <- (ms_r - ms_e) / (ms_r + (k - 1) * ms_e)
  icc3_k <- (ms_r - ms_e) / ms_r
  f3 <- ms_r / ms_e
  p3 <- p2
  ic_icc3_1 <- ic_icc2_1
  ic_icc3_k <- ic_icc2_k

  # Classificação qualitativa (Koo & Li, 2016)
  classif_fn <- function(val) {
    if (is.na(val)) return("-")
    if (val < 0.50) "Pobre"
    else if (val < 0.75) "Moderada"
    else if (val < 0.90) "Boa"
    else "Excelente"
  }

  tb_icc <- data.frame(
    Modelo = c("ICC(1,1)", "ICC(1,k)", "ICC(2,1)", "ICC(2,k)", "ICC(3,1)", "ICC(3,k)"),
    Tipo_Avaliacao = c("Individual", "Média Juízes", "Individual", "Média Juízes", "Individual", "Média Juízes"),
    Natureza = c("1 Fator Aleatório", "1 Fator Aleatório",
                 "Concordância Absoluta", "Concordância Absoluta",
                 "Consistência Relativa", "Consistência Relativa"),
    Estimativa = c(icc1_1, icc1_k, icc2_1, icc2_k, icc3_1, icc3_k),
    IC_Inf = c(ic_icc1_1[1], ic_icc1_k[1], ic_icc2_1[1], ic_icc2_k[1], ic_icc3_1[1], ic_icc3_k[1]),
    IC_Sup = c(ic_icc1_1[2], ic_icc1_k[2], ic_icc2_1[2], ic_icc2_k[2], ic_icc3_1[2], ic_icc3_k[2]),
    F_stat = c(f1, f1, f2, f2, f3, f3),
    p_valor = c(p1, p1, p2, p2, p3, p3),
    Classificacao = sapply(c(icc1_1, icc1_k, icc2_1, icc2_k, icc3_1, icc3_k), classif_fn),
    stringsAsFactors = FALSE
  )

  if (tipo != "todos") {
    pref <- toupper(tipo)
    tb_icc <- tb_icc[grepl(pref, tb_icc$Modelo, ignore.case = TRUE), ]
  }

  # ── Saída no console ─────────────────────────────────────────────────────────
  .print_titulo("COEFICIENTE DE CORRELAÇÃO INTRACLASSE (ICC)")

  cat(sprintf("  Sujeitos Avaliados (n): %d  |  Avaliadores/Juízes (k): %d
", n, k))
  cat(sprintf("  Nível de Confiança: %.0f%%

", confianca * 100))

  .print_topico("TABELA DE COEFICIENTES DE CONFIABILIDADE")
  w_icc <- c(10, 14, 22, 10, 18, 12, 14)
  cat(.th_sep(sum(w_icc)), "
")
  cat("  ",
      .th_pad("Modelo", w_icc[1], "left"),
      .th_pad("Medida", w_icc[2], "center"),
      .th_pad("Natureza", w_icc[3], "left"),
      .th_pad("ICC", w_icc[4], "center"),
      .th_pad(sprintf("IC (%.0f%%)", confianca*100), w_icc[5], "center"),
      .th_pad("p-valor", w_icc[6], "center"),
      .th_pad("Confiabilidade", w_icc[7], "center"),
      "
", sep = "")
  cat(.th_sep(sum(w_icc)), "
")

  for (i in 1:nrow(tb_icc)) {
    ic_str <- sprintf("[%.2f; %.2f]", tb_icc$IC_Inf[i], tb_icc$IC_Sup[i])
    cat("  ",
        .th_pad(tb_icc$Modelo[i], w_icc[1], "left"),
        .th_pad(tb_icc$Tipo_Avaliacao[i], w_icc[2], "center"),
        .th_pad(tb_icc$Natureza[i], w_icc[3], "left"),
        .th_pad(formatC(tb_icc$Estimativa[i], format = "f", digits = 3, decimal.mark = ","), w_icc[4], "center"),
        .th_pad(ic_str, w_icc[5], "center"),
        .th_pad(.fmt_p(tb_icc$p_valor[i]), w_icc[6], "center"),
        .th_pad(tb_icc$Classificacao[i], w_icc[7], "center"),
        "
", sep = "")
  }
  cat(.th_sep(sum(w_icc)), "

")

  # Nota explicativa
  cat("  * Critério de Koo & Li (2016): <0,50 Pobre | 0,50-0,75 Moderada | 0,75-0,90 Boa | >=0,90 Excelente.
")
  cat("  * Para avaliar se juízes atribuem notas idênticas, recomenda-se ICC(2,1) [Concordância Absoluta].
")
  .print_rodape()

  p_icc <- NULL
  if (grafico) {
    nomes_cols <- colnames(mat)
    if (is.null(nomes_cols)) nomes_cols <- paste0("Juiz_", 1:k)
    df_long <- data.frame(
      Sujeito   = rep(1:n, times = k),
      Avaliador = factor(rep(nomes_cols, each = n), levels = nomes_cols),
      Escore    = as.vector(mat)
    )

    p_icc <- ggplot2::ggplot(df_long, ggplot2::aes(x = Sujeito, y = Escore, color = Avaliador, group = Avaliador)) +
      ggplot2::geom_line(linewidth = 0.9, alpha = 0.8) +
      ggplot2::geom_point(size = 2.2) +
      ggplot2::scale_x_continuous(breaks = seq(1, n, by = max(1, floor(n/10)))) +
      ggplot2::labs(
        title    = "Perfis das Medições por Avaliador",
        subtitle = sprintf("Concordância Geral ICC(2,1) = %.3f (%s)", icc2_1, classif_fn(icc2_1)),
        x = "Sujeito / Caso Avaliado", y = "Escore Atribuído"
      )

    if (exists("meu_tema", mode = "function")) {
      p_icc <- p_icc + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p_icc <- p_icc + tema_estatR(estilo = 1, grade = "dupla")
    } else {
      p_icc <- p_icc + ggplot2::theme_minimal()
    }
    suppressMessages(print(p_icc))
  }

  res <- list(
    tabela  = tb_icc,
    anova   = data.frame(
      Fonte = c("Entre Sujeitos", "Entre Juízes", "Resíduo / Erro"),
      SQ    = c(ss_rows, ss_cols, ss_error),
      GL    = c(df_rows, df_cols, df_error),
      QM    = c(ms_r, ms_c, ms_e)
    ),
    grafico = p_icc
  )
  class(res) <- c("estatR_icc", "list")
  return(invisible(res))
}
