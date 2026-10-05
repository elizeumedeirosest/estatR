# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: DIAGNÓSTICO
# Funções para inspecionar estrutura e valores ausentes do banco de dados
# ─────────────────────────────────────────────────────────────────────────────

# ── MOTOR DE TABELAS INTELIGENTE ──────────────────────────────────────────────
.print_tabela_estatR <- function(df, align = NULL) {
  if (nrow(df) == 0) return(invisible(NULL))
  
  nomes_cols <- names(df)
  n_cols <- length(nomes_cols)
  
  if (is.null(align)) {
    align <- rep("center", n_cols)
    align[1] <- "left" # Por padrao, a primeira coluna (Variavel) alinha a esquerda
  }
  
  # Calcula largura de cada coluna baseada no cabecalho e nos dados
  larguras <- sapply(1:n_cols, function(i) {
    max(nchar(as.character(nomes_cols[i])), max(nchar(as.character(df[[i]])), na.rm = TRUE))
  })
  
  # Padding extra
  pad <- 3
  larguras <- larguras + pad
  
  # Helper de alinhamento
  .pad_str <- function(s, w, al) {
    s <- as.character(s)
    sp <- max(0, w - nchar(s))
    if (al == "left") return(paste0(" ", s, strrep(" ", sp - 1)))
    if (al == "right") return(paste0(strrep(" ", sp - 1), s, " "))
    # Center
    left_pad <- floor(sp / 2)
    right_pad <- ceiling(sp / 2)
    paste0(strrep(" ", left_pad), s, strrep(" ", right_pad))
  }
  
  # Construir cabecalho
  hdr <- paste(sapply(1:n_cols, function(i) .pad_str(nomes_cols[i], larguras[i], align[i])), collapse = "")
  
  cat("  ", strrep("\u2500", nchar(hdr)), "\n", sep = "")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", strrep("\u2500", nchar(hdr)), "\n", sep = "")
  
  # Construir linhas
  for (r in 1:nrow(df)) {
    linha <- paste(sapply(1:n_cols, function(i) .pad_str(df[r, i], larguras[i], align[i])), collapse = "")
    cat("  ", linha, "\n", sep = "")
  }
  cat("  ", strrep("\u2500", nchar(hdr)), "\n\n", sep = "")
}

# ── FUNÇÕES AUXILIARES ────────────────────────────────────────────────────────
.tipo_traduzido <- function(vec) {
  if (is.integer(vec)) return("Integer (Inteiro)")
  if (is.numeric(vec)) return("Double (Decimal)")
  if (is.factor(vec))  return("Factor (Fator)")
  if (is.character(vec)) return("Character (Texto)")
  if (is.logical(vec)) return("Logical (Lógico)")
  if (inherits(vec, c("Date", "POSIXt"))) return("Date (Data)")
  return(class(vec)[1])
}

.truncar_niveis <- function(vec, max_char = 30) {
  if (!is.factor(vec) && !is.character(vec)) return("-")
  lvls <- if (is.factor(vec)) levels(vec) else unique(na.omit(vec))
  if (length(lvls) == 0) return("-")
  
  txt <- paste(lvls, collapse = ", ")
  if (nchar(txt) > max_char) {
    txt <- paste0(substr(txt, 1, max_char - 3), "...")
  }
  return(txt)
}

.detectar_suspeitos <- function(vec) {
  if (is.numeric(vec)) {
    suspeitos <- c(-99, -999, -9999, 99, 999, 9999)
    encontrados <- intersect(vec, suspeitos)
    if (length(encontrados) > 0) {
      contagem <- sum(vec %in% encontrados, na.rm = TRUE)
      # Pega o primeiro suspeito pra nomear a coluna de forma resumida
      nome_col <- sprintf("Detectado (%s)", encontrados[1])
      return(list(nome = nome_col, count = contagem))
    }
  } else {
    suspeitos <- c("", " ", "?", "-", "NA", "N/A", "null", "NULL")
    vec_char <- as.character(vec)
    encontrados <- intersect(trimws(vec_char), suspeitos)
    if (length(encontrados) > 0) {
      contagem <- sum(trimws(vec_char) %in% encontrados, na.rm = TRUE)
      nome_col <- sprintf("Detectado ('%s')", encontrados[1])
      if (encontrados[1] == "" || encontrados[1] == " ") nome_col <- "Detectado (vazio)"
      return(list(nome = nome_col, count = contagem))
    }
  }
  return(NULL)
}

.calc_out_resumo <- function(vec) {
  vec_c <- vec[!is.na(vec)]
  n <- length(vec_c)
  if (n < 4) return(c("-", "-"))
  q1 <- quantile(vec_c, 0.25)
  q3 <- quantile(vec_c, 0.75)
  iqr <- q3 - q1
  li <- q1 - 1.5 * iqr
  ls <- q3 + 1.5 * iqr
  n_out <- sum(vec_c < li | vec_c > ls)
  if (n_out == 0) return(c("-", "-"))
  pct_out <- sprintf("%.1f%%", (n_out / n) * 100)
  c(as.character(n_out), pct_out)
}

.resumo_variaveis <- function(dados, tipo = c("numerica", "categorica"), 
                              completo = TRUE, ausentes = FALSE) {
  tipo <- match.arg(tipo)
  
  if (tipo == "numerica") {
    cols <- names(dados)[sapply(dados, is.numeric)]
  } else {
    cols <- names(dados)[sapply(dados, function(x) is.factor(x) || is.character(x) || is.logical(x))]
  }
  
  if (length(cols) == 0) return(NULL)
  
  n_total <- nrow(dados)
  df_resumo <- data.frame("Variável" = cols, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!ausentes) {
    df_resumo$Tipo <- sapply(dados[cols], .tipo_traduzido)
  }
  
  if (completo || ausentes) {
    df_resumo[["Total (N)"]] <- n_total
    
    nas <- sapply(dados[cols], function(x) sum(is.na(x)))
    df_resumo$NA_oficial <- nas
    names(df_resumo)[names(df_resumo) == "NA_oficial"] <- "NA"
    df_resumo[["NA (%)"]] <- sprintf("%.1f%%", (nas / n_total) * 100)
    
    # Detecção de suspeitos ocultos
    suspeitos_detectados <- lapply(dados[cols], .detectar_suspeitos)
    nomes_susp <- unique(unlist(lapply(suspeitos_detectados, function(x) x$nome)))
    if (length(nomes_susp) > 0) {
      for (nm in nomes_susp) {
        df_resumo[[nm]] <- sapply(suspeitos_detectados, function(x) {
          if (!is.null(x) && x$nome == nm) as.character(x$count) else "-"
        })
      }
    }
  }

  if (completo) {
    # Valores Únicos para todos os tipos (antes de Níveis/Outliers)
    df_resumo[["Únicos"]] <- sapply(dados[cols], function(x) length(unique(na.omit(x))))
  }

  if (!ausentes && tipo == "categorica") {
    df_resumo[["Níveis (Levels)"]] <- sapply(dados[cols], .truncar_niveis)
  }
  
  if (completo && tipo == "numerica") {
    outs <- lapply(dados[cols], .calc_out_resumo)
    df_resumo[["Outliers"]] <- sapply(outs, `[`, 1)
    df_resumo[["Out (%)"]]  <- sapply(outs, `[`, 2)
  }
  
  return(df_resumo)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Diagnóstico Completo do Banco de Dados
#' @description Gera um painel completo mesclando informações de estrutura, 
#' valores ausentes, outliers e valores únicos, alertando para duplicatas.
#' @param dados Data frame ou vetor para análise.
#' @return Retorna o objeto original invisivelmente.
#' @export
diagnostico <- function(dados) {
  is_vetor <- FALSE
  obj_original <- dados
  
  if (!is.data.frame(dados)) {
    var_nome <- sub(".*\\$", "", deparse(substitute(dados)))
    if (var_nome == "" || var_nome == "dados") var_nome <- "Variavel"
    dados <- data.frame(dados, stringsAsFactors = FALSE)
    names(dados) <- var_nome
    is_vetor <- TRUE
  }
  
  n_linhas <- nrow(dados)
  n_cols   <- ncol(dados)
  
  if (is_vetor) {
    .print_titulo(sprintf("DIAGNÓSTICO DE VARIÁVEL — %s", names(dados)[1]))
  } else {
    .print_titulo("DIAGNÓSTICO DO BANCO DE DADOS")
    cat(sprintf("  Total de observações (linhas): %d\n", n_linhas))
    cat(sprintf("  Total de variáveis (colunas):  %d\n\n", n_cols))
  }
  
  # Numéricas
  df_num <- .resumo_variaveis(dados, "numerica", completo = TRUE)
  if (!is.null(df_num)) {
    pct_num <- (nrow(df_num) / n_cols) * 100
    titulo_sec <- if (is_vetor) "RESUMO DA VARIÁVEL" else sprintf("NUMÉRICAS (%.1f%%)", pct_num)
    .print_topico(titulo_sec)
    .print_tabela_estatR(df_num)
  }
  
  # Categóricas
  df_cat <- .resumo_variaveis(dados, "categorica", completo = TRUE)
  if (!is.null(df_cat)) {
    pct_cat <- (nrow(df_cat) / n_cols) * 100
    titulo_sec <- if (is_vetor) "RESUMO DA VARIÁVEL" else sprintf("CATEGÓRICAS E FATORES (%.1f%%)", pct_cat)
    .print_topico(titulo_sec)
    
    align_cat <- rep("center", ncol(df_cat))
    align_cat[1] <- "left"
    if ("Níveis (Levels)" %in% names(df_cat)) align_cat[which(names(df_cat) == "Níveis (Levels)")] <- "left"
    
    .print_tabela_estatR(df_cat, align = align_cat)
  }
  
  # Duplicatas (Somente se for Data Frame com > 1 coluna)
  if (!is_vetor && n_cols > 1) {
    alertas <- character(0)
    
    # Busca 1: Colunas exatas idênticas
    for (i in 1:(n_cols - 1)) {
      for (j in (i + 1):n_cols) {
        if (identical(dados[[i]], dados[[j]])) {
          alertas <- c(alertas, sprintf("  [!] : '%s' e '%s' possuem dados exatos idênticos.", 
                                        names(dados)[i], names(dados)[j]))
        }
      }
    }
    
    # Busca 2: Correlação > 0.99 para numéricas
    num_cols <- names(dados)[sapply(dados, is.numeric)]
    if (length(num_cols) > 1) {
      for (i in 1:(length(num_cols) - 1)) {
        for (j in (i + 1):length(num_cols)) {
          v1 <- dados[[num_cols[i]]]
          v2 <- dados[[num_cols[j]]]
          # Pula se for identical pois já pegou na busca 1
          if (!identical(v1, v2)) {
            comp <- complete.cases(v1, v2)
            if (sum(comp) > 2) {
              r <- cor(v1[comp], v2[comp])
              if (!is.na(r) && abs(r) > 0.99) {
                alertas <- c(alertas, sprintf("  [!] : '%s' e '%s' altamente correlacionadas (r = %.3f).", 
                                              num_cols[i], num_cols[j], r))
              }
            }
          }
        }
      }
    }
    
    if (length(alertas) > 0) {
      .print_topico("VARIÁVEIS DUPLICADAS")
      cat(paste(alertas, collapse = "\n"), "\n\n")
    }
  }
  
  .print_rodape()
  invisible(obj_original)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Estrutura do Banco de Dados
#' @description Retorna apenas os tipos de dados de cada variável.
#' @param dados Data frame ou vetor para análise.
#' @return Retorna o objeto original invisivelmente.
#' @export
estrutura <- function(dados) {
  is_vetor <- FALSE
  obj_original <- dados
  
  if (!is.data.frame(dados)) {
    var_nome <- sub(".*\\$", "", deparse(substitute(dados)))
    if (var_nome == "" || var_nome == "dados") var_nome <- "Variavel"
    dados <- data.frame(dados, stringsAsFactors = FALSE)
    names(dados) <- var_nome
    is_vetor <- TRUE
  }
  
  n_linhas <- nrow(dados)
  n_cols   <- ncol(dados)
  
  if (is_vetor) {
    .print_titulo(sprintf("ESTRUTURA DE VARIÁVEL — %s", names(dados)[1]))
  } else {
    .print_titulo("ESTRUTURA DO BANCO DE DADOS")
    cat(sprintf("  Total de observações (linhas): %d\n", n_linhas))
    cat(sprintf("  Total de variáveis (colunas):  %d\n\n", n_cols))
  }
  
  df_num <- .resumo_variaveis(dados, "numerica", completo = FALSE)
  if (!is.null(df_num)) {
    titulo_sec <- if (is_vetor) "NUMÉRICA" else sprintf("NUMÉRICAS (%.1f%%)", (nrow(df_num) / n_cols) * 100)
    .print_topico(titulo_sec)
    .print_tabela_estatR(df_num)
  }
  
  df_cat <- .resumo_variaveis(dados, "categorica", completo = FALSE)
  if (!is.null(df_cat)) {
    titulo_sec <- if (is_vetor) "CATEGÓRICA / FATOR" else sprintf("CATEGÓRICAS E FATORES (%.1f%%)", (nrow(df_cat) / n_cols) * 100)
    .print_topico(titulo_sec)
    
    align_cat <- rep("center", ncol(df_cat))
    align_cat[1] <- "left"
    if ("Níveis (Levels)" %in% names(df_cat)) align_cat[which(names(df_cat) == "Níveis (Levels)")] <- "left"
    
    .print_tabela_estatR(df_cat, align = align_cat)
  }
  
  .print_rodape()
  invisible(obj_original)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Análise de Valores Ausentes
#' @description Inspeciona procurando por NAs oficiais e valores
#' que indicam respostas vazias ou erros ocultos (como -99 ou textos vazios).
#' @param dados Data frame ou vetor para análise.
#' @return Retorna o objeto original invisivelmente.
#' @export
valores_ausentes <- function(dados) {
  is_vetor <- FALSE
  obj_original <- dados
  
  if (!is.data.frame(dados)) {
    var_nome <- sub(".*\\$", "", deparse(substitute(dados)))
    if (var_nome == "" || var_nome == "dados") var_nome <- "Variavel"
    dados <- data.frame(dados, stringsAsFactors = FALSE)
    names(dados) <- var_nome
    is_vetor <- TRUE
  }
  
  n_cols <- ncol(dados)
  
  if (is_vetor) {
    .print_titulo(sprintf("VALORES AUSENTES — %s", names(dados)[1]))
  } else {
    .print_titulo("ANÁLISE DE VALORES AUSENTES")
  }
  
  df_num <- .resumo_variaveis(dados, "numerica", completo = FALSE, ausentes = TRUE)
  if (!is.null(df_num)) {
    titulo_sec <- if (is_vetor) "NUMÉRICA" else sprintf("NUMÉRICAS (%.1f%%)", (nrow(df_num) / n_cols) * 100)
    .print_topico(titulo_sec)
    .print_tabela_estatR(df_num)
  }
  
  df_cat <- .resumo_variaveis(dados, "categorica", completo = FALSE, ausentes = TRUE)
  if (!is.null(df_cat)) {
    titulo_sec <- if (is_vetor) "CATEGÓRICA / FATOR" else sprintf("CATEGÓRICAS E FATORES (%.1f%%)", (nrow(df_cat) / n_cols) * 100)
    .print_topico(titulo_sec)
    .print_tabela_estatR(df_cat)
  }
  
  .print_rodape()
  invisible(obj_original)
}
