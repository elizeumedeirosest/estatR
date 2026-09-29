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
    if (tipo == "categorica") {
      df_resumo[["Níveis (Levels)"]] <- sapply(dados[cols], .truncar_niveis)
    }
  }
  
  if (completo || ausentes) {
    df_resumo[["Total (N)"]] <- n_total
    
    nas <- sapply(dados[cols], function(x) sum(is.na(x)))
    df_resumo$NA_oficial <- nas
    names(df_resumo)[names(df_resumo) == "NA_oficial"] <- "NA"
    
    df_resumo[["NA (%)"]] <- sprintf("%.1f%%", (nas / n_total) * 100)
    
    # Detecção de suspeitos
    suspeitos_detectados <- lapply(dados[cols], .detectar_suspeitos)
    
    # Encontrar todos os nomes de colunas de suspeitos encontrados
    nomes_susp <- unique(unlist(lapply(suspeitos_detectados, function(x) x$nome)))
    
    if (length(nomes_susp) > 0) {
      for (nm in nomes_susp) {
        df_resumo[[nm]] <- sapply(suspeitos_detectados, function(x) {
          if (!is.null(x) && x$nome == nm) as.character(x$count) else "-"
        })
      }
    }
  }
  
  return(df_resumo)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Diagnóstico Completo do Banco de Dados
#' @description Gera um painel completo mesclando informações de estrutura e 
#' valores ausentes por tipo de variável, incluindo detecção automática de 
#' NAs ocultos (como -99 ou vazio).
#' @param dados Data frame para análise.
#' @return Retorna o data frame invisivelmente.
#' @export
diagnostico <- function(dados) {
  if (!is.data.frame(dados)) stop("'dados' deve ser um data frame.")
  
  n_linhas <- nrow(dados)
  n_cols   <- ncol(dados)
  
  .print_titulo("DIAGNÓSTICO DO BANCO DE DADOS")
  
  cat(sprintf("  Total de observações (linhas): %d\n", n_linhas))
  cat(sprintf("  Total de variáveis (colunas):  %d\n\n", n_cols))
  
  # Numéricas
  df_num <- .resumo_variaveis(dados, "numerica", completo = TRUE)
  if (!is.null(df_num)) {
    pct_num <- (nrow(df_num) / n_cols) * 100
    .print_topico(sprintf("NUMÉRICAS (%.1f%%)", pct_num))
    .print_tabela_estatR(df_num)
  }
  
  # Categóricas
  df_cat <- .resumo_variaveis(dados, "categorica", completo = TRUE)
  if (!is.null(df_cat)) {
    pct_cat <- (nrow(df_cat) / n_cols) * 100
    .print_topico(sprintf("CATEGÓRICAS E FATORES (%.1f%%)", pct_cat))
    
    # Alinhamento especifico para a tabela categorica (Niveis alinhado a esquerda)
    align_cat <- rep("center", ncol(df_cat))
    align_cat[1] <- "left"
    if ("Níveis (Levels)" %in% names(df_cat)) {
      align_cat[which(names(df_cat) == "Níveis (Levels)")] <- "left"
    }
    
    .print_tabela_estatR(df_cat, align = align_cat)
  }
  
  .print_rodape()
  invisible(dados)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Estrutura do Banco de Dados
#' @description Retorna apenas os tipos de dados de cada variável.
#' @param dados Data frame para análise.
#' @return Retorna o data frame invisivelmente.
#' @export
estrutura <- function(dados) {
  if (!is.data.frame(dados)) stop("'dados' deve ser um data frame.")
  
  n_linhas <- nrow(dados)
  n_cols   <- ncol(dados)
  
  .print_titulo("ESTRUTURA DO BANCO DE DADOS")
  
  cat(sprintf("  Total de observações (linhas): %d\n", n_linhas))
  cat(sprintf("  Total de variáveis (colunas):  %d\n\n", n_cols))
  
  df_num <- .resumo_variaveis(dados, "numerica", completo = FALSE)
  if (!is.null(df_num)) {
    .print_topico(sprintf("NUMÉRICAS (%.1f%%)", (nrow(df_num) / n_cols) * 100))
    .print_tabela_estatR(df_num)
  }
  
  df_cat <- .resumo_variaveis(dados, "categorica", completo = FALSE)
  if (!is.null(df_cat)) {
    .print_topico(sprintf("CATEGÓRICAS E FATORES (%.1f%%)", (nrow(df_cat) / n_cols) * 100))
    
    align_cat <- rep("center", ncol(df_cat))
    align_cat[1] <- "left"
    if ("Níveis (Levels)" %in% names(df_cat)) align_cat[which(names(df_cat) == "Níveis (Levels)")] <- "left"
    
    .print_tabela_estatR(df_cat, align = align_cat)
  }
  
  .print_rodape()
  invisible(dados)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Análise de Valores Ausentes
#' @description Inspeciona o banco procurando por NAs oficiais e valores
#' que indicam respostas vazias ou erros ocultos (como -99 ou textos vazios).
#' @param dados Data frame para análise.
#' @return Retorna o data frame invisivelmente.
#' @export
valores_ausentes <- function(dados) {
  if (!is.data.frame(dados)) stop("'dados' deve ser um data frame.")
  
  n_cols <- ncol(dados)
  
  .print_titulo("ANÁLISE DE VALORES AUSENTES")
  
  df_num <- .resumo_variaveis(dados, "numerica", completo = FALSE, ausentes = TRUE)
  if (!is.null(df_num)) {
    .print_topico(sprintf("NUMÉRICAS (%.1f%%)", (nrow(df_num) / n_cols) * 100))
    .print_tabela_estatR(df_num)
  }
  
  df_cat <- .resumo_variaveis(dados, "categorica", completo = FALSE, ausentes = TRUE)
  if (!is.null(df_cat)) {
    .print_topico(sprintf("CATEGÓRICAS E FATORES (%.1f%%)", (nrow(df_cat) / n_cols) * 100))
    .print_tabela_estatR(df_cat)
  }
  
  .print_rodape()
  invisible(dados)
}
